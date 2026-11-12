"""Unit tests of the python/hscp package. No ROOT needed:  python3 -m unittest discover -s tests -v"""

import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "python"))
from hscp import config, driver, selector  # noqa: E402


class ConfigTests(unittest.TestCase):
    DOC = {
        "prod_dir": "/p/",
        "binning": {"ptcut": 70, "etabins": 96, "ihbins": 400, "pbins": 1000, "massbins": 400, "fpixbins": 20},
        "samples": {
            "A": {"list": "L/a.txt", "type": "data", "dataset": "ds1"},
            "B": {"file": "/abs/b.root", "type": "mc", "dataset": "ds1", "binning": {"etabins": 12}},
            "C": {"file": "root://host//store/c.root", "type": "signal", "dataset": "ds2", "mass": 700},
        },
    }

    def _cfg(self, doc=None, **kw):
        os.environ.pop("HSCP_PROD_DIR", None)
        return config.Config(doc or self.DOC, "test.json", **kw)

    def test_repository_file_is_valid(self):
        cfg = config.load(prod_dir="/x/")
        self.assertIn("TestMET2024", cfg.names)
        self.assertEqual(cfg.binning("TestMET2024"), config.Binning(70.0, 96, 400, 1000, 400, 20))
        self.assertEqual(cfg.source("TestMET2024"), ("list", "/x/JetMET2024/testJetMET.txt"))
        self.assertEqual(cfg.type("Stau_Run3_MET_247"), "signal")
        self.assertIn("JetMET2024", cfg.datasets())
        self.assertEqual(cfg.dataset("JetMET2024C"), "JetMET2024")
        self.assertEqual(len(cfg.select(datasets=["JetMET2024"])), 7)

    def test_prod_dir_precedence(self):
        self.assertEqual(self._cfg().prod_dir, "/p/")
        self.assertEqual(self._cfg(prod_dir="/arg").prod_dir, "/arg/")
        os.environ["HSCP_PROD_DIR"] = "/env"
        try:
            self.assertEqual(config.Config(self.DOC, "t").prod_dir, "/env/")
            self.assertEqual(config.Config(self.DOC, "t", prod_dir="/arg/").prod_dir, "/arg/")
        finally:
            del os.environ["HSCP_PROD_DIR"]

    def test_paths_and_overrides(self):
        cfg = self._cfg()
        self.assertEqual(cfg.source("A"), ("list", "/p/L/a.txt"))
        self.assertEqual(cfg.source("B"), ("file", "/abs/b.root"))          # absolute: no prod_dir
        self.assertEqual(cfg.source("C"), ("file", "root://host//store/c.root"))
        self.assertEqual(cfg.binning("A"), cfg.default_binning)
        self.assertEqual(cfg.binning("B").etabins, 12)
        self.assertEqual(cfg.binning("B").ihbins, 400)
        self.assertEqual(cfg.files("B"), ["/abs/b.root"])

    def test_select_by_pattern_and_dataset(self):
        cfg = self._cfg()
        self.assertEqual(cfg.select(["A", "C"]), ["A", "C"])
        self.assertEqual(cfg.select(["[AB]"]), ["A", "B"])
        self.assertEqual(cfg.select(datasets=["ds1"]), ["A", "B"])
        self.assertEqual(cfg.select(["C"], ["ds1"]), ["A", "B", "C"])      # file order, no duplicates
        with self.assertRaises(config.ConfigError):
            cfg.select(["Nothing*"])
        with self.assertRaises(config.ConfigError):
            cfg.select(datasets=["ds9"])

    def test_file_list_reading(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "L"))
            with open(os.path.join(d, "L", "a.txt"), "w") as f:
                f.write("/a/1.root\n\n# comment\n/a/2.root\n")
            cfg = self._cfg(prod_dir=d)
            self.assertEqual(cfg.files("A"), ["/a/1.root", "/a/2.root"])
            with open(os.path.join(d, "L", "a.txt"), "w") as f:
                f.write("\n")
            with self.assertRaises(config.ConfigError):
                cfg.files("A")
        with self.assertRaises(config.ConfigError):
            self._cfg(prod_dir="/nowhere").files("A")

    def test_validation(self):
        import copy
        for mutate in (
            lambda d: d["samples"]["A"].pop("type"),
            lambda d: d["samples"]["A"].update(type="mcc"),
            lambda d: d["samples"]["A"].update(file="x.root"),             # list and file
            lambda d: d["samples"]["A"].pop("list"),
            lambda d: d["samples"]["A"].update(unknown=1),
            lambda d: d["samples"]["B"]["binning"].update(etabins="many"),
            lambda d: d["binning"].pop("ptcut"),
            lambda d: d.pop("samples"),
        ):
            doc = copy.deepcopy(self.DOC)
            mutate(doc)
            with self.assertRaises(config.ConfigError):
                self._cfg(doc)

    def test_load_reports_bad_json(self):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
            f.write("{not json")
        try:
            with self.assertRaises(config.ConfigError):
                config.load(f.name)
        finally:
            os.remove(f.name)

    def test_full_path(self):
        self.assertEqual(config.full_path("A/b.txt", "/p/"), "/p/A/b.txt")
        self.assertEqual(config.full_path("/data/x.root", "/p/"), "/data/x.root")
        self.assertEqual(config.full_path("root://eos//x.root", "/p/"), "root://eos//x.root")


class SelectorTests(unittest.TestCase):
    def _write(self, doc):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
            if isinstance(doc, str):
                f.write(doc)
            else:
                json.dump(doc, f)
        self.addCleanup(os.remove, f.name)
        return f.name

    def test_repository_selection_file(self):
        sels = selector.read_selections()
        self.assertEqual(len(sels), 12)
        self.assertTrue(all(selector.IDENT.match(s.label) for s in sels))
        for s in sels:                                   # the common cuts come first, then the specific ones
            self.assertEqual(s.cuts[0], "Flag_allMETFilters[0] == true")
            self.assertEqual(s.cuts[-1][:12], "ptOverptErr[")

    def test_common_cuts_object_form_and_disabled(self):
        doc = {"common": ["Pt[i] > 50"],
               "selections": {"a": ["Eta[i] < 1"],
                              "b": {"cuts": ["Eta[i] < 2"], "comment": "object form"},
                              "c": {"cuts": ["Eta[i] < 3"], "enabled": False},
                              "_note": "ignored"}}
        self.assertEqual(selector.read_selections(self._write(doc)),
                         [("a", ["Pt[i] > 50", "Eta[i] < 1"]), ("b", ["Pt[i] > 50", "Eta[i] < 2"])])
        self.assertEqual(selector.expression(["Pt[i] > 50", "Eta[i] < 1"]), "(Pt[i] > 50) && (Eta[i] < 1)")

    def test_rejects_bad_files(self):
        bad = ('{"selections": {"a": ["x"]',                                   # not JSON
               {"selections": {"bad-label": ["x"]}},
               {"selections": {"a": []}},
               {"selections": {"a": ["x", 3]}},
               {"selections": {"a": {"enabled": True}}},                      # no cuts
               {"selections": {"a": {"cuts": ["x"], "typo": 1}}},
               {"selections": {"a": ["x"]}, "commmon": ["y"]},
               {"selections": {"a": {"cuts": ["x"], "enabled": False}}},      # nothing active
               {"common": ["x"]},                                             # no selections
               ["a"])
        for doc in bad:
            with self.assertRaises(selector.SelectionError, msg=repr(doc)):
                selector.read_selections(self._write(doc))

    def test_generate_inserts_each_selection(self):
        sels = [("selA", ["Pt[i] > 50.0"]), ("selB", ["Eta[i] < 1.0", "Pt[i] > 60"])]
        with tempfile.TemporaryDirectory() as d:
            h, c = selector.generate(sels, d)
            with open(h) as fh, open(c) as fc:
                hdr, src = fh.read(), fc.read()
            for label in ("selA", "selB"):
                self.assertIn(f"bool PassHSCPpresel_{label}(int hscpIndex);", hdr)
                self.assertIn(f"bool HSCPSelector::PassHSCPpresel_{label}(int i)", src)
                self.assertIn(f'selLabels_.push_back("{label}");', src)
            self.assertIn("i>=(int)Pt.GetSize()", src)
            self.assertIn("return (Pt[i] > 50.0);", src)
            self.assertIn("return (Eta[i] < 1.0) && (Pt[i] > 60);", src)

    def test_write_then_read_round_trip(self):
        sels = [selector.Selection("a", ["x"]), selector.Selection("b", ["y", "z"])]
        path = self._write("")
        selector.write_selections(sels, path)
        self.assertEqual(selector.read_selections(path), sels)


class DriverTests(unittest.TestCase):
    def test_macro_call_argument_order(self):
        b = config.Binning(70.0, 96, 400, 1000, 400, 20)
        self.assertEqual(driver.macro_call("D", "", "out.root", True, b),
                         'macro.cc("D","","out.root",true,70,96,400,1000,400,20)')
        self.assertEqual(driver.macro_call("D", "/l.txt", "o.root", False, config.Binning(72.5, 1, 2, 3, 4, 5)),
                         'macro.cc("D","/l.txt","o.root",false,72.5,1,2,3,4,5)')
        self.assertEqual(driver.root_command("D", "", "o.root", True, b)[:4], ["root", "-l", "-b", "-q"])

    def test_build_loads_libtools_before_compiling(self):
        cmd = driver.build_command()
        self.assertIn('gSystem->Load("../libTools.so")', cmd)
        self.assertLess(cmd.index('gSystem->Load("../libTools.so")'), [i for i, c in enumerate(cmd) if "CompileMacro" in c][0])

    def test_runjob_uses_the_same_argument_order(self):
        """condor/runJob.sh encodes the macro.cc call by hand: keep it aligned with driver.macro_call."""
        with open(os.path.join(os.path.dirname(HERE), "condor", "runJob.sh")) as f:
            script = f.read()
        self.assertIn('macro.cc(\\"$SAMPLE\\",\\"$CHUNK_ABS\\",\\"$OUT_ABS\\",$NOSEL,$PTCUT,$ETABINS,$IHBINS,$PBINS,$MASSBINS,$FPIXBINS)', script)


if __name__ == "__main__":
    unittest.main()
