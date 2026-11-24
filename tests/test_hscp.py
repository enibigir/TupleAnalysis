"""Unit tests of the python/hscp package. No ROOT needed:  python3 -m unittest discover -s tests -v"""

import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "python"))
from hscp import config, driver, plots, proxy, selector, stack, storage, weight  # noqa: E402


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
        self.assertIn("TestJetMET2024C", cfg.names)
        self.assertEqual(cfg.binning("TestJetMET2024C"), config.Binning(70.0, 96, 400, 1000, 400, 20))
        self.assertEqual(cfg.source("TestJetMET2024C"), ("list", "/x/input/TestJetMET2024C.txt"))
        for name in cfg.names:          # every simulated sample can be weighted (python/Weight.py)
            if cfg.type(name) != "data":
                self.assertIsNotNone(cfg.xsec(name), name)
        self.assertEqual(cfg.type("Gluino_Run3_MET_madgraph_1600"), "signal")
        self.assertIn("JetMET2024", cfg.datasets())
        self.assertEqual(cfg.dataset("JetMET2024C"), "JetMET2024")
        self.assertEqual(len(cfg.select(datasets=["JetMET2024"])), 1)

    def test_empty_prod_dir_is_the_repository(self):
        doc = dict(self.DOC, prod_dir="")
        self.assertEqual(self._cfg(doc).prod_dir, config.REPO + "/")
        self.assertEqual(self._cfg(doc, prod_dir="/arg").prod_dir, "/arg/")
        self.assertEqual(self._cfg(doc).source("A"), ("list", config.REPO + "/L/a.txt"))
        missing = {k: v for k, v in self.DOC.items() if k != "prod_dir"}
        with self.assertRaises(config.ConfigError):
            self._cfg(missing)

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

    def test_describe_exit(self):
        self.assertEqual(driver.describe_exit(0), "exit code 0")
        self.assertIn("no output file", driver.describe_exit(driver.NO_OUTPUT))
        self.assertIn("signal 9 (SIGKILL)", driver.describe_exit(-9))
        self.assertIn("memory", driver.describe_exit(-9))
        self.assertIn("segmentation", driver.describe_exit(-11))
        self.assertIn("processing failed", driver.describe_exit(3))
        self.assertIn("signal 9 (SIGKILL) in a wrapper", driver.describe_exit(137))
        # the root launcher returns the signal of root.exe as its exit code
        self.assertIn("root.exe was killed by signal 9 (SIGKILL)", driver.describe_exit(9))
        self.assertIn("memory", driver.describe_exit(9))
        self.assertIn("segmentation", driver.describe_exit(11))
        self.assertEqual(driver.describe_exit(1), "exit code 1")

    def test_build_loads_libtools_before_compiling(self):
        cmd = driver.build_command()
        self.assertIn('gSystem->Load("../libTools.so")', cmd)
        self.assertLess(cmd.index('gSystem->Load("../libTools.so")'), [i for i, c in enumerate(cmd) if "CompileMacro" in c][0])

    def test_runjob_uses_the_same_argument_order(self):
        """condor/runJob.sh encodes the macro.cc call by hand: keep it aligned with driver.macro_call."""
        with open(os.path.join(os.path.dirname(HERE), "condor", "runJob.sh")) as f:
            script = f.read()
        self.assertIn('macro.cc(\\"$SAMPLE\\",\\"$CHUNK_ABS\\",\\"$OUT_ABS\\",$NOSEL,$PTCUT,$ETABINS,$IHBINS,$PBINS,$MASSBINS,$FPIXBINS)', script)


class ProxyTests(unittest.TestCase):
    LOCAL = ["/data/a.root", "file:///data/b.root"]
    REMOTE = LOCAL + ["davs://host//store/c.root"]

    def setUp(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        f.close()
        self.addCleanup(os.remove, f.name)
        self.file = f.name
        self._env = os.environ.pop("X509_USER_PROXY", None)
        self.addCleanup(lambda: os.environ.__setitem__("X509_USER_PROXY", self._env) if self._env
                        else os.environ.pop("X509_USER_PROXY", None))

    def test_remote_detection(self):
        self.assertFalse(any(proxy.is_remote(f) for f in self.LOCAL))
        self.assertTrue(proxy.is_remote("root://eos//x.root"))

    def test_auto(self):
        long = lambda p: 100 * 3600   # noqa: E731
        self.assertEqual(proxy.resolve("auto", self.LOCAL, 3600, long), (None, None))
        os.environ["X509_USER_PROXY"] = self.file
        self.assertEqual(proxy.resolve("auto", self.REMOTE, 3600, long), (self.file, 100 * 3600))
        os.environ["X509_USER_PROXY"] = self.file + ".missing"
        with self.assertRaisesRegex(proxy.ProxyError, "voms-proxy-init"):
            proxy.resolve("auto", self.REMOTE, 3600, long)

    def test_none_and_explicit_path(self):
        self.assertEqual(proxy.resolve("none", self.REMOTE, 3600), (None, None))
        self.assertEqual(proxy.resolve(self.file, self.LOCAL, 3600, lambda p: None), (self.file, None))
        with self.assertRaises(proxy.ProxyError):
            proxy.resolve(self.file + ".missing", self.LOCAL, 3600, lambda p: None)

    def test_lifetime(self):
        with self.assertRaisesRegex(proxy.ProxyError, "1 h 30 min left"):
            proxy.resolve(self.file, [], 3600, lambda p: 5400)          # runtime + 1 h = 2 h needed
        self.assertEqual(proxy.resolve(self.file, [], 3600, lambda p: 7200)[1], 7200)

    def test_sub_line(self):
        self.assertEqual(proxy.sub_line(None), "")
        self.assertEqual(proxy.sub_line("/tmp/x509up_u1"), "x509userproxy         = /tmp/x509up_u1\n")
        self.assertIn("not checked", proxy.describe("/p", None))


class WeightTests(unittest.TestCase):
    DOC = {
        "prod_dir": "/p/", "luminosity": 1000,
        "binning": ConfigTests.DOC["binning"],
        "samples": {
            "Q1": {"file": "q1.root", "type": "mc", "dataset": "QCD", "xsec": 20},
            "Q2": {"file": "q2.root", "type": "mc", "dataset": "QCD", "xsec": 3, "sum_weights": 600},
            "Q3": {"file": "q3.root", "type": "mc", "dataset": "QCD"},
            "D1": {"file": "d1.root", "type": "data", "dataset": "Data"},
            "S1": {"file": "s1.root", "type": "signal", "dataset": "Sig", "xsec": 0.5},
        },
    }

    def setUp(self):
        os.environ.pop("HSCP_PROD_DIR", None)
        self.cfg = config.Config(self.DOC, "test.json")
        self.read = []

    def _sumw(self, s):
        self.read.append(s)
        return 400.0

    def test_config_values(self):
        self.assertEqual(self.cfg.luminosity, 1000.0)
        self.assertEqual(self.cfg.xsec("Q1"), 20.0)
        self.assertIsNone(self.cfg.xsec("Q3"))
        self.assertEqual(self.cfg.sum_weights("Q2"), 600.0)
        import copy
        for bad in ({"xsec": -1}, {"xsec": "1"}, {"sum_weights": 0}, {"xsec": True}):
            doc = copy.deepcopy(self.DOC)
            doc["samples"]["Q1"].update(bad)
            with self.assertRaises(config.ConfigError, msg=bad):
                config.Config(doc, "t")
        doc = dict(self.DOC, luminosity=0)
        with self.assertRaises(config.ConfigError):
            config.Config(doc, "t")

    def test_factors(self):
        entries = weight.plan(self.cfg, ["Q1", "Q2"], 1000.0, self._sumw)
        self.assertEqual([e.factor for e in entries], [20 * 1000 / 400, 3 * 1000 / 600])
        self.assertEqual([e.source for e in entries], [weight.NORM_HISTO, "datasets.json"])
        self.assertEqual(self.read, ["Q1"])                       # sum_weights of datasets.json wins
        self.assertEqual(weight.plan(self.cfg, ["S1"], 2.0, self._sumw)[0].factor, 0.5 * 2 / 400)
        self.assertIn("luminosity: 1000 pb^-1", weight.table(entries, 1000.0))

    def test_data_is_added_only(self):
        entries = weight.plan(self.cfg, ["D1"], None, self._sumw)
        self.assertEqual([(e.type, e.factor) for e in entries], [("data", 1.0)])
        self.assertEqual(self.read, [])

    def test_errors(self):
        with self.assertRaisesRegex(weight.WeightError, "data and simulation"):
            weight.plan(self.cfg, ["D1", "Q1"], 1000.0, self._sumw)
        with self.assertRaisesRegex(weight.WeightError, "no 'xsec'.*Q3"):
            weight.plan(self.cfg, ["Q1", "Q3"], 1000.0, self._sumw)
        with self.assertRaisesRegex(weight.WeightError, "no luminosity"):
            weight.plan(self.cfg, ["Q1"], None, self._sumw)
        with self.assertRaisesRegex(weight.WeightError, "Q1: sum of weights 0"):
            weight.plan(self.cfg, ["Q1"], 1000.0, lambda s: 0.0)

    def test_missing_inputs(self):
        with tempfile.TemporaryDirectory() as d:
            paths = {s: os.path.join(d, f"{s}.root") for s in ("Q1", "Q2")}
            open(paths["Q1"], "w").close()
            with self.assertRaisesRegex(weight.WeightError, "missing input.*--skip-missing"):
                weight.available(paths)
            self.assertEqual(weight.available(paths, skip_missing=True), (["Q1"], [paths["Q2"]]))
            os.remove(paths["Q1"])
            with self.assertRaisesRegex(weight.WeightError, "none of the input files"):
                weight.available(paths, skip_missing=True)

    def test_default_output(self):
        self.assertEqual(weight.default_output(self.cfg, ["Q1", "Q2"], "m"), os.path.join("m", "QCD_weighted.root"))
        self.assertIsNone(weight.default_output(self.cfg, ["Q1", "S1"], "m"))


class StorageTests(unittest.TestCase):
    def test_resolve(self):
        os.environ["HSCP_TEST_USER"] = "me"
        try:
            self.assertEqual(storage.resolve("root://cmseos.fnal.gov//store/user/$HSCP_TEST_USER/hscp/", "qcd"),
                             "root://cmseos.fnal.gov//store/user/me/hscp/qcd")
        finally:
            del os.environ["HSCP_TEST_USER"]
        self.assertEqual(storage.resolve("/data/x/", "t"), "/data/x/t")
        self.assertEqual(storage.resolve("rel", "t"), os.path.join(os.getcwd(), "rel", "t"))
        with self.assertRaisesRegex(storage.StorageError, "double slash"):
            storage.resolve("root://cmseos.fnal.gov/store/user/me", "t")
        with self.assertRaisesRegex(storage.StorageError, "undefined variable"):
            storage.resolve("/data/$HSCP_UNDEFINED_VARIABLE", "t")

    def test_join_and_split(self):
        self.assertEqual(storage.join("root://h//store/a", "f.root"), "root://h//store/a/f.root")
        self.assertEqual(storage.join("/d", "f.root"), "/d/f.root")
        self.assertEqual(storage.split_root_url("root://cmseos.fnal.gov:1094//store/user/me/t"),
                         ("root://cmseos.fnal.gov:1094", "/store/user/me/t"))

    def test_parsers(self):
        xrdfs = ("-rw-r--r-- me us_cms 7062264 2026-11-17 10:00:00 /store/user/me/t/A_g0_c0.root\n"
                 "drwxr-xr-x me us_cms 0 2026-11-17 10:00:00 /store/user/me/t/sub\n"
                 "-r-- 2026-11-17 10:00:01 512 /store/user/me/t/A_g0_c1.root\n")
        self.assertEqual(storage.parse_xrdfs(xrdfs), {"A_g0_c0.root": 7062264, "A_g0_c1.root": 512})
        gfal = ("-rw-r--r--   1 0     0      7062264 Nov 17 10:00 A_g0_c0.root\n"
                "drwxr-xr-x   1 0     0            0 Nov 17 10:00 sub\n")
        self.assertEqual(storage.parse_gfal(gfal), {"A_g0_c0.root": 7062264})

    def test_location_file(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(storage.read_location(d), os.path.join(d, "out"))
            with open(os.path.join(d, storage.LOCATION_FILE), "w") as f:
                f.write("root://h//store/x/t\n")
            self.assertEqual(storage.read_location(d), "root://h//store/x/t")
            self.assertEqual(storage.listing(os.path.join(d, "nothing")), {})

    def test_resolve_follows_symlinks(self):
        with tempfile.TemporaryDirectory() as d:
            real = os.path.realpath(os.path.join(d, "real"))
            os.mkdir(real)
            os.symlink(real, os.path.join(d, "link"))
            self.assertEqual(storage.resolve(os.path.join(d, "link"), "t"), os.path.join(real, "t"))

    def test_lpc_shared_filesystems(self):
        lpc = "cmslpc301.fnal.gov"
        self.assertIsNone(storage.unshared_for_scheduler("/uscms_data/d3/me/TupleAnalysis", lpc))
        self.assertIsNone(storage.unshared_for_scheduler("/uscms/home/me/TupleAnalysis", lpc))
        for bad in ("/tmp/me/TupleAnalysis", "/storage/local/data1/me", "/eos/uscms/store/user/me", "/uscms_database"):
            self.assertRegex(storage.unshared_for_scheduler(bad, lpc) or "", "hold code 14", bad)
        self.assertIsNone(storage.unshared_for_scheduler("/tmp/me", "lxplus912.cern.ch"))   # LPC rule only on the LPC
        self.assertIsNone(storage.unshared_for_scheduler("/tmp/me", "node.example.org"))   # other sites: no check

    def test_lxplus_refuses_eos(self):
        for host in ("lxplus912.cern.ch", "lxplus9"):
            self.assertIsNone(storage.unshared_for_scheduler("/afs/cern.ch/user/m/me/TupleAnalysis", host))
            for bad in ("/eos/user/m/me/TupleAnalysis", "/eos/home-m/me/hscp/t", "/eos"):
                self.assertRegex(storage.unshared_for_scheduler(bad, host) or "", "EOS mount", bad)
        self.assertIsNone(storage.unshared_for_scheduler("/eosx/me", "lxplus912.cern.ch"))
        self.assertIsNone(storage.unshared_for_scheduler("/eos/user/m/me", "laptop.local"))


class SubmitTests(unittest.TestCase):
    """condor/submit.py run from a symlink to the repository writes the real paths in hscp.sub."""
    SAMPLE = "TestGluino_Run3_MET_madgraph_1600"

    def test_submit_file_uses_real_paths(self):
        import shutil
        import subprocess
        repo = os.path.realpath(config.REPO)
        tag = f"_unittest_symlink_{os.getpid()}"
        work = os.path.join(repo, "condor", "work", tag)
        with tempfile.TemporaryDirectory() as d:
            link = os.path.join(d, "checkout")
            os.symlink(repo, link)
            try:
                p = subprocess.run([sys.executable, os.path.join(link, "condor", "submit.py"), "--samples", self.SAMPLE,
                                    "--tag", tag, "--proxy", "none"], capture_output=True, text=True, cwd=d)
                self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
                self.assertIn("is a symlink", p.stdout)
                with open(os.path.join(work, "hscp.sub")) as f:
                    sub = f.read()
                with open(os.path.join(work, "jobs.txt")) as f:
                    jobs = f.read()
            finally:
                shutil.rmtree(work, ignore_errors=True)
        self.assertIn(f"initialdir            = {work}\n", sub)
        self.assertIn(f"executable            = {repo}/condor/runJob.sh\n", sub)
        self.assertNotIn(link, sub + jobs)
        self.assertIn(os.path.join(work, "chunks", f"{self.SAMPLE}_c0.txt"), jobs)
        self.assertIn(f"log                   = {work}/log/condor.log\n", sub)

    def test_log_dir(self):
        import shutil
        import subprocess
        repo = os.path.realpath(config.REPO)
        tag = f"_unittest_logdir_{os.getpid()}"
        work = os.path.join(repo, "condor", "work", tag)
        submit = os.path.join(repo, "condor", "submit.py")
        with tempfile.TemporaryDirectory() as d:
            logs = os.path.join(os.path.realpath(d), tag)
            try:
                p = subprocess.run([sys.executable, submit, "--samples", self.SAMPLE, "--tag", tag, "--proxy", "none",
                                    "--log-dir", d], capture_output=True, text=True)
                self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
                with open(os.path.join(work, "hscp.sub")) as f:
                    sub = f.read()
                with open(os.path.join(work, "log_dir.txt")) as f:
                    recorded = f.read().strip()
                self.assertTrue(os.path.isdir(logs))
                self.assertFalse(os.path.exists(os.path.join(work, "log")))
            finally:
                shutil.rmtree(work, ignore_errors=True)
            self.assertEqual(recorded, logs)
            self.assertIn(f"log                   = {logs}/condor.log\n", sub)
            self.assertIn(f"error                 = {logs}/$(sample)_g$(group)_c$(chunkid).err\n", sub)
            self.assertIn(f"logs:        {logs}", p.stdout)
        # a URL is refused before anything is written
        p = subprocess.run([sys.executable, submit, "--samples", self.SAMPLE, "--tag", tag, "--proxy", "none",
                            "--log-dir", "root://eosuser.cern.ch//eos/user/x/xx/logs"], capture_output=True, text=True)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("--log-dir must be a local directory", p.stderr)
        self.assertFalse(os.path.exists(work))



class MergeTests(unittest.TestCase):
    """condor/merge.py --check reads the outputs from --input-dir instead of output_dir.txt."""

    def test_input_dir(self):
        import subprocess
        merge = os.path.join(os.path.realpath(config.REPO), "condor", "merge.py")
        with tempfile.TemporaryDirectory() as work, tempfile.TemporaryDirectory() as other:
            with open(os.path.join(work, "jobs.txt"), "w") as f:
                f.write("A /c/A_c0.txt 0 A_g0_c0.root 1 0\nA /c/A_c1.txt 0 A_g0_c1.root 1 1\n")
            with open(os.path.join(work, storage.LOCATION_FILE), "w") as f:
                f.write(os.path.join(work, "nowhere") + "\n")
            for name in ("A_g0_c0.root", "A_g0_c1.root"):
                with open(os.path.join(other, name), "wb") as f:
                    f.write(b"x" * 2048)
            run = lambda *extra: subprocess.run([sys.executable, merge, work, "--check", *extra],
                                                capture_output=True, text=True)
            p = run()                                   # output_dir.txt: nothing there
            self.assertEqual(p.returncode, 1, p.stdout + p.stderr)
            self.assertIn("2 / 2 missing", p.stdout)
            p = run("--input-dir", other)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn(f"job outputs: {os.path.realpath(other)}", p.stdout.replace(other, os.path.realpath(other)))
            self.assertIn("complete", p.stdout)
            self.assertFalse(os.path.exists(os.path.join(work, "resubmit.txt")))
            p = run("--input-dir", os.path.join(other, "nope"))
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("is not a directory", p.stderr)



class PlotsTests(unittest.TestCase):
    DOC = {"output_dir": "out/plots",
           "inputs": {"data": "output/d.root", "qcd": "/abs/q.root", "sig": "output/S_%d.root"},
           "plots": [{"function": "nPV", "args": ["qcd", "data"]},
                     {"name": "Y_MET", "function": "Y", "args": [True, 2, 0.5, {"str": 'a"b'}, "sig"]}]}
    MACRO = """// comment
TString OutPath(const char *fmt, ...) { return ""; }
void nPV(const char *qcd = "../output/x.root",
         const char *data = "y.root") {
    TFile *f = new TFile(qcd, "READ");
    TH1F *h = (TH1F*)f->Get("METanalysis_Eta2p4_nPVgood");
    // TH1F *old = (TH1F*)f->Get("Old_name");
    TH1F *g = (TH1F*)f->Get(Form("%s_nPV", "x"));
}
void Y(bool a, int b, double c = 1, std::string d = "x,y", const char *e = "s") {
    h = f->Get("A"); h2 = f->Get("A");
}
"""

    def cfg(self, **change):
        return plots.PlotConfig(dict(self.DOC, **change), "t.json", base="/repo")

    def test_config(self):
        c = self.cfg()
        self.assertEqual(c.output_dir, "/repo/out/plots")
        self.assertEqual(c.inputs["data"], "/repo/output/d.root")
        self.assertEqual(c.inputs["qcd"], "/abs/q.root")
        self.assertEqual([p.name for p in c.plots], ["nPV", "Y_MET"])
        self.assertEqual(c.call(c.plots[0]), 'nPV("/abs/q.root", "/repo/output/d.root")')
        self.assertEqual(c.call(c.plots[1]), 'Y(true, 2, 0.5, "a\\"b", "/repo/output/S_%d.root")')
        self.assertEqual([p.name for p in c.select(["Y*"])], ["Y_MET"])
        with self.assertRaisesRegex(plots.PlotError, "no plot named zz"):
            c.select(["zz"])
        # null: an empty path, the sample is left out
        c = self.cfg(plots=[{"function": "nPV", "args": [None, "data"]}])
        self.assertEqual(c.call(c.plots[0]), 'nPV("", "/repo/output/d.root")')
        self.assertEqual(c.input_names(c.plots[0]), ["data"])
        self.assertEqual(c.legend_labels(c.plots[0], plots.signature(self.MACRO, "nPV")), {})

    def test_config_errors(self):
        for change, msg in ((dict(plots=[{"function": "f", "args": ["nope"]}]), "'nope' is not an input"),
                            (dict(plots=[{"function": "f"}, {"function": "f"}]), "used twice: f"),
                            (dict(plots="x"), "must be a list"),
                            (dict(plots=[{"function": "f", "args": [{"text": "x"}]}]), "literal string"),
                            (dict(plot=[]), "unknown key")):
            with self.assertRaisesRegex(plots.PlotError, msg, msg=change):
                self.cfg(**change)

    def test_labels(self):
        c = self.cfg(inputs={"data": {"file": "d.root", "label": "Data 2024"}, "qcd": "q.root",
                             "sig": {"file": "S_%d.root", "label": 'g "1600"'}},
                     lumi_label="7.98 fb^{-1}")
        self.assertEqual(c.inputs["data"], "/repo/d.root")
        self.assertEqual(c.labels, {"data": "Data 2024", "sig": 'g "1600"'})
        params = plots.signature(self.MACRO, "nPV")
        self.assertEqual([plots.param_name(q) for q in params], ["qcd", "data"])
        # keyed by the parameter receiving the input: here input "data" goes to parameter "data"
        self.assertEqual(c.legend_labels(c.plots[0], params), {"data": "Data 2024", "lumi": "7.98 fb^{-1}"})
        y = plots.signature(self.MACRO, "Y")
        self.assertEqual(c.legend_labels(c.plots[1], y), {"e": 'g "1600"', "lumi": "7.98 fb^{-1}"})
        line = " ".join(plots.root_command(c, c.plots[1], c.legend_labels(c.plots[1], y)))
        self.assertIn('gLabels["e"] = "g \\"1600\\"";', line)
        for bad in ({"data": {"label": "x"}}, {"data": {"file": "d", "style": 2}}, {"data": {"file": "d", "color": 2.5}}, {"data": 3}):
            with self.assertRaisesRegex(plots.PlotError, "must be a path", msg=bad):
                self.cfg(inputs=dict(self.DOC["inputs"], **bad))

    def test_signature_and_checks(self):
        self.assertEqual([p.has_default for p in plots.signature(self.MACRO, "nPV")], [True, True])
        params = plots.signature(self.MACRO, "Y")
        self.assertEqual(len(params), 5)                       # the comma inside "x,y" is not a separator
        self.assertEqual([p.has_default for p in params], [False, False, True, True, True])
        self.assertIsNone(plots.signature(self.MACRO, "Nope"))
        c = self.cfg(plots=[{"function": "nPV", "args": ["qcd", "data", "qcd"]},
                            {"name": "y1", "function": "Y", "args": [True]},
                            {"name": "y2", "function": "Y", "args": [True, 1]},
                            {"function": "Nope"}])
        errors = plots.check_calls(c, c.plots, self.MACRO)
        self.assertEqual(len(errors), 3)
        self.assertIn("nPV: nPV takes 0 to 2 argument(s), 3 given", errors[0])
        self.assertIn("y1: Y takes 2 to 5 argument(s), 1 given", errors[1])
        self.assertIn("no function Nope", errors[2])

    def test_histogram_names(self):
        self.assertEqual(plots.histogram_names(self.MACRO, "nPV"), ["METanalysis_Eta2p4_nPVgood"])
        self.assertEqual(plots.histogram_names(self.MACRO, "Y"), ["A"])

    def test_inputs_and_command(self):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "S_1600.root"), "w").close()
            open(os.path.join(d, "d.root"), "w").close()
            c = plots.PlotConfig(dict(self.DOC, inputs={"data": "d.root", "qcd": "q.root", "sig": "S_%d.root"}), "t", base=d)
            self.assertEqual(c.missing_inputs(c.plots[0]), [("qcd", os.path.join(d, "q.root"))])
            self.assertEqual(c.missing_inputs(c.plots[1]), [])           # the pattern matches S_1600.root
        argv = plots.root_command(self.cfg(), self.cfg().plots[0])
        self.assertEqual(argv[:4], ["root", "-l", "-b", "-q"])
        line = " ".join(argv)
        self.assertIn('LoadMacro("CombineHistos.C")', line)
        self.assertIn('gOutDir = "/repo/out/plots";', line)
        self.assertIn('gInterpreter->Declare(R"HSCP(void hscp_plot() { nPV("/abs/q.root", "/repo/output/d.root"); })HSCP")', line)

    def test_repository_plots_match_the_macro(self):
        """Every plot of cfg/plots.json calls a function of the macro with an accepted number of arguments."""
        c = plots.load()
        with open(c.macro) as f:
            macro = f.read()
        self.assertEqual(plots.check_calls(c, c.plots, macro), [])
        self.assertIn("TString OutPath(const char *fmt, ...)", macro)
        self.assertIn("const char *Label(const char *key, const char *fallback)", macro)
        self.assertIn("TFile *OpenOptional(const char *path)", macro)
        # the optional inputs are opened with OpenOptional, so that null (an empty path) leaves them out
        for name, params in (("Cutflows", ["QCD", "TTbar", "TTbarSemiLep", "Wjets"]), ("nPV", ["qcd", "ttbar", "wjets"]),
                             ("nHSCP", ["qcd", "ttbar", "wjets"]), ("FpixelPlot", ["wjets"]),
                             ("FpixSlices", ["wjets", "ttbar", "qcd"]), ("CompareMaping", ["ttbar", "wjets", "qcd"]),
                             ("CompareKinematics", ["gluino1400", "gluino2000", "gluino2600"])):
            body = plots.function_body(macro, name)
            self.assertIn("OpenOptional(", body, name)
            for p in params:
                self.assertNotIn(f"new TFile({p},", body, f"{name}: {p}")
        # every label key of a converted function, Label("key", ...) or legKeys, is one of its parameters or "lumi"
        import re
        for name in ("Cutflows", "nPV", "nHSCP", "FpixelPlot", "Acceptance_EtaSlice_DataMC", "CompareKinematics",
                     "YieldAfterSF"):
            body = plots.function_body(macro, name)
            keys = set(re.findall(r'Label\("(\w+)"', body))
            leg = re.search(r"legKeys\[\] = \{([^}]*)\}", body)
            if leg:
                keys |= set(re.findall(r'"(\w+)"', leg.group(1)))
            names = {plots.param_name(q) for q in plots.signature(macro, name)} | {"lumi"}
            self.assertTrue(keys, name)
            self.assertLessEqual(keys, names, name)
        for name in ("FpixelPlot", "CompareMaping", "nPV", "nHSCP", "FpixSlices", "Acceptance_EtaSlice_DataMC",
                     "CompareKinematics", "ShowSignalEfficiency", "YieldAfterSF", "Cutflows"):
            params = plots.signature(macro, name)
            self.assertIsNotNone(params, name)
            body = macro[macro.index(f"void {name}("):]
            body = body[:body.index("\nvoid ", 10)]
            live = [l for l in body.splitlines() if "../outputDisplay" in l and not l.strip().startswith("//")]
            self.assertEqual(live, [], name)



class StackTests(unittest.TestCase):
    DOC = {"inputs": {"data": {"file": "d.root", "label": "Data"}, "qcd": {"file": "q.root", "color": "kGreen-4"},
                      "tt": {"file": "t.root", "color": 4}, "sig": {"file": "s.root", "label": "g 1600"},
                      "pattern": "S_%d.root"},
           "stack": {"data": "data", "backgrounds": ["qcd", "tt"], "signals": ["sig"]}}

    def cfg(self, **change):
        return plots.PlotConfig(dict(self.DOC, **change), "t.json", base="/repo")

    def test_config(self):
        c = self.cfg()
        self.assertEqual(c.plots, [])                               # "plots" may be absent
        self.assertEqual(c.colors, {"qcd": "kGreen-4", "tt": 4})
        self.assertEqual(c.stack, {"data": "data", "backgrounds": ["qcd", "tt"], "signals": ["sig"]})
        for bad, msg in (({"backgrounds": ["nope"]}, "'nope' is not an input"), ({"data": 3}, "input name or null"),
                         ({"backgrounds": "qcd"}, "list of input names"), ({"other": []}, "expected")):
            with self.assertRaisesRegex(plots.PlotError, msg, msg=bad):
                self.cfg(stack=bad)

    def test_samples(self):
        c = self.cfg()
        self.assertEqual(stack.samples(c), stack.Samples("data", ["qcd", "tt"], ["sig"]))
        # the command line wins: '' for no data, lists for the others
        self.assertEqual(stack.samples(c, data="", backgrounds=["tt"]), stack.Samples(None, ["tt"], ["sig"]))
        self.assertEqual(stack.split_names(" qcd, tt ,"), ["qcd", "tt"])
        self.assertIsNone(stack.split_names(None))
        for kw, msg in ((dict(backgrounds=["nope"]), "not an input"), (dict(signals=["pattern"]), "pattern"),
                        (dict(backgrounds=["qcd", "qcd"]), "given twice")):
            with self.assertRaisesRegex(plots.PlotError, msg, msg=kw):
                stack.samples(c, **kw)
        with self.assertRaisesRegex(plots.PlotError, "no sample"):
            stack.samples(self.cfg(stack=None))

    def test_helpers(self):
        keys = ["A_Fpix", "A_eta", "B_Fpix", "EventCutflow"]
        self.assertEqual(stack.expand(["*_Fpix", "EventCutflow"], keys), ["A_Fpix", "B_Fpix", "EventCutflow"])
        with self.assertRaisesRegex(plots.PlotError, "no histogram named C_\\*"):
            stack.expand(["A_*", "C_*"], keys)
        self.assertEqual(stack.parse_color("kRed+2"), ("kRed", 2))
        self.assertEqual(stack.parse_color("kGreen - 4"), ("kGreen", -4))
        self.assertEqual(stack.parse_color(" kBlack"), ("kBlack", 0))
        self.assertEqual(stack.parse_color(7), 7)
        with self.assertRaises(plots.PlotError):
            stack.parse_color("red")
        cols = stack.colors(self.cfg(), ["qcd", "tt", "x", "y"], ["kBlue", "kRed"])
        self.assertEqual(cols, {"qcd": ("kGreen", -4), "tt": 4, "x": ("kBlue", 0), "y": ("kRed", 0)})
        self.assertEqual(stack.file_name("a/b c*d"), "a_b_c_d")
        self.assertEqual(stack.legend_label(self.cfg(), "sig", 100), "g 1600 #times 100")
        self.assertEqual(stack.legend_label(self.cfg(), "qcd"), "qcd")

    def test_repository_stack_block(self):
        c = plots.load()
        if c.stack:
            stack.samples(c)                                         # names are inputs, no pattern


if __name__ == "__main__":
    unittest.main()
