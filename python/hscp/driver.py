"""
Running macros/macro.cc: the one place that knows its argument order.

    macro(dataset, fileList, outFile, fillNosel, ptcut, etabins, ihbins, pbins, massbins, fpixbins)

Used by launch.py locally; condor/runJob.sh calls root with the same arguments in the job.
"""

import glob
import os
import subprocess

from .config import MACROS_DIR, Binning

SELECTOR_SOURCE = "HSCPSelector.C"
TOOLS_LIB = "../libTools.so"          # relative to macros/, as in macro.cc


def _num(x):
    return str(int(x)) if float(x).is_integer() else str(x)


def macro_call(dataset, file_list, out_file, fill_nosel, binning):
    """The C++ call passed to root, e.g. macro.cc("X","","X.root",true,70,96,400,1000,400,20)."""
    b = Binning(*binning)
    return ('macro.cc("{}","{}","{}",{},{},{},{},{},{},{})'.format(
        dataset, file_list, out_file, "true" if fill_nosel else "false",
        _num(b.ptcut), b.etabins, b.ihbins, b.pbins, b.massbins, b.fpixbins))


def root_command(dataset, file_list, out_file, fill_nosel, binning):
    """argv of the root process running one dataset (to be started in the macros/ directory)."""
    return ["root", "-l", "-b", "-q", macro_call(dataset, file_list, out_file, fill_nosel, binning)]


def build_command():
    """argv compiling the selector once with ACLiC, started in macros/ (exit code 1 if it fails).
    libTools.so is loaded first, as macro.cc does, so that ACLiC links the selector against it."""
    load = 'gSystem->Load("{}")'.format(TOOLS_LIB)
    compile_ = 'gSystem->Exit(gSystem->CompileMacro("{}","k") ? 0 : 1)'.format(SELECTOR_SOURCE)
    return ["root", "-l", "-b", "-q", "-e", load, "-e", compile_]


def build_selector(macros_dir=MACROS_DIR, log=None):
    """Compile the selector in macros_dir. Returns the exit code."""
    with open(log, "w") if log else open(os.devnull, "w") as out:
        return subprocess.run(build_command(), cwd=macros_dir, stdout=out, stderr=subprocess.STDOUT).returncode


def clean_build(macros_dir=MACROS_DIR):
    """Remove the ACLiC products of the selector."""
    for f in glob.glob(os.path.join(macros_dir, "HSCPSelector_C*")):
        os.remove(f)


def run_dataset(dataset, file_list, out_file, fill_nosel, binning, log, macros_dir=MACROS_DIR, threads=None):
    """Run one dataset; stdout and stderr go to log. Returns the exit code (non-zero also if no output was written)."""
    env = dict(os.environ)
    if threads is not None:
        env["HSCP_NTHREADS"] = str(threads)
    with open(log, "w") as out:
        rc = subprocess.run(root_command(dataset, file_list, out_file, fill_nosel, binning),
                            cwd=macros_dir, env=env, stdout=out, stderr=subprocess.STDOUT).returncode
    if rc == 0 and not (os.path.isfile(out_file) and os.path.getsize(out_file) > 0):
        rc = 2
    return rc
