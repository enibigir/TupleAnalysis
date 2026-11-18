"""
Where the outputs of the Condor jobs go (condor/submit.py --out-dir, condor/merge.py).

A location is a local directory or a URL:
    /uscms_data/d1/me/hscp/<tag>                  local: Condor transfers the outputs there
    root://cmseos.fnal.gov//store/user/me/hscp/<tag>  remote: each job copies its output itself
                                                  (xrdcp for root://, gfal-copy for any other URL)
Remote directories are listed with "xrdfs ls -l" (root://) or "gfal-ls -l" (other URLs).
"""

import os
import subprocess

DEFAULT_SUBDIR = "out"            # condor/work/<tag>/out when no --out-dir is given
LOCATION_FILE = "output_dir.txt"  # written by submit.py in the working directory, read by merge.py


class StorageError(RuntimeError):
    pass


def is_url(path):
    return "://" in path


def resolve(out_dir, tag):
    """The output location of a submission: <out_dir>/<tag>, with $VARS and ~ expanded; local paths made absolute."""
    out_dir = os.path.expanduser(os.path.expandvars(out_dir))
    if "$" in out_dir:
        raise StorageError(f"undefined variable in --out-dir {out_dir}")
    if is_url(out_dir):
        if out_dir.startswith("root://") and "//" not in out_dir[len("root://"):]:
            raise StorageError(f"{out_dir}: an XRootD URL needs a double slash before the path, "
                               "e.g. root://cmseos.fnal.gov//store/user/...")
        return out_dir.rstrip("/") + "/" + tag
    return os.path.join(os.path.abspath(out_dir), tag)


def join(location, name):
    return location.rstrip("/") + "/" + name if is_url(location) else os.path.join(location, name)


def split_root_url(url):
    """root://host[:port]//path -> ("root://host[:port]", "/path")."""
    rest = url[len("root://"):]
    host, _, path = rest.partition("/")
    return "root://" + host, "/" + path.lstrip("/")


def parse_xrdfs(text):
    """{file name: size} from "xrdfs ls -l": the path is the last field, the size the last all-digit one before it."""
    out = {}
    for line in text.splitlines():
        fields = line.split()
        if len(fields) < 2 or fields[0].startswith("d"):
            continue
        sizes = [f for f in fields[:-1] if f.isdigit()]
        if sizes:
            out[fields[-1].rstrip("/").rsplit("/", 1)[-1]] = int(sizes[-1])
    return out


def parse_gfal(text):
    """{file name: size} from "gfal-ls -l": mode, links, uid, gid, size, date (3 fields), name."""
    out = {}
    for line in text.splitlines():
        fields = line.split()
        if len(fields) < 9 or fields[0].startswith("d") or not fields[4].isdigit():
            continue
        out[fields[-1].rstrip("/").rsplit("/", 1)[-1]] = int(fields[4])
    return out


def _missing_dir(stderr):
    s = stderr.lower()
    return "no such file" in s or "not found" in s or "does not exist" in s


def listing(location, timeout=600):
    """{file name: size in bytes} of the files in location ({} if the directory does not exist yet)."""
    if not is_url(location):
        if not os.path.isdir(location):
            return {}
        return {e.name: e.stat().st_size for e in os.scandir(location) if e.is_file()}
    if location.startswith("root://"):
        host, path = split_root_url(location)
        cmd, parse, env = ["xrdfs", host, "ls", "-l", path], parse_xrdfs, None
    else:
        import MakeFileList                    # python/ is in sys.path of the scripts using this module
        cmd, parse, env = ["gfal-ls", "-l", location], parse_gfal, MakeFileList.gfal_environment()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
    except FileNotFoundError:
        raise StorageError(f"{cmd[0]} not found in the PATH (source env.sh)")
    except subprocess.TimeoutExpired:
        raise StorageError(f"{' '.join(cmd)} timed out after {timeout} s")
    if p.returncode != 0:
        if _missing_dir(p.stderr + p.stdout):
            return {}
        raise StorageError(f"{' '.join(cmd)} failed:\n{(p.stderr or p.stdout).strip()}")
    return parse(p.stdout)


def read_location(work):
    """The output location recorded by submit.py, or <work>/out for a working directory made before it was."""
    path = os.path.join(work, LOCATION_FILE)
    if os.path.isfile(path):
        with open(path) as f:
            return f.read().strip()
    return os.path.join(work, DEFAULT_SUBDIR)
