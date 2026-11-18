#!/usr/bin/env python3
"""
Write the list of the ROOT files of a production (output_1.root, output_2.root, ...) with gfal-ls.

    python3 python/MakeFileList.py <dir> [<dir> ...] [-o list.txt] [-a]

    # a CRAB task: every output_<N>.root of its 0000/, 0001/, ... sub-directories
    python3 python/MakeFileList.py \\
        HSCP-Gluino_Par-M-1600_TuneCP5_13p6TeV_madgraphMLM-pythia8/Analysis_HSCP-Gluino_Par-M-1600_CodeV19p6 \\
        -o lists/Gluino_1600_V19p6.txt

<dir> is a URL (davs://, root://, file://, ...) or a path relative to the production directory
(--prod-dir, or $HSCP_PROD_DIR). It is searched recursively
(--depth levels, default 3: task/timestamp/0000), skipping the log/ and failed/ directories of CRAB.
The list holds one full URL per line, sorted by directory then by file number, and can be used as
the 'list' of a sample in cfg/datasets.json (give its absolute path there).
-o replaces an existing list; with -a the files are added at its end instead, leaving out the ones
it already holds.

Needs gfal-ls in the PATH and, for grid storage, a valid proxy (voms-proxy-init -voms cms).
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

SKIP_DIRS = {"log", "failed"}


class GfalError(RuntimeError):
    pass


PYTHON_VARS = ("PYTHONHOME", "PYTHONPATH")
# what a Python interpreter prints when PYTHONHOME belongs to another Python
FOREIGN_HOME = ("No module named 'encodings'", "init_fs_encoding")
_gfal = {}


def _interpreter(script):
    """The interpreter a script will run, from its first line (following "env"), or None."""
    try:
        with open(script, "rb") as f:
            first = f.readline(512).decode("utf-8", "replace").strip()
    except (OSError, TypeError):
        return None, ""
    if not first.startswith("#!"):
        return None, "(not a script: a compiled program or launcher)"
    words = first[2:].split()
    if not words:
        return None, first
    if os.path.basename(words[0]) == "env":
        names = [w for w in words[1:] if not w.startswith("-") and "=" not in w]
        return (shutil.which(names[0]) if names else None), first
    return words[0], first


def _clean(env):
    return {k: v for k, v in env.items() if k not in PYTHON_VARS}


def gfal_environment():
    """Environment for gfal-ls: ours, or ours without PYTHONHOME and PYTHONPATH.

    An LCG view (env.sh) sets PYTHONHOME and PYTHONPATH for its own Python. They are only valid for
    that interpreter: a gfal-ls running on another one (the system Python, "#!/usr/bin/python", or
    "#!/usr/bin/env python" falling through to it) would look for its standard library in the LCG tree
    and stop with "No module named 'encodings'". So the two variables are removed when the interpreter
    of gfal-ls lies outside PYTHONHOME, and kept otherwise: a gfal-ls of the LCG view itself may need
    that PYTHONPATH. gfal_ls retries without them if this guess was wrong.
    """
    if "env" not in _gfal:
        env = dict(os.environ)
        path = shutil.which("gfal-ls")
        interp, first = _interpreter(path)
        _gfal.update(path=path, first=first, interp=interp, env=env)
        home = env.get("PYTHONHOME")
        if home and interp and "python" in os.path.basename(interp):
            inside = os.path.realpath(interp).startswith(os.path.realpath(home).rstrip("/") + "/")
            if not inside:
                _strip(f"its interpreter {interp} is not the Python of PYTHONHOME")
    return _gfal["env"]


def _strip(why):
    _gfal["env"] = _clean(os.environ)
    _gfal["stripped"] = True
    print(f"MakeFileList: calling {_gfal.get('path')} without PYTHONHOME and PYTHONPATH ({why})", file=sys.stderr)


def _run_gfal(url, timeout):
    try:
        return subprocess.run(["gfal-ls", "-l", url], capture_output=True, text=True, timeout=timeout,
                              env=gfal_environment())
    except subprocess.TimeoutExpired:
        raise GfalError(f"gfal-ls timed out after {timeout} s on {url}")


def gfal_ls(url, timeout):
    """Entries of a directory as (name, is_dir), from 'gfal-ls -l'."""
    p = _run_gfal(url, timeout)
    if p.returncode != 0 and not _gfal.get("stripped") and any(s in p.stderr for s in FOREIGN_HOME):
        _strip("it stopped on the PYTHONHOME of another Python")
        p = _run_gfal(url, timeout)
    if p.returncode != 0:
        where = (f"\n(gfal-ls is {_gfal.get('path')}, first line {_gfal.get('first')!r}, "
                 f"interpreter {_gfal.get('interp')}, PYTHONHOME {'removed' if _gfal.get('stripped') else 'kept'})")
        raise GfalError(f"gfal-ls failed on {url}:\n{p.stderr.strip() or p.stdout.strip()}{where}")
    entries = []
    for line in p.stdout.splitlines():
        fields = line.split()
        if len(fields) < 2:
            continue
        # mode, links, uid, gid, size, date (3 fields), name; only the mode and the name are used
        entries.append((fields[-1], fields[0].startswith("d")))
    return entries


def find_files(url, pattern, depth, timeout, verbose=False):
    """Full URLs of the files below url whose name matches pattern, sorted by directory then by number."""
    url = url.rstrip("/")
    if verbose:
        print(f"  listing {url}", file=sys.stderr)
    files, subdirs = [], []
    for name, is_dir in gfal_ls(url, timeout):
        base = name.rstrip("/").rsplit("/", 1)[-1]
        if is_dir:
            if base not in SKIP_DIRS and base not in (".", ".."):
                subdirs.append(base)
            continue
        m = pattern.fullmatch(base)
        if m:
            # gfal-ls on a file prints the URL it was given: keep it as it is
            number = int(m.group(1)) if m.groups() else 0
            files.append((number, base, name if "://" in name else f"{url}/{base}"))
    found = [full for _, _, full in sorted(files)]
    if depth > 0:
        for sub in sorted(subdirs):
            found += find_files(f"{url}/{sub}", pattern, depth - 1, timeout, verbose)
    return found


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", metavar="dir", help="URL, or path relative to the production directory")
    ap.add_argument("-o", "--output", help="file list to write (default: standard output)")
    ap.add_argument("-a", "--append", action="store_true",
                    help="add to the --output list instead of replacing it; files already in it are not repeated")
    ap.add_argument("--prod-dir", help="production directory for relative paths (default: $HSCP_PROD_DIR)")
    ap.add_argument("--pattern", default=r"output_(\d+)\.root",
                    help="regular expression of the file names; its first group is the file number (default: %(default)s)")
    ap.add_argument("--depth", type=int, default=3, help="levels of sub-directories searched (default: %(default)s)")
    ap.add_argument("--timeout", type=int, default=120, help="seconds allowed to each gfal-ls call (default: %(default)s)")
    ap.add_argument("-v", "--verbose", action="store_true", help="print every directory listed")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.append and not args.output:
        sys.exit("MakeFileList: --append needs --output")
    if not shutil.which("gfal-ls"):
        sys.exit("MakeFileList: gfal-ls not found in the PATH")
    try:
        pattern = re.compile(args.pattern)
    except re.error as e:
        sys.exit(f"MakeFileList: bad --pattern ({e})")

    prod = args.prod_dir or os.environ.get("HSCP_PROD_DIR")
    files = []
    for d in args.dirs:
        if "://" not in d and not d.startswith("/"):
            if not prod:
                sys.exit(f"MakeFileList: '{d}' is a relative path: give --prod-dir or set HSCP_PROD_DIR, or pass a URL")
            d = prod.rstrip("/") + "/" + d
        elif d.startswith("/"):
            d = "file://" + d
        try:
            found = find_files(d, pattern, args.depth, args.timeout, args.verbose)
        except GfalError as e:
            sys.exit(f"MakeFileList: {e}")
        if not found:
            sys.exit(f"MakeFileList: no file matching '{args.pattern}' below {d}")
        print(f"{len(found):5d} files  {d}", file=sys.stderr)
        files += found

    if args.output:
        out_dir = os.path.dirname(os.path.abspath(args.output))
        os.makedirs(out_dir, exist_ok=True)
        existing = ""
        if args.append and os.path.isfile(args.output):
            with open(args.output) as f:
                existing = f.read()
            known = {line.strip() for line in existing.splitlines()}
            new = [f for f in files if f not in known]
            print(f"{len(files) - len(new)} files already in {args.output}", file=sys.stderr)
            files = new
        with open(args.output, "a" if args.append else "w") as f:
            if existing and not existing.endswith("\n"):
                f.write("\n")
            f.write("".join(name + "\n" for name in files))
        print(f"{len(files)} files {'added' if args.append else 'written'} to {args.output}", file=sys.stderr)
    else:
        sys.stdout.write("".join(name + "\n" for name in files))


if __name__ == "__main__":
    main()
