"""
Grid proxy for the Condor jobs (condor/submit.py --proxy).

    auto  (default) a proxy is given to the jobs if one of their input files is on remote storage
    none  no proxy
    PATH  this proxy file

Condor copies the proxy named by x509userproxy to the worker node and sets X509_USER_PROXY there,
which is what ROOT (davix, xrootd) reads: runJob.sh needs nothing else.
"""

import os
import shutil
import subprocess

MARGIN = 3600            # seconds of proxy lifetime required beyond --runtime


class ProxyError(RuntimeError):
    pass


def is_remote(path):
    """True for a URL that needs credentials (davs://, root://, ...), False for a local path or file://."""
    return "://" in path and not path.startswith("file://")


def default_path():
    """Where voms-proxy-init writes the proxy: $X509_USER_PROXY, else /tmp/x509up_u<uid>."""
    return os.environ.get("X509_USER_PROXY") or f"/tmp/x509up_u{os.getuid()}"


def timeleft(path):
    """Seconds of lifetime left on the proxy, or None if voms-proxy-info is not available."""
    exe = shutil.which("voms-proxy-info")
    if not exe:
        return None
    try:
        p = subprocess.run([exe, "--file", path, "--timeleft"], capture_output=True, text=True, timeout=60)
    except subprocess.TimeoutExpired:
        raise ProxyError(f"voms-proxy-info timed out on {path}")
    try:
        return int(p.stdout.split()[-1])
    except (ValueError, IndexError):
        raise ProxyError(f"voms-proxy-info cannot read {path}: {(p.stderr or p.stdout).strip()}")


def resolve(mode, files, runtime, timeleft_fn=timeleft):
    """(path, seconds left) of the proxy to give the jobs, (None, None) if none is needed.

    seconds left is None when voms-proxy-info is not available (the lifetime is then not checked).
    Raises ProxyError if a proxy is needed but missing, or expires before runtime + MARGIN.
    """
    if mode == "none":
        return None, None
    if mode == "auto":
        if not any(is_remote(f) for f in files):
            return None, None
        path = default_path()
    else:
        path = mode
    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.isfile(path):
        raise ProxyError(f"no grid proxy at {path}: run  voms-proxy-init -voms cms -valid 192:00"
                         + ("  or pass --proxy none" if mode == "auto" else ""))
    left = timeleft_fn(path)
    if left is not None and left < runtime + MARGIN:
        raise ProxyError(f"the proxy {path} has {left // 3600} h {left % 3600 // 60} min left, less than "
                         f"--runtime + 1 h: renew it with  voms-proxy-init -voms cms -valid 192:00")
    return path, left


def sub_line(path):
    """The line of hscp.sub giving the proxy to the jobs ('' for none)."""
    return f"x509userproxy         = {path}\n" if path else ""


def describe(path, left):
    """One line for the submit.py summary."""
    if not path:
        return "none"
    if left is None:
        return f"{path} (lifetime not checked: voms-proxy-info not found)"
    return f"{path} ({left // 3600} h left)"
