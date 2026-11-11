#!/bin/bash
# Provide ROOT without CMSSW: use the ROOT already in the PATH, otherwise source the
# LCG view from CVMFS matching the OS (el8 or el9) and the architecture.
#
#   source env.sh            returns non-zero if no ROOT could be set up
#   HSCP_LCG_RELEASE=LCG_107 source env.sh     picks another LCG release
#
# Used by compile.sh, the launch*.sh scripts and condor/runJob.sh.

hscp_setup_root() {
    if command -v root >/dev/null 2>&1 && command -v root-config >/dev/null 2>&1; then
        echo "[env] using installed ROOT $(root-config --version) at $(command -v root)"
        return 0
    fi

    local lcg_release="${HSCP_LCG_RELEASE:-LCG_106}"
    local arch os_id os_like os_major platform setup

    arch=$(uname -m)                         # x86_64 or aarch64
    if [[ -r /etc/os-release ]]; then
        # shellcheck disable=SC1091
        . /etc/os-release
        os_id="${ID:-}"; os_like="${ID_LIKE:-}"; os_major="${VERSION_ID%%.*}"
    else
        echo "[env] no ROOT in PATH and no /etc/os-release: cannot pick an LCG view" >&2
        return 1
    fi

    # el8 / el9 covers RHEL, AlmaLinux, Rocky and CentOS Stream
    if [[ "$os_id" == "rhel" || "$os_like" == *rhel* || "$os_like" == *fedora* ]]; then
        case "$os_major" in
            9) platform="${arch}-el9-gcc13-opt" ;;
            8) platform="${arch}-el8-gcc13-opt" ;;
            *) echo "[env] unsupported EL major version: $os_major" >&2; return 1 ;;
        esac
    else
        echo "[env] unsupported OS: ID=$os_id ID_LIKE=$os_like" >&2
        return 1
    fi

    setup="/cvmfs/sft.cern.ch/lcg/views/${lcg_release}/${platform}/setup.sh"
    if [[ ! -r "$setup" ]]; then
        echo "[env] LCG view not found: $setup (is CVMFS mounted?)" >&2
        return 1
    fi

    echo "[env] sourcing $setup"
    # shellcheck disable=SC1090
    source "$setup"

    if ! command -v root >/dev/null 2>&1; then
        echo "[env] ROOT still not available after sourcing $setup" >&2
        return 1
    fi
    echo "[env] using LCG ROOT $(root-config --version)"
}

hscp_setup_root
