#!/usr/bin/env python3
"""
Select the appropriate bwa-mem2 binary based on CPU vendor.

Prints the name of the binary to stdout:
- 'bwa-mem2'      on Intel CPUs (lets it pick the best precompiled binary).
- 'bwa-mem2.avx2' on AMD CPUs   (compatible with Zen3/Zen4).

Exits with an error on non-x86 (e.g. ARM) or unknown CPUs.
"""

import sys


def err(msg):
    print(msg, file=sys.stderr)


def fatal(msg, code=1):
    err(msg)
    sys.exit(code)


def get_cpu_vendor():
    """Read /proc/cpuinfo and return the vendor_id string, or None."""
    try:
        with open("/proc/cpuinfo", "r") as f:
            for line in f:
                if line.startswith("vendor_id"):
                    # Format: "vendor_id       : GenuineIntel"
                    return line.split(":", 1)[1].strip()
    except FileNotFoundError:
        return None
    return None


def select_binary():
    """Return the name of the bwa-mem2 binary to use for this CPU."""
    vendor = get_cpu_vendor()

    if vendor == "GenuineIntel":
        return "bwa-mem2"
    elif vendor == "AuthenticAMD":
        return "bwa-mem2.avx2"
    else:
        err(
            "Fatal error: failed to resolve CPU architecture, neither Intel nor AMD detected."
        )
        err("  ├── Are you trying to run bwa-mem2 on a non-x86 CPU (ARM-based) cpu?")
        fatal("  └── Only Intel and AMD CPUs are supported, exiting now!")


def main():
    print(select_binary())


if __name__ == "__main__":
    main()
