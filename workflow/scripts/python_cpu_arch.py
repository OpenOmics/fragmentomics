#!/usr/bin/env python3
"""
Detects CPU instruction set support and prints the name of the best
bwa-mem2 precompiled binary for this CPU.

Selection priority (best to worst):
    avx2 -> sse4.2 -> sse4.1 -> bwa-mem2 (upstream dispatcher shim)

AVX-512 is intentionally NOT considered. The upstream bwa-mem2.avx512bw
binary has known instability including segfaults on both Intel and AMD
CPUs, and additional issues on AMD Zen where the CPU reports avx512bw
support but the binary uses Intel-tuned code paths that fail at runtime.
See bwa-mem2 GitHub issues #50, #112, #160, #199, #254.

If none of the specific SIMD variants match (very unusual on modern
x86_64), fall back to the plain 'bwa-mem2' binary, which is the upstream
shim that performs its own runtime dispatch.
"""

import sys


def err(msg):
    print(msg, file=sys.stderr)


def fatal(msg):
    err(msg)
    sys.exit(1)


def read_cpu_flags():
    """Return the set of CPU flags from /proc/cpuinfo."""
    try:
        with open("/proc/cpuinfo", "r") as f:
            for line in f:
                if line.startswith("flags"):
                    _, _, values = line.partition(":")
                    return set(values.split())
    except OSError as e:
        fatal(f"Fatal error: unable to read /proc/cpuinfo: {e}")
    return set()


def select_binary(flags):
    """Pick the best bwa-mem2 variant the CPU can run."""
    # Ordered from most to least preferred.
    candidates = [
        ("avx2", "bwa-mem2.avx2"),
        ("sse4_2", "bwa-mem2.sse42"),
        ("sse4_1", "bwa-mem2.sse41"),
    ]
    for flag, binary in candidates:
        if flag in flags:
            return binary
    # Ultimate fallback: the upstream shim will make its own decision.
    return "bwa-mem2"


def main():
    flags = read_cpu_flags()

    if not flags:
        fatal(
            "Fatal error: failed to resolve CPU architecture, "
            "no CPU flags detected in /proc/cpuinfo."
        )

    # x86-only sanity check: if none of the baseline x86 flags are present,
    # we're almost certainly on ARM or another non-x86 architecture.
    if not flags & {"sse", "sse2", "sse4_1", "sse4_2", "avx", "avx2"}:
        err(
            "Fatal error: failed to resolve CPU architecture, "
            "neither Intel nor AMD detected."
        )
        err("  ├── Are you trying to run bwa-mem2 on a non-x86 CPU (ARM-based) cpu?")
        fatal("  └── Only Intel and AMD CPUs are supported, exiting now!")

    print(select_binary(flags))


if __name__ == "__main__":
    main()
