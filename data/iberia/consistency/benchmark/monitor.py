"""Run a command; record its wall time and peak memory (Linux).

    python monitor.py OUT.json COMMAND [ARGS...]

Two measures of peak memory are recorded:

- ``max_single_process_mb``: the largest resident set any one process of the command reached, over its
  whole lifetime, from the kernel's accounting (``getrusage(RUSAGE_CHILDREN)``, which covers every
  descendant that was waited for). Exact, but per process.
- ``peak_tree_rss_mb``: the largest sum of resident sets over the command's process tree, sampled from
  /proc every 0.2 s, from the moment it starts. Counts processes that run at the same time (two
  Snakemake jobs), but can miss a peak shorter than the interval.

``peak_rss_mb``, which the gate uses, is the larger of the two.
"""

import json
import os
import resource
import subprocess
import sys
import time

INTERVAL = 0.2


def tree_rss(root):
    children = {}
    for pid in filter(str.isdigit, os.listdir("/proc")):
        try:
            with open(f"/proc/{pid}/stat") as h:
                parts = h.read().rsplit(")", 1)[1].split()
            children.setdefault(int(parts[1]), []).append(int(pid))
        except OSError:
            pass
    total, stack = 0, [root]
    while stack:
        pid = stack.pop()
        stack.extend(children.get(pid, []))
        try:
            with open(f"/proc/{pid}/statm") as h:
                total += int(h.read().split()[1]) * os.sysconf("SC_PAGE_SIZE")
        except OSError:
            pass
    return total


def main():
    out, cmd = sys.argv[1], sys.argv[2:]
    start, peak = time.time(), 0
    proc = subprocess.Popen(cmd)
    while proc.poll() is None:
        peak = max(peak, tree_rss(proc.pid))
        time.sleep(INTERVAL)
    single = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024  # kB on Linux
    result = {
        "command": " ".join(cmd), "seconds": round(time.time() - start, 1),
        "peak_tree_rss_mb": round(peak / 2**20), "max_single_process_mb": round(single / 2**20),
        "peak_rss_mb": round(max(peak, single) / 2**20), "sample_interval_s": INTERVAL,
        "returncode": proc.returncode,
    }
    with open(out, "w") as handle:
        json.dump(result, handle, indent=1)
    return proc.returncode


if __name__ == "__main__":
    sys.exit(main())
