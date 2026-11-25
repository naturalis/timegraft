"""Run a command; record its wall time and peak memory (Linux and macOS).

    python monitor.py OUT.json COMMAND [ARGS...]

Two measures of peak memory are recorded:

- ``max_single_process_mb``: the largest resident set any one process of the command reached, over its
  whole lifetime, from the kernel's accounting (``getrusage(RUSAGE_CHILDREN)``, which covers every
  descendant that was waited for). Exact, but per process.
- ``peak_tree_rss_mb``: the largest sum of resident sets over the command's process tree, sampled from
  /proc every 0.2 s (on macOS, which has no /proc, from ``ps``), from the moment it starts. Counts processes that run at the same time (two
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
# ru_maxrss is in kilobytes on Linux and in bytes on macOS
MAXRSS_UNIT = 1 if sys.platform == "darwin" else 1024


def tree_rss(root):
    """Sum of the resident sets of ``root`` and its descendants, in bytes."""
    if not os.path.isdir("/proc"):
        return tree_rss_ps(root)
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


def tree_rss_ps(root):
    """As tree_rss, from ``ps`` (rss in kilobytes on macOS and Linux), for systems without /proc."""
    try:
        lines = subprocess.run(["ps", "-A", "-o", "pid=,ppid=,rss="], capture_output=True, text=True,
                               check=True).stdout.splitlines()
    except (OSError, subprocess.CalledProcessError):
        return 0
    children, rss = {}, {}
    for line in lines:
        fields = line.split()
        if len(fields) != 3 or not all(f.isdigit() for f in fields):
            continue
        pid, ppid, kb = map(int, fields)
        children.setdefault(ppid, []).append(pid)
        rss[pid] = kb * 1024
    total, stack = 0, [root]
    while stack:
        pid = stack.pop()
        stack.extend(children.get(pid, []))
        total += rss.get(pid, 0)
    return total


def main():
    out, cmd = sys.argv[1], sys.argv[2:]
    start, peak = time.time(), 0
    proc = subprocess.Popen(cmd)
    while proc.poll() is None:
        peak = max(peak, tree_rss(proc.pid))
        time.sleep(INTERVAL)
    single = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * MAXRSS_UNIT
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
