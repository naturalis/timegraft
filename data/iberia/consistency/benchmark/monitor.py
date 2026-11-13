"""Run a command; record wall time and the peak resident memory of its process tree (Linux /proc)."""
import json, os, subprocess, sys, time

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

out, cmd = sys.argv[1], sys.argv[2:]
start, peak = time.time(), 0
proc = subprocess.Popen(cmd)
while proc.poll() is None:
    peak = max(peak, tree_rss(proc.pid))
    time.sleep(2)
json.dump({"command": " ".join(cmd), "seconds": round(time.time() - start, 1), "peak_rss_mb": round(peak / 2**20),
           "returncode": proc.returncode}, open(out, "w"), indent=1)
