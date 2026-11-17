"""Merging API caches into the benchmark cache (#87)."""

import gzip
import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "data/iberia/consistency/benchmark/merge_caches.py"
ANSWER = {"v": 2, "age": 12.0, "api_version": "20260820", "found": [1, 2]}


def _run(out, *caches):
    return subprocess.run([sys.executable, str(SCRIPT), str(out), *map(str, caches)], capture_output=True, text=True)


def _cache(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))
    return path


def test_merges_and_writes_the_same_bytes(tmp_path):
    a = _cache(tmp_path / "a.json", {"mrca:1+2": ANSWER, "col:x": {"usage": 1}})
    b = _cache(tmp_path / "b.json", {"mrca:1+2": ANSWER, "mrca:3+4": {**ANSWER, "age": 5.0}})
    assert _run(tmp_path / "m.json.gz", a, b).returncode == 0
    merged = json.loads(gzip.open(tmp_path / "m.json.gz").read())
    assert set(merged) == {"mrca:1+2", "mrca:3+4", "col:x"}
    first = (tmp_path / "m.json.gz").read_bytes()
    assert _run(tmp_path / "m.json.gz", b, a, tmp_path / "m.json.gz").returncode == 0
    assert (tmp_path / "m.json.gz").read_bytes() == first  # order of the inputs and time do not matter
    assert _run(tmp_path / "elsewhere.json.gz", a, b).returncode == 0
    assert (tmp_path / "elsewhere.json.gz").read_bytes() == first  # nor the name of the output


def test_refuses_conflicts_perturbed_copies_and_old_formats(tmp_path):
    a = _cache(tmp_path / "a.json", {"mrca:1+2": ANSWER})
    other = _cache(tmp_path / "b.json", {"mrca:1+2": {**ANSWER, "age": 13.0}})
    result = _run(tmp_path / "m.json", a, other)
    assert result.returncode != 0 and "differently" in result.stderr
    noise = _cache(tmp_path / "exp/noise/s1/cache.json", {"mrca:1+2": ANSWER})
    result = _run(tmp_path / "m.json", a, noise)
    assert result.returncode != 0 and "perturbed" in result.stderr
    for legacy in ({"age": None, "found": [1]}, {"v": 1, "age": 12.0, "api_version": "20260820", "found": [1, 2]}):
        # the version marker decides, as in timegraft.timetree.is_legacy(), not the other fields
        old = _cache(tmp_path / "old.json", {"mrca:1+2": legacy})
        result = _run(tmp_path / "m.json", old)
        assert result.returncode != 0 and "format 2" in result.stderr
