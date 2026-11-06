"""Small HTTP helpers with an on-disk JSON cache, so each query is made only once."""

from __future__ import annotations

import http.client
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

USER_AGENT = "timegraft (https://github.com/naturalis/timegraft)"


class ApiError(RuntimeError):
    def __init__(self, status: int, body: str):
        super().__init__(f"HTTP {status}: {body[:300]}")
        self.status = status
        self.body = body


def request_json(url: str, payload: dict | None = None, retries: int = 5, timeout: int = 120):
    data = json.dumps(payload).encode() if payload is not None else None
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    if data is not None:
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers, method="POST" if data is not None else "GET")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as err:
            body = err.read().decode(errors="replace")
            if err.code >= 500 and attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            raise ApiError(err.code, body) from None
        except (urllib.error.URLError, http.client.HTTPException, OSError):
            # connection resets and dropped tunnels surface as several exception types
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue
            raise


class JsonCache:
    """A dict persisted as one JSON file."""

    def __init__(self, path: str | Path | None):
        self.path = Path(path) if path else None
        self.data: dict = {}
        if self.path and self.path.exists():
            self.data = json.loads(self.path.read_text())

    def get(self, key: str):
        return self.data.get(key)

    def __contains__(self, key: str) -> bool:
        return key in self.data

    def set(self, key: str, value) -> None:
        self.data[key] = value

    def save(self) -> None:
        if self.path:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.path.write_text(json.dumps(self.data, indent=1, sort_keys=True))
