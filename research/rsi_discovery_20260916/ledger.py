"""Content-addressed evidence and append-only, crash-visible cost records."""
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import resource
import sqlite3
import time


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


class Ledger:
    def __init__(self, directory):
        self.root = Path(directory).resolve()
        (self.root / "objects").mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(self.root / "history.sqlite")
        self.db.executescript("""
            PRAGMA journal_mode=WAL;
            CREATE TABLE IF NOT EXISTS events (
                seq INTEGER PRIMARY KEY, kind TEXT NOT NULL,
                record TEXT NOT NULL, previous TEXT NOT NULL, hash TEXT NOT NULL UNIQUE);
            CREATE TRIGGER IF NOT EXISTS no_update BEFORE UPDATE ON events
                BEGIN SELECT RAISE(ABORT, 'append-only history'); END;
            CREATE TRIGGER IF NOT EXISTS no_delete BEFORE DELETE ON events
                BEGIN SELECT RAISE(ABORT, 'append-only history'); END;
        """)

    def blob(self, content):
        raw = content.encode() if isinstance(content, str) else content
        key = hashlib.sha256(raw).hexdigest()
        path = self.root / "objects" / key
        if path.exists():
            if path.read_bytes() != raw:
                raise ValueError("Evidence hash collision or corrupted object")
        else:
            with path.open("xb") as stream:
                stream.write(raw)
        return key

    def read_blob(self, key):
        if len(key) != 64 or any(c not in "0123456789abcdef" for c in key):
            raise ValueError("Invalid content address")
        raw = (self.root / "objects" / key).read_bytes()
        if hashlib.sha256(raw).hexdigest() != key:
            raise ValueError("Corrupted evidence")
        return raw

    def append(self, kind, **record):
        record = {"timestamp_unix": time.time(), **record}
        encoded = canonical(record).decode()
        with self.db:
            self.db.execute("BEGIN IMMEDIATE")
            last = self.db.execute("SELECT hash FROM events ORDER BY seq DESC LIMIT 1").fetchone()
            previous = last[0] if last else "0" * 64
            key = digest([kind, record, previous])
            self.db.execute("INSERT INTO events(kind,record,previous,hash) VALUES(?,?,?,?)",
                            (kind, encoded, previous, key))
        return key

    def events(self, kind=None):
        query = "SELECT kind,record,previous,hash FROM events"
        args = ()
        if kind:
            query += " WHERE kind=?"
            args = (kind,)
        return [{"kind": k, **json.loads(r), "event_hash": h, "previous": p}
                for k, r, p, h in self.db.execute(query + " ORDER BY seq", args)]

    def audit(self):
        previous = "0" * 64
        for k, r, p, h in self.db.execute("SELECT kind,record,previous,hash FROM events ORDER BY seq"):
            if p != previous or h != digest([k, json.loads(r), p]):
                raise ValueError("Broken history chain")
            previous = h
        for path in (self.root / "objects").iterdir():
            self.read_blob(path.name)
        starts = {e["event_hash"] for e in self.events("operation_started")}
        ends = {e["operation"] for e in self.events("operation_finished")}
        attempts = {e["id"] for e in self.events("attempt_started")}
        completed = {e["id"] for e in self.events("attempt_finished")}
        return {"chain_valid": True, "unfinished_operations": sorted(starts - ends),
                "unfinished_attempts": sorted(attempts - completed),
                "events": self.db.execute("SELECT count(*) FROM events").fetchone()[0]}

    @contextmanager
    def measure(self, label, **metadata):
        start = time.monotonic()
        before = resource.getrusage(resource.RUSAGE_SELF)
        child = resource.getrusage(resource.RUSAGE_CHILDREN)
        op = self.append("operation_started", label=label, metadata=metadata)
        result = {"status": "completed"}
        try:
            yield result
        except BaseException as error:
            result.update(status="failed", error=f"{type(error).__name__}: {error}")
            raise
        finally:
            after = resource.getrusage(resource.RUSAGE_SELF)
            after_child = resource.getrusage(resource.RUSAGE_CHILDREN)
            self.append("operation_finished", operation=op, label=label, metadata=metadata,
                        wall_seconds=time.monotonic() - start,
                        self_cpu_seconds=after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
                        child_cpu_seconds=after_child.ru_utime + after_child.ru_stime - child.ru_utime - child.ru_stime,
                        result=result)

    def costs(self):
        rows = self.events("operation_finished")
        # Operations are deliberately non-nested; a sum counts each stage once.
        return {"metered_wall_seconds": sum(r["wall_seconds"] for r in rows),
                "self_cpu_seconds": sum(r["self_cpu_seconds"] for r in rows),
                "child_cpu_seconds": sum(r["child_cpu_seconds"] for r in rows),
                "operations": len(rows), "failed_operations": sum(r["result"]["status"] in ("failed", "refuted", "timeout") for r in rows),
                "model_usage": [r["result"].get("usage") for r in rows if r["label"] == "model_proposal"],
                "money_charged_by_runner": 0,
                "subscription_marginal_price": None}
