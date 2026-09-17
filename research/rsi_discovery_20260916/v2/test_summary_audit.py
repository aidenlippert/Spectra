"""The report must distinguish inaccessible evidence from verified evidence."""
import pytest

from research.rsi_discovery_20260916.ledger import Ledger
from research.rsi_discovery_20260916.v2 import summarize


def test_complete_audit_preserves_unfinished_operations(tmp_path):
    ledger = Ledger(tmp_path)
    ledger.blob("evidence")
    started = ledger.append("operation_started", label="unfinished")
    result = summarize.audit_available_evidence(ledger)
    assert result["chain_valid"]
    assert result["all_blob_bytes_verified"]
    assert result["blobs_verified"] == 1
    assert result["unfinished_operations"] == [started]
    ledger.db.close()


def test_offloaded_objects_are_never_counted_as_verified(tmp_path, monkeypatch):
    ledger = Ledger(tmp_path)
    key = ledger.blob("evidence")
    monkeypatch.setattr(summarize, "is_offloaded", lambda path: True)
    monkeypatch.setattr(ledger, "read_blob", lambda key: pytest.fail("offloaded read"))
    result = summarize.audit_available_evidence(ledger)
    assert result["chain_valid"]
    assert not result["all_blob_bytes_verified"]
    assert result["blobs_verified"] == 0
    assert [row["object"] for row in result["unavailable_blobs"]] == [key]
    ledger.db.close()


def test_corrupt_resident_object_refuses_audit(tmp_path):
    ledger = Ledger(tmp_path)
    key = ledger.blob("evidence")
    (tmp_path / "objects" / key).write_text("changed")
    with pytest.raises(ValueError, match="Corrupted evidence"):
        summarize.audit_available_evidence(ledger)
    ledger.db.close()


def test_broken_chain_refuses_audit(tmp_path):
    ledger = Ledger(tmp_path)
    ledger.append("evidence", value=1)
    ledger.db.execute("DROP TRIGGER no_update")
    ledger.db.execute("UPDATE events SET record='{}'")
    with pytest.raises(ValueError, match="Broken history chain"):
        summarize.audit_available_evidence(ledger)
    ledger.db.close()
