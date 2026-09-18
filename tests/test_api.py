"""API integration tests using FastAPI's TestClient."""
import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as test_client:
        yield test_client


PAYLOAD = {
    "primers": [
        {"name": "fwd", "sequence": "ATGCATGC"},
        {"name": "rev", "sequence": "GCATGCAT"},
        {"name": "probe", "sequence": "NNNNACGT"},
    ],
    "max_consecutive_threshold": 4,
    "three_prime_threshold": 3,
}


class TestCreateBatch:
    def test_create_returns_full_report(self, client):
        resp = client.post("/api/v1/batches", json=PAYLOAD)
        assert resp.status_code == 201
        data = resp.json()
        assert data["batch_id"] >= 1
        assert data["primer_count"] == 3
        assert data["pair_count"] == 6  # 3 self + 3 hetero
        assert data["max_consecutive_threshold"] == 4
        assert data["three_prime_threshold"] == 3
        assert len(data["results"]) == 6
        hetero = next(
            r
            for r in data["results"]
            if r["primer_a"] == "fwd" and r["primer_b"] == "rev"
        )
        assert hetero["kind"] == "hetero"
        assert hetero["max_consecutive"] == 8
        assert hetero["three_prime_consecutive"] == 8
        assert hetero["flagged"] is True
        assert hetero["reasons"] == [
            "MAX_CONSECUTIVE_COMPLEMENTARITY",
            "THREE_PRIME_END_COMPLEMENTARITY",
        ]
        assert data["flagged_count"] >= 1

    def test_sequences_normalized_to_uppercase(self, client):
        payload = {"primers": [{"name": "p", "sequence": "atgc"}]}
        resp = client.post("/api/v1/batches", json=payload)
        assert resp.status_code == 201
        assert resp.json()["primers"][0]["sequence"] == "ATGC"

    def test_default_thresholds(self, client):
        payload = {"primers": [{"name": "p", "sequence": "ATGC"}]}
        resp = client.post("/api/v1/batches", json=payload)
        assert resp.status_code == 201
        data = resp.json()
        assert data["max_consecutive_threshold"] == 4
        assert data["three_prime_threshold"] == 3

    def test_invalid_iupac_symbol_rejected(self, client):
        payload = {"primers": [{"name": "bad", "sequence": "ATGCX"}]}
        resp = client.post("/api/v1/batches", json=payload)
        assert resp.status_code == 422

    def test_duplicate_names_rejected(self, client):
        payload = {
            "primers": [
                {"name": "p", "sequence": "ATGC"},
                {"name": "p", "sequence": "GCAT"},
            ]
        }
        resp = client.post("/api/v1/batches", json=payload)
        assert resp.status_code == 422

    def test_empty_primer_list_rejected(self, client):
        resp = client.post("/api/v1/batches", json={"primers": []})
        assert resp.status_code == 422

    def test_invalid_threshold_rejected(self, client):
        payload = {
            "primers": [{"name": "p", "sequence": "ATGC"}],
            "max_consecutive_threshold": 0,
        }
        resp = client.post("/api/v1/batches", json=payload)
        assert resp.status_code == 422


class TestGetBatch:
    def test_roundtrip(self, client):
        created = client.post("/api/v1/batches", json=PAYLOAD).json()
        resp = client.get(f"/api/v1/batches/{created['batch_id']}")
        assert resp.status_code == 200
        fetched = resp.json()
        assert fetched == created

    def test_not_found(self, client):
        resp = client.get("/api/v1/batches/9999")
        assert resp.status_code == 404

    def test_list_batches(self, client):
        client.post("/api/v1/batches", json=PAYLOAD)
        client.post(
            "/api/v1/batches",
            json={"primers": [{"name": "solo", "sequence": "ATGC"}]},
        )
        resp = client.get("/api/v1/batches")
        assert resp.status_code == 200
        summaries = resp.json()
        assert len(summaries) == 2
        assert summaries[0]["batch_id"] > summaries[1]["batch_id"]
        assert summaries[0]["pair_count"] == 1  # newest first: single primer


class TestHealth:
    def test_health(self, client):
        resp = client.get("/health")
        assert resp.status_code == 200
        assert resp.json() == {"status": "ok"}
