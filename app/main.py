"""PCR primer-dimer screening API."""
from __future__ import annotations

import json
from collections.abc import AsyncIterator
from contextlib import asynccontextmanager
from datetime import datetime, timezone

from fastapi import FastAPI, HTTPException

from app.db import get_conn, init_db
from app.dimer import screen_primers
from app.schemas import (
    BatchRequest,
    BatchResponse,
    BatchSummary,
    PairResult,
    PrimerInput,
)


@asynccontextmanager
async def lifespan(_: FastAPI) -> AsyncIterator[None]:
    init_db()
    yield


app = FastAPI(
    title="PCR Primer Dimer Screening API",
    version="1.0.0",
    description=(
        "Screens named IUPAC DNA primers for self-dimers and pairwise "
        "heterodimers by enumerating every antiparallel alignment offset, "
        "scoring consecutive and 3'-end complementarity, and reporting the "
        "most dangerous alignment per pair."
    ),
    lifespan=lifespan,
)


@app.get("/health", tags=["meta"])
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post(
    "/api/v1/batches",
    status_code=201,
    response_model=BatchResponse,
    tags=["screening"],
)
def create_batch(request: BatchRequest) -> BatchResponse:
    """Screen a batch of primers and persist the results.

    Every primer is checked against itself (self-dimer) and against every
    other primer (heterodimer).  The full per-pair report is returned and
    stored for later retrieval via ``GET /api/v1/batches/{batch_id}``.
    """
    pairs = [(p.name, p.sequence) for p in request.primers]
    results = screen_primers(
        pairs,
        request.max_consecutive_threshold,
        request.three_prime_threshold,
    )
    created_at = datetime.now(timezone.utc).isoformat()
    with get_conn() as conn:
        cursor = conn.execute(
            "INSERT INTO batches (created_at, max_consecutive_threshold,"
            " three_prime_threshold, primer_count) VALUES (?, ?, ?, ?)",
            (
                created_at,
                request.max_consecutive_threshold,
                request.three_prime_threshold,
                len(request.primers),
            ),
        )
        batch_id = int(cursor.lastrowid)
        conn.executemany(
            "INSERT INTO primers (batch_id, name, sequence) VALUES (?, ?, ?)",
            [(batch_id, p.name, p.sequence) for p in request.primers],
        )
        conn.executemany(
            "INSERT INTO pair_results (batch_id, primer_a, primer_b, kind,"
            " offset, max_consecutive, three_prime_consecutive, flagged,"
            " reasons) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [
                (
                    batch_id,
                    r["primer_a"],
                    r["primer_b"],
                    r["kind"],
                    r["offset"],
                    r["max_consecutive"],
                    r["three_prime_consecutive"],
                    int(r["flagged"]),
                    json.dumps(r["reasons"]),
                )
                for r in results
            ],
        )
    return BatchResponse(
        batch_id=batch_id,
        created_at=created_at,
        max_consecutive_threshold=request.max_consecutive_threshold,
        three_prime_threshold=request.three_prime_threshold,
        primer_count=len(request.primers),
        pair_count=len(results),
        flagged_count=sum(1 for r in results if r["flagged"]),
        primers=request.primers,
        results=[PairResult(**r) for r in results],
    )


@app.get(
    "/api/v1/batches",
    response_model=list[BatchSummary],
    tags=["screening"],
)
def list_batches() -> list[BatchSummary]:
    """List all submitted batches, newest first."""
    with get_conn() as conn:
        rows = conn.execute(
            "SELECT b.id, b.created_at, b.primer_count,"
            " COUNT(r.id) AS pair_count,"
            " COALESCE(SUM(r.flagged), 0) AS flagged_count"
            " FROM batches b"
            " LEFT JOIN pair_results r ON r.batch_id = b.id"
            " GROUP BY b.id"
            " ORDER BY b.id DESC"
        ).fetchall()
    return [
        BatchSummary(
            batch_id=row["id"],
            created_at=row["created_at"],
            primer_count=row["primer_count"],
            pair_count=row["pair_count"],
            flagged_count=row["flagged_count"],
        )
        for row in rows
    ]


@app.get(
    "/api/v1/batches/{batch_id}",
    response_model=BatchResponse,
    tags=["screening"],
)
def get_batch(batch_id: int) -> BatchResponse:
    """Retrieve the stored screening report for a batch."""
    with get_conn() as conn:
        batch = conn.execute(
            "SELECT * FROM batches WHERE id = ?", (batch_id,)
        ).fetchone()
        if batch is None:
            raise HTTPException(status_code=404, detail="batch not found")
        primer_rows = conn.execute(
            "SELECT name, sequence FROM primers WHERE batch_id = ? ORDER BY id",
            (batch_id,),
        ).fetchall()
        result_rows = conn.execute(
            "SELECT * FROM pair_results WHERE batch_id = ? ORDER BY id",
            (batch_id,),
        ).fetchall()
    results = [
        PairResult(
            primer_a=row["primer_a"],
            primer_b=row["primer_b"],
            kind=row["kind"],
            offset=row["offset"],
            max_consecutive=row["max_consecutive"],
            three_prime_consecutive=row["three_prime_consecutive"],
            flagged=bool(row["flagged"]),
            reasons=json.loads(row["reasons"]),
        )
        for row in result_rows
    ]
    return BatchResponse(
        batch_id=batch["id"],
        created_at=batch["created_at"],
        max_consecutive_threshold=batch["max_consecutive_threshold"],
        three_prime_threshold=batch["three_prime_threshold"],
        primer_count=batch["primer_count"],
        pair_count=len(results),
        flagged_count=sum(1 for r in results if r.flagged),
        primers=[PrimerInput(**dict(row)) for row in primer_rows],
        results=results,
    )
