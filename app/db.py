"""SQLite persistence layer.

The database file location is read from the ``DATABASE_PATH`` environment
variable at call time (default ``primer_dimer.db`` in the working
directory).  The schema is created idempotently at application startup.
"""
from __future__ import annotations

import os
import sqlite3
from collections.abc import Iterator
from contextlib import contextmanager

SCHEMA = """
CREATE TABLE IF NOT EXISTS batches (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    created_at TEXT NOT NULL,
    max_consecutive_threshold INTEGER NOT NULL,
    three_prime_threshold INTEGER NOT NULL,
    primer_count INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS primers (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    batch_id INTEGER NOT NULL REFERENCES batches(id),
    name TEXT NOT NULL,
    sequence TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS pair_results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    batch_id INTEGER NOT NULL REFERENCES batches(id),
    primer_a TEXT NOT NULL,
    primer_b TEXT NOT NULL,
    kind TEXT NOT NULL,
    offset INTEGER NOT NULL,
    max_consecutive INTEGER NOT NULL,
    three_prime_consecutive INTEGER NOT NULL,
    flagged INTEGER NOT NULL,
    reasons TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_pair_results_batch ON pair_results(batch_id);
CREATE INDEX IF NOT EXISTS idx_primers_batch ON primers(batch_id);
"""


def db_path() -> str:
    return os.environ.get("DATABASE_PATH", "primer_dimer.db")


@contextmanager
def get_conn() -> Iterator[sqlite3.Connection]:
    conn = sqlite3.connect(db_path())
    conn.row_factory = sqlite3.Row
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()


def init_db() -> None:
    """Create the database file (and parent directory) and schema."""
    path = db_path()
    if path != ":memory:":
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with get_conn() as conn:
        conn.executescript(SCHEMA)
