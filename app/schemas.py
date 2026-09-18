"""Pydantic request/response models for the API."""
from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, field_validator

from app.dimer import IUPAC_BASES


class PrimerInput(BaseModel):
    name: str = Field(min_length=1, max_length=128, examples=["fwd_01"])
    sequence: str = Field(min_length=1, max_length=256, examples=["ATGCRTAGCT"])

    @field_validator("sequence")
    @classmethod
    def normalize_sequence(cls, value: str) -> str:
        value = value.strip().upper()
        invalid = sorted(set(value) - set(IUPAC_BASES))
        if invalid:
            raise ValueError(
                f"invalid IUPAC base symbols: {', '.join(invalid)}"
            )
        return value


class BatchRequest(BaseModel):
    primers: list[PrimerInput] = Field(min_length=1, max_length=500)
    max_consecutive_threshold: int = Field(
        default=4,
        ge=1,
        le=100,
        description="Flag pairs whose longest consecutive complementarity reaches this length.",
    )
    three_prime_threshold: int = Field(
        default=3,
        ge=1,
        le=100,
        description="Flag pairs whose 3'-end consecutive complementarity reaches this length.",
    )

    @field_validator("primers")
    @classmethod
    def unique_names(cls, value: list[PrimerInput]) -> list[PrimerInput]:
        names = [p.name for p in value]
        duplicates = sorted({n for n in names if names.count(n) > 1})
        if duplicates:
            raise ValueError(f"duplicate primer names: {', '.join(duplicates)}")
        return value


class PairResult(BaseModel):
    primer_a: str
    primer_b: str
    kind: Literal["self", "hetero"]
    offset: int = Field(
        description="Alignment offset of the most dangerous alignment: index in "
        "reversed primer B that pairs with primer A[0]. Ties are broken by the "
        "smallest offset."
    )
    max_consecutive: int = Field(
        description="Longest run of consecutive complementary bases in that alignment."
    )
    three_prime_consecutive: int = Field(
        description="Longest complementary run involving either primer's aligned 3' end."
    )
    flagged: bool
    reasons: list[str]


class BatchResponse(BaseModel):
    batch_id: int
    created_at: str
    max_consecutive_threshold: int
    three_prime_threshold: int
    primer_count: int
    pair_count: int
    flagged_count: int
    primers: list[PrimerInput]
    results: list[PairResult]


class BatchSummary(BaseModel):
    batch_id: int
    created_at: str
    primer_count: int
    pair_count: int
    flagged_count: int
