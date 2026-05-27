"""Unit tests for biomarker enrichment logic."""
import pandas as pd
import pytest

from pipeline.stages.biomarker.biomarker import _compute_enrichment, _STRONGLY_DEP


def _make_dep_records(effects: list[float], lineage: str = "CNS/Brain") -> list[dict]:
    return [
        {
            "ModelID": f"ACH-{i:06d}",
            "gene_effect": e,
            "OncotreeLineage": lineage,
            "OncotreePrimaryDisease": "Test",
            "OncotreeSubtype": "Test",
            "CCLEName": f"Cell{i}",
        }
        for i, e in enumerate(effects)
    ]


def _make_mutations(model_ids: list[str], gene: str = "TP53") -> pd.DataFrame:
    return pd.DataFrame({
        "ModelID": model_ids,
        "Hugo_Symbol": gene,
        "Variant_Classification": "Missense_Mutation",
        "isDeleterious": True,
    })


# ── enrichment direction ───────────────────────────────────────────────────────

class TestComputeEnrichment:
    def test_gene_enriched_in_dependent_lines(self) -> None:
        # 6 strongly dependent, 4 non-dependent; mutation in 5 of 6 dependent only
        effects = [_STRONGLY_DEP - 0.1] * 6 + [0.0] * 4
        dep_records = _make_dep_records(effects)
        # Mutate first 5 of the 6 dependent lines
        dep_model_ids = [r["ModelID"] for r in dep_records[:5]]
        mutations = _make_mutations(dep_model_ids)

        results = _compute_enrichment(dep_records, mutations, "CNS/Brain", min_lines=3)

        assert len(results) == 1
        row = results.iloc[0]
        assert row["gene"] == "TP53"
        assert row["odds_ratio"] > 1.0
        assert row["enrichment_direction"] == "enriched_in_dependent"

    def test_gene_below_min_lines_excluded(self) -> None:
        # Only 2 dependent lines with mutation — below min_lines=3
        effects = [_STRONGLY_DEP - 0.1] * 6 + [0.0] * 4
        dep_records = _make_dep_records(effects)
        dep_model_ids = [r["ModelID"] for r in dep_records[:2]]
        mutations = _make_mutations(dep_model_ids)

        results = _compute_enrichment(dep_records, mutations, "CNS/Brain", min_lines=3)
        assert len(results) == 0

    def test_gene_exactly_at_min_lines_included(self) -> None:
        # Exactly min_lines=3 dependent lines with mutation
        effects = [_STRONGLY_DEP - 0.1] * 6 + [0.0] * 4
        dep_records = _make_dep_records(effects)
        dep_model_ids = [r["ModelID"] for r in dep_records[:3]]
        mutations = _make_mutations(dep_model_ids)

        results = _compute_enrichment(dep_records, mutations, "CNS/Brain", min_lines=3)
        assert len(results) == 1
        assert results.iloc[0]["n_dependent_with_mut"] == 3

    def test_no_genes_pass_filter_returns_empty_df(self) -> None:
        effects = [_STRONGLY_DEP - 0.1] * 6 + [0.0] * 4
        dep_records = _make_dep_records(effects)
        # Mutation only in 1 dependent line — below min_lines=3
        mutations = _make_mutations([dep_records[0]["ModelID"]])

        results = _compute_enrichment(dep_records, mutations, "CNS/Brain", min_lines=3)
        assert len(results) == 0
        assert "gene" in results.columns

    def test_invalid_lineage_raises_value_error(self) -> None:
        dep_records = _make_dep_records([_STRONGLY_DEP - 0.1] * 3 + [0.0] * 3)
        mutations = _make_mutations(["ACH-000000"])

        with pytest.raises(ValueError, match="not found in DepMap cache"):
            _compute_enrichment(dep_records, mutations, "InvalidLineage", min_lines=3)

    def test_results_sorted_by_p_value_ascending(self) -> None:
        # Two genes: one with strong enrichment, one with weaker
        effects = [_STRONGLY_DEP - 0.1] * 8 + [0.0] * 2
        dep_records = _make_dep_records(effects)
        dep_ids = [r["ModelID"] for r in dep_records[:8]]
        nondep_ids = [r["ModelID"] for r in dep_records[8:]]

        # TP53: enriched in all 8 dependent — very low p
        muts_tp53 = _make_mutations(dep_ids, gene="TP53")
        # CDKN2A: enriched in 3 dependent + 1 non-dependent — weaker signal
        muts_cdkn2a = _make_mutations(dep_ids[:3] + nondep_ids[:1], gene="CDKN2A")
        mutations = pd.concat([muts_tp53, muts_cdkn2a], ignore_index=True)

        results = _compute_enrichment(dep_records, mutations, "CNS/Brain", min_lines=3)
        assert len(results) >= 2
        p_values = results["p_value"].tolist()
        assert p_values == sorted(p_values), "Results not sorted by p_value ascending"

    def test_odds_ratio_and_counts_correct(self) -> None:
        # 4 strongly dep, 4 non-dep; mutation in all 4 dep, 0 non-dep
        effects = [_STRONGLY_DEP - 0.1] * 4 + [0.0] * 4
        dep_records = _make_dep_records(effects)
        dep_ids = [r["ModelID"] for r in dep_records[:4]]
        mutations = _make_mutations(dep_ids)

        results = _compute_enrichment(dep_records, mutations, "CNS/Brain", min_lines=3)
        row = results.iloc[0]
        assert row["n_dependent_with_mut"] == 4
        assert row["n_nondependent_with_mut"] == 0
        assert row["odds_ratio"] == float("inf") or row["odds_ratio"] > 10
