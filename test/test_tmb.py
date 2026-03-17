# coding: utf-8
"""
Integration tests for pytmb.tmb.calculate_tmb using the real WES test VCFs.
"""

import pytest

from pytmb.config import loadConfig
from pytmb.tmb import calculate_tmb

EFF_GENOME_SIZE = 33_280_000  # same value used in test/wes/test.sh


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _run_tmb(vcf_path, annovar_config, mutect2_config, **kwargs):
    """Convenience wrapper around calculate_tmb with sensible defaults."""
    db_flags = loadConfig(annovar_config)
    caller_flags = loadConfig(mutect2_config)
    return calculate_tmb(
        vcf_path=vcf_path,
        db_flags=db_flags,
        caller_flags=caller_flags,
        eff_genome_size=EFF_GENOME_SIZE,
        vaf=0.05,
        min_depth=20,
        min_alt_depth=2,
        filter_low_qual=True,
        filter_non_coding=True,
        filter_syn=True,
        filter_polym=True,
        polym_db="1k,gnomad",
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Result structure
# ---------------------------------------------------------------------------


class TestCalculateTmbResultStructure:
    def test_result_is_dict(self, vcf1, annovar_config, mutect2_config):
        result = _run_tmb(vcf1, annovar_config, mutect2_config)
        assert isinstance(result, dict)

    def test_result_has_required_keys(self, vcf1, annovar_config, mutect2_config):
        result = _run_tmb(vcf1, annovar_config, mutect2_config)
        for key in ("tmb", "var_counter", "var_ni", "var_tmb",
                    "eff_genome_size", "filter_stats"):
            assert key in result, f"Missing key: {key}"

    def test_tmb_is_non_negative_float(self, vcf1, annovar_config, mutect2_config):
        result = _run_tmb(vcf1, annovar_config, mutect2_config)
        assert isinstance(result["tmb"], float)
        assert result["tmb"] >= 0.0

    def test_eff_genome_size_matches_input(self, vcf1, annovar_config, mutect2_config):
        result = _run_tmb(vcf1, annovar_config, mutect2_config)
        assert result["eff_genome_size"] == EFF_GENOME_SIZE

    def test_var_counters_are_non_negative_ints(self, vcf1, annovar_config, mutect2_config):
        result = _run_tmb(vcf1, annovar_config, mutect2_config)
        for key in ("var_counter", "var_ni", "var_tmb"):
            assert isinstance(result[key], int), f"{key} is not int"
            assert result[key] >= 0, f"{key} is negative"

    def test_var_tmb_leq_var_counter(self, vcf1, annovar_config, mutect2_config):
        result = _run_tmb(vcf1, annovar_config, mutect2_config)
        assert result["var_tmb"] <= result["var_counter"]

    def test_filter_stats_is_dict(self, vcf1, annovar_config, mutect2_config):
        result = _run_tmb(vcf1, annovar_config, mutect2_config)
        assert isinstance(result["filter_stats"], dict)


# ---------------------------------------------------------------------------
# Both VCF files are processable
# ---------------------------------------------------------------------------


class TestCalculateTmbBothVcfs:
    @pytest.mark.parametrize("vcf_fixture", ["vcf1", "vcf2"])
    def test_vcf_runs_without_error(self, vcf_fixture, request,
                                    annovar_config, mutect2_config):
        vcf_path = request.getfixturevalue(vcf_fixture)
        result = _run_tmb(vcf_path, annovar_config, mutect2_config)
        assert result["var_counter"] > 0

    def test_tmb_formula(self, vcf1, annovar_config, mutect2_config):
        """TMB == round(var_tmb / (eff_genome_size / 1e6), 2)."""
        result = _run_tmb(vcf1, annovar_config, mutect2_config)
        expected = round(result["var_tmb"] / (EFF_GENOME_SIZE / 1e6), 2)
        assert result["tmb"] == expected


# ---------------------------------------------------------------------------
# Filter behaviour
# ---------------------------------------------------------------------------


class TestCalculateTmbFilters:
    def test_no_filters_gives_higher_count(self, vcf1, annovar_config, mutect2_config):
        """Disabling all filters should yield ≥ as many TMB variants as
        with strict filtering enabled."""
        db_flags = loadConfig(annovar_config)
        caller_flags = loadConfig(mutect2_config)
        result_lenient = calculate_tmb(
            vcf_path=vcf1,
            db_flags=db_flags,
            caller_flags=caller_flags,
            eff_genome_size=EFF_GENOME_SIZE,
        )
        result_strict = _run_tmb(vcf1, annovar_config, mutect2_config)
        assert result_lenient["var_tmb"] >= result_strict["var_tmb"]

    def test_filter_indels_reduces_count(self, vcf1, annovar_config, mutect2_config):
        """Enabling indel filtering should yield ≤ TMB variants vs no indel filter."""
        db_flags = loadConfig(annovar_config)
        caller_flags = loadConfig(mutect2_config)

        no_indel_filter = calculate_tmb(
            vcf_path=vcf1,
            db_flags=db_flags,
            caller_flags=caller_flags,
            eff_genome_size=EFF_GENOME_SIZE,
            filter_indels=False,
        )
        with_indel_filter = calculate_tmb(
            vcf_path=vcf1,
            db_flags=db_flags,
            caller_flags=caller_flags,
            eff_genome_size=EFF_GENOME_SIZE,
            filter_indels=True,
        )
        assert with_indel_filter["var_tmb"] <= no_indel_filter["var_tmb"]
