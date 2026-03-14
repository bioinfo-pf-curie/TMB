# coding: utf-8
"""
Shared pytest fixtures for the pyTMB test suite.
"""

import os
import tempfile
import pytest

# ---------------------------------------------------------------------------
# Paths (relative to repo root, resolved at collection time)
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


@pytest.fixture(scope="session")
def repo_root():
    """Absolute path to the repository root."""
    return REPO_ROOT


@pytest.fixture(scope="session")
def config_dir(repo_root):
    """Path to the config/ directory."""
    return os.path.join(repo_root, "config")


@pytest.fixture(scope="session")
def wes_dir(repo_root):
    """Path to the test/wes/ directory."""
    return os.path.join(repo_root, "test", "wes")


@pytest.fixture(scope="session")
def panel_dir(repo_root):
    """Path to the test/panel/ directory."""
    return os.path.join(repo_root, "test", "panel")


@pytest.fixture(scope="session")
def annovar_config(config_dir):
    """Path to the ANNOVAR 102015 YAML config."""
    return os.path.join(config_dir, "annovar_102015.yml")


@pytest.fixture(scope="session")
def mutect2_config(config_dir):
    """Path to the Mutect2 caller YAML config."""
    return os.path.join(config_dir, "mutect2.yml")


@pytest.fixture(scope="session")
def snpeff_config(config_dir):
    """Path to the snpEff YAML config."""
    return os.path.join(config_dir, "snpeff.yml")


@pytest.fixture(scope="session")
def vcf1(wes_dir):
    """Path to the first test WES VCF (Mutect2 + snpEff)."""
    return os.path.join(wes_dir, "test1_mutect2_snpeff.vcf.gz")


@pytest.fixture(scope="session")
def vcf2(wes_dir):
    """Path to the second test WES VCF (Mutect2 + snpEff)."""
    return os.path.join(wes_dir, "test2_mutect2_snpeff.vcf.gz")


@pytest.fixture(scope="session")
def panel_bed(panel_dir):
    """Path to the panel BED file."""
    return os.path.join(panel_dir, "dragon_design_v2.bed")


@pytest.fixture()
def tmp_bed(tmp_path):
    """Return a helper that writes a small BED file and returns its path."""
    def _make_bed(lines):
        p = tmp_path / "test.bed"
        p.write_text("\n".join(lines) + "\n")
        return str(p)
    return _make_bed
