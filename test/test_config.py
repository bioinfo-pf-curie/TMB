# coding: utf-8
"""
Tests for pytmb.config (loadConfig).
"""

import os
import pytest
import yaml

from pytmb.config import loadConfig


class TestLoadConfig:
    def test_load_annovar_config(self, annovar_config):
        """loadConfig returns a dict for a valid YAML file."""
        cfg = loadConfig(annovar_config)
        assert isinstance(cfg, dict)

    def test_annovar_config_has_required_keys(self, annovar_config):
        """The ANNOVAR config has all expected top-level keys."""
        cfg = loadConfig(annovar_config)
        for key in ("isCoding", "isSplicing", "isNonCoding", "isSynonymous",
                    "isNonSynonymous", "cancerDb", "polymDb"):
            assert key in cfg, f"Missing key: {key}"

    def test_load_mutect2_config(self, mutect2_config):
        """Mutect2 caller config loads and contains freq/depth/altDepth."""
        cfg = loadConfig(mutect2_config)
        assert "freq" in cfg
        assert "depth" in cfg
        assert "altDepth" in cfg

    def test_load_nonexistent_file_raises(self, tmp_path):
        """loadConfig raises FileNotFoundError for a missing file."""
        with pytest.raises(FileNotFoundError):
            loadConfig(str(tmp_path / "does_not_exist.yml"))

    def test_load_invalid_yaml_raises(self, tmp_path):
        """loadConfig raises an exception for malformed YAML."""
        bad = tmp_path / "bad.yml"
        bad.write_text("key: [unclosed bracket\n")
        with pytest.raises(Exception):
            loadConfig(str(bad))

    def test_roundtrip(self, tmp_path):
        """A dict written to YAML and read back via loadConfig is identical."""
        data = {"a": 1, "b": [2, 3], "c": {"nested": True}}
        p = tmp_path / "roundtrip.yml"
        p.write_text(yaml.dump(data))
        assert loadConfig(str(p)) == data
