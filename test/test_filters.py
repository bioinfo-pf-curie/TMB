# coding: utf-8
"""
Tests for pytmb.filters (subsetINFO, infoTag2dl, info2dl,
isAnnotatedAs, isPolym, isCancerHotspot).
"""

import pytest

from pytmb.filters import (
    subsetINFO,
    infoTag2dl,
    info2dl,
    isAnnotatedAs,
    isPolym,
    isCancerHotspot,
)

# ---------------------------------------------------------------------------
# subsetINFO
# ---------------------------------------------------------------------------


class TestSubsetINFO:
    def test_subset_dict_keeps_present_keys(self):
        annot = {"a": 1, "b": 2, "c": 3}
        result = subsetINFO(annot, ["a", "c"])
        assert result == {"a": 1, "c": 3}

    def test_subset_dict_ignores_missing_keys(self):
        annot = {"a": 1}
        result = subsetINFO(annot, ["a", "z"])
        assert result == {"a": 1}

    def test_subset_dict_all_missing_returns_empty(self):
        annot = {"a": 1}
        result = subsetINFO(annot, ["x", "y"])
        assert result == {}

    def test_subset_list_of_dicts(self):
        annot = [{"a": 1, "b": 2}, {"b": 3, "c": 4}]
        result = subsetINFO(annot, ["b"])
        assert result == [{"b": 2}, {"b": 3}]

    def test_subset_list_drops_empty_dicts(self):
        annot = [{"a": 1}, {"b": 2}]
        result = subsetINFO(annot, ["b"])
        # dict for first element has no "b" → dropped
        assert result == [{"b": 2}]


# ---------------------------------------------------------------------------
# infoTag2dl
# ---------------------------------------------------------------------------


class TestInfoTag2dl:
    def test_none_returns_none(self):
        assert infoTag2dl(None) is None

    def test_single_annotation(self):
        result = infoTag2dl("A|B|C")
        assert result == [{0: "A", 1: "B", 2: "C"}]

    def test_multiple_annotations(self):
        result = infoTag2dl("A|B,C|D")
        assert len(result) == 2
        assert result[0] == {0: "A", 1: "B"}
        assert result[1] == {0: "C", 1: "D"}


# ---------------------------------------------------------------------------
# info2dl
# ---------------------------------------------------------------------------


class TestInfo2dl:
    def test_none_returns_none(self):
        assert info2dl(None) is None

    def test_dict_wrapped_in_list(self):
        d = {"key": "val"}
        result = info2dl(d)
        assert result == [{"key": "val"}]

    def test_list_length_is_one(self):
        result = info2dl({"a": 1, "b": 2})
        assert len(result) == 1


# ---------------------------------------------------------------------------
# isAnnotatedAs
# ---------------------------------------------------------------------------


class TestIsAnnotatedAs:
    """isAnnotatedAs(v, infos, flags, sep) – v is unused."""

    def _infos(self, raw):
        return infoTag2dl(raw)

    def test_match_returns_true(self):
        # Simulates snpEff-style ANN: positional index 1 = effect
        infos = [{0: "gene1", 1: "exonic&splicing", 2: "HIGH"}]
        flags = {1: ["exonic"]}
        assert isAnnotatedAs(None, infos, flags, sep="&") is True

    def test_no_match_returns_false(self):
        infos = [{0: "gene1", 1: "intronic", 2: "LOW"}]
        flags = {1: ["exonic"]}
        assert isAnnotatedAs(None, infos, flags, sep="&") is False

    def test_multiple_flags_any_match(self):
        infos = [{0: "gene1", 1: "splicing", 2: "MODERATE"}]
        flags = {1: ["exonic", "splicing"]}
        assert isAnnotatedAs(None, infos, flags, sep="&") is True

    def test_multiple_annots_first_matches(self):
        infos = [
            {0: "gene1", 1: "intronic", 2: "LOW"},
            {0: "gene2", 1: "exonic", 2: "HIGH"},
        ]
        flags = {1: ["exonic"]}
        assert isAnnotatedAs(None, infos, flags, sep="&") is True


# ---------------------------------------------------------------------------
# isPolym
# ---------------------------------------------------------------------------


class TestIsPolym:
    def test_above_threshold_returns_true(self):
        infos = {"gnomAD_exome_ALL": 0.05}
        assert isPolym(None, infos, ["gnomAD_exome_ALL"], 0.001) is True

    def test_below_threshold_returns_false(self):
        infos = {"gnomAD_exome_ALL": 0.0001}
        assert isPolym(None, infos, ["gnomAD_exome_ALL"], 0.001) is False

    def test_missing_key_returns_false(self):
        infos = {}
        assert isPolym(None, infos, ["gnomAD_exome_ALL"], 0.001) is False

    def test_dot_value_returns_false(self):
        infos = {"gnomAD_exome_ALL": "."}
        assert isPolym(None, infos, ["gnomAD_exome_ALL"], 0.001) is False

    def test_none_value_returns_false(self):
        infos = {"gnomAD_exome_ALL": None}
        assert isPolym(None, infos, ["gnomAD_exome_ALL"], 0.001) is False

    def test_tuple_any_above_threshold(self):
        infos = {"gnomAD_exome_ALL": (None, 0.02, 0.0001)}
        assert isPolym(None, infos, ["gnomAD_exome_ALL"], 0.01) is True

    def test_tuple_all_below_threshold(self):
        infos = {"gnomAD_exome_ALL": (0.0001, 0.0002)}
        assert isPolym(None, infos, ["gnomAD_exome_ALL"], 0.01) is False


# ---------------------------------------------------------------------------
# isCancerHotspot
# ---------------------------------------------------------------------------


class TestIsCancerHotspot:
    def test_present_and_non_empty_returns_true(self):
        infos = {"cosmic86": "COSM12345"}
        assert isCancerHotspot(None, infos, ["cosmic86"]) is True

    def test_dot_value_returns_false(self):
        infos = {"cosmic86": "."}
        assert isCancerHotspot(None, infos, ["cosmic86"]) is False

    def test_none_value_returns_false(self):
        infos = {"cosmic86": None}
        assert isCancerHotspot(None, infos, ["cosmic86"]) is False

    def test_missing_key_returns_false(self):
        infos = {}
        assert isCancerHotspot(None, infos, ["cosmic86"]) is False
