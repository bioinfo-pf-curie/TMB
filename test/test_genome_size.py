# coding: utf-8
"""
Tests for pytmb.genome_size (getEffGenomeSizeFromBed,
getEffGenomeSizeFromMosdepth).
"""

import pytest

from pytmb.genome_size import getEffGenomeSizeFromBed, getEffGenomeSizeFromMosdepth


# ---------------------------------------------------------------------------
# getEffGenomeSizeFromBed
# ---------------------------------------------------------------------------


class TestGetEffGenomeSizeFromBed:
    def test_basic_intervals(self, tmp_bed):
        """Sum of (end - start) across all intervals."""
        bed = tmp_bed(["chr1\t0\t100", "chr1\t200\t350"])
        assert getEffGenomeSizeFromBed(bed) == 250  # 100 + 150

    def test_single_interval(self, tmp_bed):
        bed = tmp_bed(["chr1\t0\t1000"])
        assert getEffGenomeSizeFromBed(bed) == 1000

    def test_skips_comment_lines(self, tmp_bed):
        bed = tmp_bed(["# header", "chr1\t0\t500"])
        assert getEffGenomeSizeFromBed(bed) == 500

    def test_skips_empty_lines(self, tmp_bed):
        bed = tmp_bed(["chr1\t0\t300", "", "chr2\t0\t200"])
        assert getEffGenomeSizeFromBed(bed) == 500

    def test_extra_columns_ignored(self, tmp_bed):
        """BED files may have > 3 columns; only cols 1-2 (start/end) matter."""
        bed = tmp_bed(["chr1\t0\t200\tsome_gene\t0\t+"])
        assert getEffGenomeSizeFromBed(bed) == 200

    def test_file_not_found_raises_systemexit(self, tmp_path):
        with pytest.raises(SystemExit):
            getEffGenomeSizeFromBed(str(tmp_path / "missing.bed"))

    def test_malformed_line_raises_systemexit(self, tmp_bed):
        """A line with fewer than 3 columns triggers sys.exit."""
        bed = tmp_bed(["chr1\t100"])  # only 2 columns
        with pytest.raises(SystemExit):
            getEffGenomeSizeFromBed(bed)

    def test_non_integer_coordinates_raise_systemexit(self, tmp_bed):
        bed = tmp_bed(["chr1\tstart\tend"])
        with pytest.raises(SystemExit):
            getEffGenomeSizeFromBed(bed)

    def test_real_panel_bed(self, panel_bed):
        """Panel BED file returns a positive integer."""
        size = getEffGenomeSizeFromBed(panel_bed)
        assert isinstance(size, int)
        assert size > 0


# ---------------------------------------------------------------------------
# getEffGenomeSizeFromMosdepth (plain BED mode)
# ---------------------------------------------------------------------------


class TestGetEffGenomeSizeFromMosdepthPlain:
    def test_plain_bed_returns_sum(self, tmp_bed):
        bed = tmp_bed(["chr1\t0\t100", "chr2\t50\t200"])
        assert getEffGenomeSizeFromMosdepth(bed, use_mosdepth=False) == 250

    def test_verbose_does_not_crash(self, tmp_bed, capsys):
        bed = tmp_bed(["chr1\t0\t100"])
        result = getEffGenomeSizeFromMosdepth(bed, use_mosdepth=False, verbose=True)
        assert result == 100
        captured = capsys.readouterr()
        assert "Total region size" in captured.out


# ---------------------------------------------------------------------------
# getEffGenomeSizeFromMosdepth (mosdepth thresholds mode)
# ---------------------------------------------------------------------------


class TestGetEffGenomeSizeFromMosdepthMosdepth:
    def _make_mosdepth_bed(self, tmp_path, lines):
        p = tmp_path / "mosdepth.bed"
        header = "#chrom\tstart\tend\tregion\tcoverage"
        p.write_text(header + "\n" + "\n".join(lines) + "\n")
        return str(p)

    def test_returns_sum_of_coverage_column(self, tmp_path):
        bed = self._make_mosdepth_bed(
            tmp_path,
            ["chr1\t0\t100\tregion1\t80", "chr1\t100\t200\tregion2\t90"],
        )
        assert getEffGenomeSizeFromMosdepth(bed, use_mosdepth=True) == 170

    def test_verbose_shows_callable_region(self, tmp_path, capsys):
        bed = self._make_mosdepth_bed(
            tmp_path,
            ["chr1\t0\t200\tregion1\t150"],
        )
        getEffGenomeSizeFromMosdepth(bed, use_mosdepth=True, verbose=True)
        captured = capsys.readouterr()
        assert "Callable region" in captured.out
