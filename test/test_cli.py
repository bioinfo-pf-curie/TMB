# coding: utf-8
"""
Tests for the pytmb.cli entry points (pyTMB, pyEffGenomeSize).
"""

import sys
import pytest

from pytmb import __version__
from pytmb.cli import run_effgenomesize, run_tmb


class TestVersion:
    @pytest.mark.parametrize("cli, prog", [
        (run_tmb, "pyTMB"),
        (run_effgenomesize, "pyEffGenomeSize"),
    ])
    def test_version_output(self, cli, prog, monkeypatch, capsys):
        """--version prints the program name and the version, no parentheses."""
        monkeypatch.setattr(sys, "argv", [prog, "--version"])
        with pytest.raises(SystemExit) as exc:
            cli.argsParse()
        assert exc.value.code == 0
        assert capsys.readouterr().out == f"{prog} {__version__}\n"
