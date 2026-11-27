import sys
from pathlib import Path
import pytest
from phacs_nmr import cli, gui


class TestCLI:
    def test_parse_args_defaults(self):
        args = cli.parse_args([])
        assert args.path is None
        assert args.noise_level == 10
        assert args.opt_noise is None
        assert args.gui is False

    def test_parse_args_single_file_and_noise(self):
        args1 = cli.parse_args(["sample.txt"])
        assert args1.path == "sample.txt"
        assert args1.noise_level == 10

        args2 = cli.parse_args(["sample.txt", "5"])
        assert args2.path == "sample.txt"
        assert args2.noise_level == 5

        args3 = cli.parse_args(["sample.txt", "-n", "5"])
        assert args3.opt_noise == 5

    def test_parse_args_gui_flag(self):
        args = cli.parse_args(["--gui"])
        assert args.gui is True

    def test_parse_args_invalid_noise(self):
        with pytest.raises(SystemExit):
            cli.parse_args(["sample.txt", "7"])

    def test_main_file_not_found(self, tmp_path):
        missing = tmp_path / "does_not_exist.txt"
        with pytest.raises(FileNotFoundError, match="Input path not found"):
            cli.main([str(missing)])

    def test_main_triggers_gui_when_no_args_or_gui_flag(self, monkeypatch):
        called = []

        def mock_gui_main():
            called.append(True)
            return "gui_launched"

        monkeypatch.setattr(gui, "main", mock_gui_main)

        # 1. No arguments
        res1 = cli.main([])
        assert res1 == "gui_launched"
        assert len(called) == 1

        # 2. Explicit --gui flag
        res2 = cli.main(["--gui"])
        assert res2 == "gui_launched"
        assert len(called) == 2

    def test_main_executes_single_file(self, monkeypatch, tmp_path):
        test_file = tmp_path / "spectrum.txt"
        test_file.write_text("1, 10.0, 0.0, 100.0\n", encoding="cp1252")

        called_with = []

        def mock_start_analysis(file_path, noise):
            called_with.append((file_path, noise))
            return ["mock_figure"]

        monkeypatch.setattr(cli, "start_analysis", mock_start_analysis)

        # 1. Default noise (10)
        res1 = cli.main([str(test_file)])
        assert res1 == ["mock_figure"]
        assert called_with[-1] == (str(test_file.resolve()), 10)

        # 2. Positional noise (5)
        res2 = cli.main([str(test_file), "5"])
        assert res2 == ["mock_figure"]
        assert called_with[-1] == (str(test_file.resolve()), 5)

        # 3. Flag noise override (-n 5)
        res3 = cli.main([str(test_file), "-n", "5"])
        assert res3 == ["mock_figure"]
        assert called_with[-1] == (str(test_file.resolve()), 5)

    def test_main_executes_directory_batch(self, monkeypatch, tmp_path, capsys):
        test_dir = tmp_path / "batch_spectra"
        test_dir.mkdir()
        (test_dir / "s1.txt").write_text("1, 10.0, 0.0, 100.0\n", encoding="cp1252")
        (test_dir / "s2.txt").write_text("2, 20.0, 0.0, 90.0\n", encoding="cp1252")

        processed = []

        def mock_start_analysis(file_path, noise):
            processed.append((file_path, noise))
            return ["fig"]

        monkeypatch.setattr(cli, "start_analysis", mock_start_analysis)

        results = cli.main([str(test_dir), "5"])
        assert len(results) == 2
        assert len(processed) == 2
        captured = capsys.readouterr()
        assert "=== Processing s1.txt ===" in captured.out
        assert "=== Processing s2.txt ===" in captured.out
