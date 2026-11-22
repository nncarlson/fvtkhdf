"""Exercise diagnostic alternatives without relying on MPI abort timing."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class AssertionCheckerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.source = self.directory / "fixture.F90"
        self.source.write_text("! write guard\n! registration guard\n! unrelated guard\n")

    def run_checker(self, line, expected, *, status=1, logfile="", markers=None):
        command = [sys.executable, str(Path(__file__).with_name("vtkhdf_expect_assertion.py")),
                   "--source", str(self.source)]
        for marker in markers or ("write guard", "registration guard"):
            command.extend(("--assertion", marker))
        command.extend(("--diagnostic-glob", "interface-assertion-R*.log", "--",
                        sys.executable, "-c",
                        "from pathlib import Path; import sys; "
                        "text = sys.argv[1]; "
                        "Path(sys.argv[3]).write_text(text) if sys.argv[3] else print(text); "
                        "sys.exit(int(sys.argv[2]))",
                        f"Panic! Assertion failed at {self.source}:{line}" if line else "",
                        str(status), logfile))
        result = subprocess.run(command, cwd=self.directory, capture_output=True, text=True)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)

    def test_accepts_either_surviving_diagnostic(self):
        for line in (1, 2):
            for logfile in ("", f"interface-assertion-R{line-1}.log"):
                with self.subTest(line=line, logfile=logfile):
                    self.run_checker(line, 0, logfile=logfile)

    def test_single_assertion_remains_strict(self):
        self.run_checker(1, 0, markers=("write guard",))
        self.run_checker(2, 1, markers=("write guard",))

    def test_rejects_unrelated_failure(self):
        self.run_checker(3, 1)

    def test_rejects_success_even_with_valid_diagnostic(self):
        self.run_checker(2, 1, status=0)

    def test_rejects_failure_without_diagnostic(self):
        self.run_checker(None, 1)

    def test_stale_rank_log_cannot_satisfy_expectation(self):
        (self.directory / "interface-assertion-R1.log").write_text(
            f"Panic! Assertion failed at {self.source}:2\n")
        self.run_checker(None, 1)


if __name__ == "__main__":
    unittest.main()
