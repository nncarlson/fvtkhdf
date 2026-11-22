"""Check local fatal termination when a test-only MPI_Abort returns."""

import argparse
from pathlib import Path
import re
import subprocess
import sys


def check(executable, scenario):
    directory = Path(scenario)
    directory.mkdir(exist_ok=True)
    sentinel = directory / "after-assertion"
    sentinel.unlink(missing_ok=True)
    try:
        result = subprocess.run(
            [str(executable), scenario], cwd=directory,
            capture_output=True, text=True, errors="replace", timeout=10)
    except (OSError, subprocess.TimeoutExpired) as error:
        sys.stderr.write(f"{scenario}: could not complete assertion test: {error}\n")
        return False

    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)
    output = result.stdout + result.stderr
    lines = output.splitlines()
    panics = [line for line in lines if re.match(
        r"^(?:\[R\d+\])?Panic! Assertion failed at ", line)]
    abort = "ASSERT_MPI_ABORT_RETURNED comm=123 code=1"
    continued = "ASSERT_AFTER" in output or sentinel.exists()
    entered = lines.count("ASSERT_BEFORE") == 1
    if scenario == "pass":
        valid = (entered and result.returncode == 0 and not panics
                 and "ASSERT_MPI_ABORT_RETURNED" not in output
                 and lines.count("ASSERT_AFTER") == 1 and sentinel.is_file()
                 and sentinel.read_text().strip() == "ASSERT_AFTER")
    else:
        valid = (entered and result.returncode > 0 and not continued
                 and len(panics) == 1 and re.fullmatch(
                     r"\[R0\]Panic! Assertion failed at .*vtkhdf_assert_test\.F90:\d+",
                     panics[0])
                 and lines.count(abort) == 1 and "error stop" in output.lower())
    print(f"{scenario}: exit={result.returncode}, post_call={continued}")
    if not valid:
        sys.stderr.write(f"{scenario}: assertion termination contract violated\n")
    return bool(valid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve()
    passed = check(executable, "pass")
    failed = check(executable, "fail")
    return 0 if passed and failed else 1


if __name__ == "__main__":
    sys.exit(main())
