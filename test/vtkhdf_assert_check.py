"""Check local fatal termination when a test-only MPI_Abort returns."""

import argparse
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve()
    status = 0
    for scenario in ("pass", "fail"):
        result = subprocess.run(
            [str(executable), scenario], capture_output=True, text=True, errors="replace")
        sys.stdout.write(result.stdout)
        sys.stderr.write(result.stderr)
        aborted = "ASSERT_MPI_ABORT_RETURNED" in result.stdout + result.stderr
        if scenario == "pass":
            valid = result.returncode == 0 and not aborted
        else:
            valid = result.returncode != 0 and aborted
        if not valid:
            sys.stderr.write(
                f"{scenario}: unexpected exit={result.returncode}, MPI_Abort called={aborted}\n")
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
