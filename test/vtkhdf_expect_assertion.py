"""Require a failing command and the diagnostic from a specific assertion."""

import argparse
from pathlib import Path
import re
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--assertion", required=True,
                        help="unique text identifying the assertion's source line")
    parser.add_argument("--diagnostic-glob",
                        help="read per-rank diagnostic files if MPI abort loses captured output")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        parser.error("a command is required after --")

    # Fypp preserves template line numbers. Resolve the location from source so
    # adding lines above the assertion does not invalidate the expected result.
    lines = [number for number, line in enumerate(args.source.read_text().splitlines(), 1)
             if args.assertion in line]
    if len(lines) != 1:
        parser.error("--assertion must identify exactly one source line")
    location = f"{args.source.name}:{lines[0]}"
    diagnostic = rf"Panic! Assertion failed at [^\r\n]*{re.escape(location)}(?!\d)"

    if args.diagnostic_glob:
        # Never accept a diagnostic left by an earlier invocation.
        for path in Path.cwd().glob(args.diagnostic_glob):
            path.unlink()
    result = subprocess.run(command, capture_output=True, text=True, errors="replace")
    diagnostics = ""
    if args.diagnostic_glob:
        diagnostics = "".join(path.read_text(errors="replace")
                              for path in sorted(Path.cwd().glob(args.diagnostic_glob)))
    sys.stdout.write(result.stdout)
    sys.stderr.write(result.stderr)
    if result.returncode == 0:
        sys.stderr.write("Expected assertion failure, but command exited successfully\n")
        return 1
    if not re.search(diagnostic, result.stdout + result.stderr + diagnostics):
        sys.stderr.write(diagnostics)
        sys.stderr.write(f"Expected assertion at {location}; command exited {result.returncode}\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
