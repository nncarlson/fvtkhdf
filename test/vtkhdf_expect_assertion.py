"""Require a failing command and a diagnostic from an allowed assertion."""

import argparse
from pathlib import Path
import re
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--assertion", required=True, action="append",
                        help="unique text identifying an allowed assertion's source line; repeat for alternatives")
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
    source_lines = args.source.read_text().splitlines()
    locations = []
    patterns = []
    for assertion in args.assertion:
        lines = [number for number, line in enumerate(source_lines, 1) if assertion in line]
        if len(lines) != 1:
            parser.error(f"--assertion {assertion!r} must identify exactly one source line")
        location = f"{args.source.name}:{lines[0]}"
        locations.append(location)
        patterns.append(rf"Panic! Assertion failed at [^\r\n]*{re.escape(location)}(?!\d)")

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
    output = result.stdout + result.stderr + diagnostics
    if not any(re.search(pattern, output) for pattern in patterns):
        sys.stderr.write(diagnostics)
        sys.stderr.write(f"Expected assertion at {' or '.join(locations)}; command exited {result.returncode}\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
