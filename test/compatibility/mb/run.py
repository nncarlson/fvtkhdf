"""Generate the small MB matrix, or read existing fixtures with a chosen ParaView."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
MODES = ("static", "fixed", "moving", "dynamic")
CONFIGURATIONS = [(kind, mode) for kind, limit in
                  (("image", 2), ("rectilinear", 3), ("structured", 3), ("unstructured", 4))
                  for mode in range(limit)]


def matrix():
    for kind, mode in CONFIGURATIONS:
        for form in ("standalone", "mb"):
            yield f"{form}-{kind}-{MODES[mode]}", form == "mb", [(kind, mode)]
    kinds = ("image", "rectilinear", "structured", "unstructured")
    yield "all-static", True, [(kind, 0) for kind in kinds]
    yield "all-fixed", True, [(kind, 1) for kind in kinds]
    moving = [("image", 1), ("rectilinear", 2), ("structured", 2), ("unstructured", 3)]
    mixed = [("image", 0), ("rectilinear", 1), ("structured", 2), ("unstructured", 3)]
    for name, blocks in (("different-modes", moving), ("mixed-modes", mixed),
                         ("unstructured-mixed", [("unstructured", 0), ("unstructured", 1)]),
                         ("rectilinear-mixed", [("rectilinear", 0), ("structured", 1)]),
                         ("structured-mixed", [("structured", 0), ("image", 1)])):
        yield name, True, blocks
        yield name + "-reversed", True, list(reversed(blocks))
    yield "unstructured-workaround", True, [("unstructured", "fixed-constant"), ("unstructured", 1)]
    yield "rectilinear-workaround", True, [("rectilinear", "fixed-constant"), ("structured", 1)]


def execute(command, log, timeout):
    try:
        with log.open("w") as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
        return result.returncode
    except subprocess.TimeoutExpired:
        with log.open("a") as stream:
            stream.write(f"\nWorker timed out after {timeout} seconds\n")
        return "timeout"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generator", type=Path, help="executable built by this project's CMakeLists.txt")
    parser.add_argument("--pvpython", help="pvpython executable for the version under test")
    parser.add_argument("--work-dir", type=Path, required=True, help="generated fixtures and manifest")
    parser.add_argument("--results-dir", type=Path, help="separate directory for this reader version's results")
    parser.add_argument("--writer-label", default="unspecified", help="installed writer version or commit")
    parser.add_argument("--read-only", action="store_true", help="reuse fixtures without rebuilding or generating")
    parser.add_argument("--generate-only", action="store_true")
    parser.add_argument("--timeout", type=int, default=90, help="timeout in seconds for each worker")
    args = parser.parse_args()
    if args.read_only and args.generate_only:
        parser.error("--read-only and --generate-only are mutually exclusive")
    if not args.read_only and not args.generator:
        parser.error("--generator is required for generation")
    if not args.generate_only and not args.pvpython:
        parser.error("--pvpython is required for readback")
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    manifest = work / "manifest.json"
    if not args.read_only:
        cases = []
        for name, composite, blocks in matrix():
            case = {"name": name, "composite": composite, "blocks": blocks,
                    "file": str(work / (name + ".vtkhdf"))}
            case_path = work / (name + ".json")
            case_path.write_text(json.dumps(case, indent=2) + "\n")
            command = [str(args.generator.resolve()), case["file"], "mb" if composite else "standalone"]
            command += [str(value) for block in blocks for value in block]
            code = execute(command, work / (name + ".generate.log"), args.timeout)
            if code != 0:
                raise RuntimeError(f"Generator failed for {name}: {code}; see {work}")
            cases.append(case)
        manifest.write_text(json.dumps({"writer": args.writer_label, "cases": cases}, indent=2) + "\n")
    fixtures = json.loads(manifest.read_text())
    # Validate the payload independently before attributing failures to a reader.
    for case in fixtures["cases"]:
        name = case["name"]
        result = work / (name + ".layout.json")
        result.unlink(missing_ok=True)
        code = execute([sys.executable, str(HERE / "check.py"), "layout", str(work / (name + ".json")), str(result)],
                       work / (name + ".layout.log"), args.timeout)
        if code != 0 or not result.exists() or json.loads(result.read_text())["status"] != "pass":
            raise RuntimeError(f"Layout check failed for {name}; see {result}")
    print(f"Validated {len(fixtures['cases'])} fixtures", flush=True)
    if args.generate_only:
        return 0
    executable = shutil.which(args.pvpython)
    if executable is None:
        parser.error(f"Cannot find pvpython: {args.pvpython}")
    results_dir = (args.results_dir or work / "results").resolve()
    results_dir.mkdir(parents=True, exist_ok=True)
    results = {"writer": fixtures["writer"], "pvpython": executable, "cases": {}}
    for case in fixtures["cases"]:
        name = case["name"]
        path = results_dir / (name + ".json")
        path.unlink(missing_ok=True)
        code = execute([executable, "--disable-registry", str(HERE / "check.py"), "reader",
                        str(work / (name + ".json")), str(path)], results_dir / (name + ".log"), args.timeout)
        if code != 0 or not path.exists():
            results["cases"][name] = {"status": "worker-failed", "exit_code": code}
        else:
            value = json.loads(path.read_text())
            if "profiles" not in value:
                value.update(status="worker-failed", exit_code=code)
            results["cases"][name] = value
        print(name + ": " + ", ".join(f"{profile}={value['status']}" for profile, value in
                                        results["cases"][name].get("profiles", {}).items()), flush=True)
    lines = [f"Writer: {results['writer']}", f"pvpython: {executable}", "",
             "| Case | Geometry | Point | Cell | Field | All | Fresh later step | Fetch |",
             "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    worker_failure = False
    for name, value in results["cases"].items():
        if value.get("status") == "worker-failed":
            worker_failure = True
        else:
            results["paraview"] = value["paraview"]
            results["vtk"] = value["vtk"]
        cells = [value.get("profiles", {}).get(profile, {}).get("status", "worker failed")
                 for profile in ("geometry", "point", "cell", "field", "all", "late", "fetch")]
        lines.append("| " + " | ".join([name] + cells) + " |")
    lines.insert(2, f"ParaView: {results.get('paraview', 'unknown')}; VTK: {results.get('vtk', 'unknown')}")
    (results_dir / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    (results_dir / "summary.md").write_text("\n".join(lines) + "\n")
    print(f"Results: {results_dir / 'summary.md'}")
    # Recorded incompatibilities are expected outcomes of this investigation;
    # generation, layout validation, and worker crashes/timeouts are harness failures.
    return int(worker_failure)


if __name__ == "__main__":
    raise SystemExit(main())
