#!/usr/bin/env python3
"""Run all assignment cases, then compare them with Marchi et al. (2009).

Build the C++ program first. Cases run in parallel; each writes results/<case>/.
Comparisons:
  1. Coupling:  SIMPLE vs PISO vs PIMPLE      (Re=100, 32x32)
  2. Mesh:      SIMPLE on 16, 32, 64 cells    (Re=100)
  3. Reynolds:  SIMPLE, Re = 100, 400, 1000   (64x64, CDS)
  4. Scheme:    hybrid vs upwind vs CDS       (Re=100, 32x32)
"""
from __future__ import annotations
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import validate
ROOT = Path(__file__).resolve().parents[1]
CASES = sorted(p.stem for p in (ROOT/'cases').glob('*.cfg'))

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe', type=Path, default=ROOT/'build/cavity')
    parser.add_argument('-j', '--jobs', type=int, default=4)
    args = parser.parse_args()
    exe = args.exe.resolve()
    if not exe.is_file():
        parser.error(f'Executable does not exist: {exe}; build first or use --exe')
    def run(case: str) -> int:
        code = subprocess.run([str(exe), str(ROOT/'cases'/f'{case}.cfg')], cwd=ROOT,
                              stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT).returncode
        print(f'{case}: exit {code}', flush=True)
        return code
    with ThreadPoolExecutor(args.jobs) as pool:
        codes = list(pool.map(run, CASES))
    summary, detail = validate.analyze(ROOT/'results')
    validate.write_table(ROOT/'results/validation_summary.csv', summary)
    validate.write_table(ROOT/'results/benchmark_point_errors.csv', detail)
    return 0 if all(c == 0 for c in codes) else 1

if __name__ == '__main__':
    raise SystemExit(main())
