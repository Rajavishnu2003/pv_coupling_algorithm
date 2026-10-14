#!/usr/bin/env python3
"""Compare converged cavity profiles with Marchi et al. (2009), Table 6.

Uses only Python's standard library. Linear interpolation at the 15 published
interior coordinates; velocities and error norms are normalised by U_lid=1.
Failed/non-steady runs are recorded but excluded from benchmark error norms.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from bisect import bisect_right
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def interpolate(x: list[float], y: list[float], point: float) -> float:
    if len(x) != len(y) or len(x) < 2 or any(b <= a for a, b in zip(x, x[1:])):
        raise ValueError("Invalid interpolation grid")
    if point < x[0] or point > x[-1]:
        raise ValueError("Reference point outside numerical profile")
    k = max(0, min(len(x)-2, bisect_right(x, point)-1))
    weight = (point-x[k])/(x[k+1]-x[k])
    return (1-weight)*y[k] + weight*y[k+1]

def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))

def analyze(results: Path) -> tuple[list[dict], list[dict]]:
    references = rows(ROOT / 'data/marchi_2009_table6_centerlines.csv')
    summaries, differences = [], []
    for file in sorted(results.glob('*/summary.json')):
        s = json.loads(file.read_text())
        row = {'case': file.parent.name, **s}
        row['momentum_inf'] = max(s['u_momentum_inf'], s['v_momentum_inf'])
        total = s['cds_faces'] + s['upwind_faces']
        row['upwind_percent'] = 100*s['upwind_faces']/total
        for component in ('u', 'v'):
            row[f'{component}_rmse'] = ''
            row[f'{component}_max_error'] = ''
            row[f'{component}_center'] = ''
            profile = rows(file.parent / f'centerline_{component}.csv')
            coord = 'y' if component == 'u' else 'x'
            x, y = [float(p[coord]) for p in profile], [float(p[component]) for p in profile]
            row[f'{component}_center'] = interpolate(x, y, 0.5)
            if not (s['completed'] and s['converged']):
                continue
            errors = []
            for ref in references:
                if float(ref['Re']) == s['Re'] and ref['component'] == component:
                    numerical = interpolate(x, y, float(ref['coordinate']))
                    error = numerical-float(ref['velocity'])
                    errors.append(error)
                    differences.append({'case':file.parent.name, 'Re':s['Re'],
                        'component':component, 'coordinate':ref['coordinate'],
                        'computed':numerical, 'reference':ref['velocity'],
                        'signed_error':error, 'absolute_error':abs(error)})
            if errors:
                if len(errors) != 15:
                    raise ValueError('Expected exactly 15 reference values per component')
                row[f'{component}_rmse'] = math.sqrt(sum(e*e for e in errors)/len(errors))
                row[f'{component}_max_error'] = max(abs(e) for e in errors)
        summaries.append(row)
    return summaries, differences

def write_table(path: Path, data: list[dict]) -> None:
    if data:
        with path.open('w', newline='') as out:
            writer = csv.DictWriter(out, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=ROOT/'results')
    args = parser.parse_args()
    summary, detail = analyze(args.results)
    write_table(args.results/'validation_summary.csv', summary)
    write_table(args.results/'benchmark_point_errors.csv', detail)
    for s in summary:
        status = ('CONVERGED' if s['converged'] else 'TRANSIENT END') if s['completed'] else 'FAILED'
        print(f"{s['case']:30s} {status:11s} steps={s['steps']:5d} "
              f"u_RMSE={s['u_rmse']} v_RMSE={s['v_rmse']}")

if __name__ == '__main__':
    main()
