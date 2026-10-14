#!/usr/bin/env python3
"""Plot the four comparisons against Marchi et al. (2009), Table 6.

Requires numpy and matplotlib. Writes PNG files to results/plots/.
"""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
RES = ROOT/'results'
DEST = RES/'plots'

with (ROOT/'data/marchi_2009_table6_centerlines.csv').open() as f:
    REFS = list(csv.DictReader(f))

def reference(re: float, component: str):
    r = [x for x in REFS if float(x['Re']) == re and x['component'] == component]
    return [float(x['coordinate']) for x in r], [float(x['velocity']) for x in r]

def profile(case: str, component: str):
    a = np.genfromtxt(RES/case/f'centerline_{component}.csv', delimiter=',', names=True)
    return a['y' if component == 'u' else 'x'], a[component]

def centerlines(cases: dict[str, str], re: float, title: str, name: str) -> None:
    """Left: u(y) at x=0.5; right: v(x) at y=0.5; one line per case."""
    fig, (au, av) = plt.subplots(1, 2, figsize=(11, 4.8))
    for label, case in cases.items():
        if not (RES/case).exists():
            continue
        y, u = profile(case, 'u'); au.plot(u, y, label=label)
        x, v = profile(case, 'v'); av.plot(x, v, label=label)
    c, val = reference(re, 'u'); au.plot(val, c, 'ko', ms=4, label='Marchi et al. (2009)')
    c, val = reference(re, 'v'); av.plot(c, val, 'ko', ms=4, label='Marchi et al. (2009)')
    au.set(xlabel='u / U_lid', ylabel='y / L', title='u at x = 0.5')
    av.set(xlabel='x / L', ylabel='v / U_lid', title='v at y = 0.5')
    for a in (au, av):
        a.grid(True, alpha=0.3); a.legend(fontsize=8)
    fig.suptitle(title); fig.tight_layout()
    fig.savefig(DEST/name, dpi=160); plt.close(fig)

def reynolds_study() -> None:
    fig, (au, av) = plt.subplots(1, 2, figsize=(11, 4.8))
    for i, re in enumerate((100, 400, 1000)):
        case = f'simple_64_Re{re}_cds'
        if not (RES/case).exists():
            continue
        col = f'C{i}'
        y, u = profile(case, 'u'); au.plot(u, y, col, label=f'Re = {re}')
        x, v = profile(case, 'v'); av.plot(x, v, col, label=f'Re = {re}')
        c, val = reference(re, 'u'); au.plot(val, c, 'o', color=col, ms=4)
        c, val = reference(re, 'v'); av.plot(c, val, 'o', color=col, ms=4)
    au.set(xlabel='u / U_lid', ylabel='y / L', title='u at x = 0.5')
    av.set(xlabel='x / L', ylabel='v / U_lid', title='v at y = 0.5')
    for a in (au, av):
        a.grid(True, alpha=0.3); a.legend(fontsize=8)
    fig.suptitle('Reynolds-number study: SIMPLE, 64 x 64, CDS (lines) vs Marchi et al. (symbols)')
    fig.tight_layout(); fig.savefig(DEST/'reynolds_study.png', dpi=160); plt.close(fig)

def coupling_convergence() -> None:
    """Residual history of each coupling algorithm."""
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for label, case in {'SIMPLE (iteration)': 'simple_32_Re100', 'PISO (time step)': 'piso_32_Re100',
                        'PIMPLE (time step)': 'pimple_32_Re100'}.items():
        path = RES/case/'residuals.csv'
        if not path.exists():
            continue
        a = np.genfromtxt(path, delimiter=',', names=True)
        ax.semilogy(a['step'], np.maximum(a['u_momentum_inf'], a['v_momentum_inf']), label=label)
    ax.set(xlabel='Outer iteration / time step', ylabel='Steady momentum residual (inf-norm)',
           title='Convergence to steady state: Re = 100, 32 x 32')
    ax.legend(); ax.grid(True, alpha=0.3); fig.tight_layout()
    fig.savefig(DEST/'coupling_convergence.png', dpi=160); plt.close(fig)

def streamlines() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))
    for ax, re in zip(axes, (100, 400, 1000)):
        path = RES/f'simple_64_Re{re}_cds/fields.csv'
        if not path.exists():
            continue
        f = np.genfromtxt(path, delimiter=',', names=True)
        x, y = np.unique(f['x']), np.unique(f['y'])
        u, v = f['u'].reshape(len(y), len(x)), f['v'].reshape(len(y), len(x))
        ax.streamplot(x, y, u, v, density=1.6, color=np.hypot(u, v), cmap='viridis', linewidth=0.8)
        ax.set(xlim=(0, 1), ylim=(0, 1), title=f'Re = {re}', xlabel='x / L', ylabel='y / L')
        ax.set_aspect('equal')
    fig.suptitle('Streamlines, 64 x 64'); fig.tight_layout()
    fig.savefig(DEST/'streamlines.png', dpi=160); plt.close(fig)

if __name__ == '__main__':
    DEST.mkdir(exist_ok=True)
    centerlines({'SIMPLE': 'simple_32_Re100', 'PISO': 'piso_32_Re100', 'PIMPLE': 'pimple_32_Re100'},
                100, 'Pressure-velocity coupling: Re = 100, 32 x 32', 'coupling_centerlines.png')
    centerlines({f'{n} x {n}': f'simple_{n}_Re100' for n in (16, 32, 64)},
                100, 'Mesh study: SIMPLE, Re = 100', 'mesh_study.png')
    centerlines({'Hybrid (CDS/upwind)': 'simple_32_Re100', 'Upwind': 'simple_32_Re100_upwind',
                 'CDS': 'simple_32_Re100_cds'},
                100, 'Convection scheme: SIMPLE, Re = 100, 32 x 32', 'scheme_comparison.png')
    reynolds_study()
    coupling_convergence()
    streamlines()
