# Lid-Driven Cavity: SIMPLE, PISO and PIMPLE

## Requirements

- A C++17 compiler (for example g++) and CMake 3.16 or later
- Python 3 with `numpy` and `matplotlib` (needed only for the plots)

## Compile

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
```

This creates the executable `build/cavity`.

## Run

Run all cases and compare them with the literature data (Marchi et al., 2009, in `data/`):

```sh
python3 scripts/run_suite.py
python3 scripts/plot_results.py
```

Or run a single case:

```sh
./build/cavity cases/simple_32_Re100.cfg
```

## Results

- `results/<case>/`: centreline velocities, residual history, full fields and a `cavity.vtk` file for ParaView
- `results/validation_summary.csv`: error of each case against the literature data
- `results/plots/`: comparison plots

## Cases (`cases/`)

| Comparison | Case files |
|---|---|
| Coupling algorithm (Re = 100, 32×32) | `simple_32_Re100`, `piso_32_Re100`, `pimple_32_Re100` |
| Mesh size (SIMPLE, Re = 100) | `simple_16_Re100`, `simple_32_Re100`, `simple_64_Re100` |
| Reynolds number (SIMPLE, 64×64) | `simple_64_Re100_cds`, `simple_64_Re400_cds`, `simple_64_Re1000_cds` |
| Convection scheme (Re = 100, 32×32) | `simple_32_Re100` (hybrid), `simple_32_Re100_upwind`, `simple_32_Re100_cds` |

To run a new case, edit a `.cfg` file. The main settings are `nx`, `ny`,
`reynolds`, `algorithm` (simple / piso / pimple), `convection` (peclet / upwind / cds),
`dt` and `output`.
