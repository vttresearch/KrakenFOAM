# KrakenFOAM

KrakenFOAM couples [OpenFOAM](https://openfoam.org) solvers to VTT's
[Kraken](https://serpent.vtt.fi/kraken/index.php?title=Main_Page)
computational framework for reactor analysis. The multiphysics driver
**Cerberus** controls the OpenFOAM run. Cerberus and OpenFOAM communicate
through a socket: Cerberus sends the time intervals to solve and the input
fields, such as power or temperature distributions. OpenFOAM sends back output
fields such as temperatures, densities and heat transfer rates.

KrakenFOAM contains:

- `libkrakenIO`: a library with the socket communication layer and a set of
  `fvModels` and `fvConstraints` that exchange field data with Cerberus.
- `krakenRun` and `krakenMultiRun`: modified versions of `foamRun` and
  `foamMultiRun`. Cerberus controls when they run.
- Tutorial cases showing how to use the models.

## Requirements

- [OpenFOAM 13](https://openfoam.org/version/13/) (openfoam.org)
- [Cerberus](https://serpent.vtt.fi/kraken/index.php?title=Main_Page) (Python
  package, part of Kraken) to run coupled simulations.
- Linux or another POSIX system with a C/C++ compiler supported by OpenFOAM.

## Installation

Source the OpenFOAM 13 environment, then build:

```sh
source /opt/openfoam13/etc/bashrc
git clone https://github.com/vttresearch/KrakenFOAM.git
cd KrakenFOAM
./Allwmake -j
```

The library `libkrakenIO.so` is installed in `$FOAM_USER_LIBBIN`. The
applications `krakenRun` and `krakenMultiRun` are installed in
`$FOAM_USER_APPBIN`. To remove the build files, run `./Allwclean`.

## Contents

```
KrakenFOAM
├── applications
│   ├── krakenRun          Single region solver driven by Cerberus
│   └── krakenMultiRun     Multi-region solver driven by Cerberus
├── src
│   └── krakenIO           Communication library and coupling models
└── tutorials
    ├── MiniCore-P         krakenMappedScalarFieldSource with krakenRun
    ├── MiniCore-T         krakenFixedTemperatureConstraint with krakenRun
    └── SingleRodCHT       krakenScalarFieldSource with krakenMultiRun
```

### Library `krakenIO`

| Class | Type | Description |
|---|---|---|
| `krakenIO` | `regIOobject` | Socket communication with Cerberus (`clientIO.c` holds the low-level C socket routines) |
| `krakenModelBase` | base class | Common interface for the `fvModels` and `fvConstraints` that exchange data with Cerberus |
| `krakenScalarFieldSource` | `fvModel` | Power source given for **every cell** of the mesh. Returns cell temperatures and densities. Single phase only. |
| `krakenMappedScalarFieldSource` | `fvModel` | Power source given for a **set of points**. Each cell in the `cellZone` is mapped to the nearest point. Returns point-wise temperatures and densities volume averaged over `outputPhases`. Multiphase only. |
| `krakenFixedTemperatureConstraint` | `fvConstraint` | Fixes the temperature of `inputPhase` to values given for a **set of points**. Each cell in the `cellZone` is mapped to the nearest point. Returns point-wise heat transfer rates (powers) from the `inputPhase` together with averaged temperatures and densities volume averaged over `outputPhases`. Multiphase only. |

### Applications

`krakenRun` and `krakenMultiRun` take the same options and `controlDict`
entries (`solver`, `regionSolvers`) as `foamRun` and `foamMultiRun`. The
difference is that the solver does not march through time on its own. It
waits for Cerberus to send instructions such as "advance to time interval",
"solve", "send/receive field" or "exit". To solve a time interval, the solver
sets the end time of the run to the end of the interval Cerberus requested and
runs the normal time loop up to that point.

## Usage

### Connection settings: `system/krakenIO`

Each coupled case needs a `system/krakenIO` dictionary:

```cpp
hostName   localhost;   // Host running Cerberus
port       2230;        // First port to try
unixSocket true;        // Use a unix socket (./.socket_<port>) instead of TCP
```

The client tries ports from `port` down to 2211 until a connection succeeds.

### Field naming

The names of the exchanged fields include the name of the model or constraint
dictionary entry:

| Pattern | Direction |
|---|---|
| `foam_if_<quantity>_<name>` | Input field, Cerberus → OpenFOAM |
| `foam_of_<quantity>_<name>` | Output field, OpenFOAM → Cerberus |
| `foam_mesh_<name>` | Mesh associated with the fields |

| Model | Inputs | Outputs |
|---|---|---|
| `krakenScalarFieldSource` | `power` [W] | `temperature` [K], `density` [kg/m³] |
| `krakenMappedScalarFieldSource` | `power` [W] | `temperature` [K], `density` [kg/m³] |
| `krakenFixedTemperatureConstraint` | `temperature` [K] | `temperature` [K], `density` [kg/m³], `power` [W] |

### Examples

Per-cell power source (`constant/fvModels`). In parallel runs the
undecomposed cell numbering is recovered from `cellProcAddressing`:

```cpp
fuel
{
    type            krakenScalarFieldSource;
    libs            ("libkrakenIO.so");
}
```

Point-mapped power source (`constant/fvModels`):

```cpp
core
{
    type            krakenMappedScalarFieldSource;
    libs            ("libkrakenIO.so");

    cellZone        all;

    inputPhase      solid;
    outputPhases    (gas liquid);

    points
    (
        #include "${FOAM_CASE}/all_cell_centres"
    );
}
```

Fixed temperature constraint (`system/fvConstraints`):

```cpp
core
{
    type            krakenFixedTemperatureConstraint;
    libs            ("libkrakenIO.so");

    cellZone        all;

    inputPhase      solid;
    outputPhases    (gas liquid);

    points
    (
        #include "${FOAM_CASE}/all_cell_centres"
    );
}
```

## Tutorials

| Case | Solver | Coupling model |
|---|---|---|
| `MiniCore-P` | `krakenRun` (`multiphaseEuler`) | `krakenMappedScalarFieldSource` |
| `MiniCore-T` | `krakenRun` (`multiphaseEuler`) | `krakenFixedTemperatureConstraint` |
| `SingleRodCHT` | `krakenMultiRun` (fluid + solid) | `krakenScalarFieldSource` |

Each case has an `input.py` script that uses the Cerberus Python package to
start OpenFOAM, send the input fields, advance time and collect the outputs.
The scripts launch OpenFOAM through the path in the `FOAMJOB_EXE_PATH`
environment variable, normally the OpenFOAM `foamJob` script:

```sh
source /opt/openfoam13/etc/bashrc
export FOAMJOB_EXE_PATH=$WM_PROJECT_DIR/bin/foamJob
cd tutorials/MiniCore-P
python3 input.py
```

Notes:

- The MiniCore cases include a precomputed point list `all_cell_centres`.
- `SingleRodCHT` builds its mesh with `./Allinit`, which `input.py` calls.
- Parallel runs are on by default in the scripts (`parallel_run`, `n_proc`).
- Run `./Allclean` to reset a case.

## Coding style

The code follows the
[OpenFOAM coding style](https://openfoam.org/dev/coding-style-guide/) used in
OpenFOAM 13. See [CONTRIBUTING.md](CONTRIBUTING.md).

## Citing

If you use KrakenFOAM in your work, please cite it as described in
[CITATION.cff](CITATION.cff).

## License

KrakenFOAM is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. See [LICENSE](LICENSE).

Parts of KrakenFOAM are derived from OpenFOAM, Copyright (C) OpenFOAM
Foundation. OpenFOAM® is a registered trademark of OpenCFD Limited, producer
and distributor of the OpenFOAM software via www.openfoam.com. KrakenFOAM is not
approved or endorsed by OpenCFD Limited or the OpenFOAM Foundation.

## Contact

Riku Tuominen, VTT Technical Research Centre of Finland Ltd —
[riku.tuominen@vtt.fi](mailto:riku.tuominen@vtt.fi)
