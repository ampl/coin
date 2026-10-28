### COIN-OR projects with CMake build support

- [Cbc](https://projects.coin-or.org/Cbc) version 2.10.13 ([source](Cbc/))
- [Ipopt](https://projects.coin-or.org/Ipopt) version 3.14.20 ([source](Ipopt/))
- [Bonmin](https://projects.coin-or.org/Bonmin) version 1.8.9 ([source](Bonmin/))
- [Couenne](https://projects.coin-or.org/Couenne) version 0.5.8 ([source](Couenne/))
- [ASL](http://netlib.org/ampl/) version 20260109 ([source](https://github.com/ampl/asl))

GitHub repository: https://github.com/ampl/coin/

All binaries are available for download from the [AMPL Open-Source Downloads](http://ampl.com/dl/open/) page.

The build pipeline packages each solver separately as
`solvers-<solver>.<platform>.<YYYYMMDD>.tar.gz`, for example
`solvers-ipopt.linux-intel64.20261027.tar.gz`. Each archive contains one
solver executable, `docs/coin-license.txt`, `docs/coin-versions.txt`, and
`docs/README.<solver>.txt`. Bonmin, Couenne, and Ipopt packages also contain
`docs/CHANGES.<solver>.md`. Changelogs use dated `YYYYMMDD` sections with
bullet lists of updates, newest first. Windows archives
for Ipopt, Bonmin, and Couenne also include `libipoptfort.dll`.

To package an existing build locally, run:

```sh
bash package-solvers.sh build dist linux-intel64
```

The platform labels used by the pipeline are `linux-intel64`, `linux-arm64`,
`macos-arm64`, and `mswin64`. Each archive uses the newest `YYYYMMDD` heading in its solver changelog
as its package version, independent of the build date. Cbc has no changelog
and uses the newest release date among the solver changelogs.

The CMake target `ampl-module` also creates these individual archives in
`<build>/ampl-module`, building the solver executables and version metadata
first. Set `-DAMPL_PLATFORM=...` when configuring to override the platform label.

Release builds use the pinned ASL submodule revision. CMake checks solver,
MUMPS, and ASL source versions against the latest changelog entries and
fails configuration if they differ.
