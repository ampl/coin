### COIN-OR projects with CMake build support

- [Cbc](https://projects.coin-or.org/Cbc) version 2.10.13, built as a dependency of Bonmin and Couenne ([source](Cbc/))
- [Ipopt](https://projects.coin-or.org/Ipopt) version 3.14.20 ([source](Ipopt/))
- [Bonmin](https://projects.coin-or.org/Bonmin) version 1.8.9 ([source](Bonmin/))
- [Couenne](https://projects.coin-or.org/Couenne) version 0.5.8 ([source](Couenne/))
- [ASL](http://netlib.org/ampl/) version 20260109 ([source](https://github.com/ampl/asl))

GitHub repository: https://github.com/ampl/coin/

All binaries are available for download from the [AMPL Open-Source Downloads](http://ampl.com/dl/open/) page.

By default, the build pipeline builds and packages Bonmin, Couenne, and Ipopt separately as
`solvers-<solver>.<platform>.<YYYYMMDD>.tar.gz`, for example
`solvers-ipopt.linux-intel64.20261027.tar.gz`. Each archive contains one
solver executable, `docs/coin-license.txt`, `docs/coin-versions.txt`, and
`docs/README.<solver>.txt`. Bonmin, Couenne, and Ipopt packages also contain
`docs/CHANGES.<solver>.md`. Changelogs use dated `YYYYMMDD` sections with
bullet lists of updates, newest first. Windows archives
for Ipopt, Bonmin, and Couenne also include `libipoptfort.dll`.

The macOS job also produces `signed-solvers-<solver>.<platform>.<YYYYMMDD>.tar.gz`
archives containing code-signed binaries. It uses the same Developer ID Application
certificate and entitlements as `misc`, with hardened runtime and a secure timestamp,
and verifies every signature before archiving. Signing operates on staged copies.
The Azure pipeline needs access to the Secure File
`Developer ID Application- AMPL Optimization Inc. (ZNNBG5892S).p12` and the secret
`certPwd` variable, as in `misc`.

To package an existing build locally, run:

```sh
bash package-solvers.sh build dist linux-intel64
```

On macOS, with the signing certificate installed in your keychain:

```sh
bash package-solvers.sh build dist macos-arm64 --signed
```

The default identity is `Developer ID Application: AMPL Optimization Inc. (ZNNBG5892S)`.
Set `MACOS_SIGNING_IDENTITY` to override it and `MACOS_SIGNING_KEYCHAIN` to select
a particular keychain. These archives contain signed binaries; notarization is
a separate step.

The platform labels used by the pipeline are `linux-intel64`, `linux-arm64`,
`macos-arm64`, and `mswin64`. Each archive uses the newest `YYYYMMDD` heading in its solver changelog
as its package version, independent of the build date.

The standalone Cbc solver is disabled by default. To enable it, configure with
`-DBUILD_CBC_EXECUTABLE=ON`; set it back to `OFF` to disable it again. The Cbc
libraries required by Bonmin and Couenne are always built. Both packaging entry
points follow the configured solver list in `<build>/coin-solvers.txt`, so old
Cbc binaries left in the build directory are excluded when the option is off.
When enabled, Cbc packages use the newest release date among the solver changelogs.

The CMake target `ampl-module` also creates these individual archives in
`<build>/ampl-module`, building the solver executables and version metadata
first. Set `-DAMPL_PLATFORM=...` when configuring to override the platform label.

Release builds use the pinned ASL submodule revision. CMake checks solver,
MUMPS, and ASL source versions against the latest changelog entries and
fails configuration if they differ.
