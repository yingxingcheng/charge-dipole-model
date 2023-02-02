This fortran demo program implements the charge-dipole interaction model with Gaussian-type density functions.
See [this link](https://link.aps.org/doi/10.1103/PhysRevB.75.045407) for more details.

The author acknowledges Prof. Alexandre Mayer for fixing the issues of the codes and providing C60 and C70 coordinates.

# Installation

One needs to install this program first by running

```shell
make
```
and it is assumed that `gfortron` compiler and library `llapack` has been installed.

# Usage

Run `./main.x {xyz_file}` to compute the molecular dipole polarizability and `xyz_file` is molecular geometry file with `xyz` format, e.g. `C60.xyz`:

```shell
./main.x C60.xyz
```

One should obtain the output as follows:

```shell
   |---------------------------Input Geometry-----------------------------------------
   All coordinates are in unit of Angstrom.
   Number of atoms   60

    C          0.001520        -1.213990        -3.310090
    C         -1.157700        -0.371730        -3.310140
    C          1.160770        -0.371770        -3.310050
    C         -0.714880         0.991030        -3.310140
    C          0.718010         0.991000        -3.310090
    C          0.001470        -2.422190        -2.563360
    C         -2.306820        -0.745050        -2.563470
    C          2.309820        -0.745140        -2.563290
    C         -1.425060         1.968520        -2.563470
    C          1.428170         1.968460        -2.563350
    C         -1.157810        -2.798820        -1.810080
    C          1.160670        -2.798860        -1.809990
    C         -2.306870        -1.963940        -1.810130
    C          2.309760        -1.964030        -1.809950
    C         -3.023270         0.241080        -1.810190
    C          3.026260         0.240960        -1.809950
    C         -2.584340         1.591890        -1.810190
    C          2.587370         1.591790        -1.809980
    C         -0.708620         2.954620        -1.810130
    C          0.711720         2.954590        -1.810080
    C         -0.715080        -3.408260        -0.591160
    C          0.717820        -3.408290        -0.591100
    C         -3.023360        -1.731120        -0.591270
    C          3.026170        -1.731240        -0.591030
    C         -3.466120        -0.368350        -0.591300
    C          3.468990        -0.368480        -0.591030
    C         -2.584370         2.345230        -0.591300
    C          2.587350         2.345130        -0.591100
    C         -1.425110         3.187440        -0.591270
    C          1.428120         3.187380        -0.591160
    C         -1.425290        -3.177480         0.617020
    C          1.427950        -3.177540         0.617130
    C         -2.584510        -2.335220         0.616970
    C          2.587210        -2.335320         0.617170
    C         -3.466150         0.378390         0.616900
    C          3.468960         0.378250         0.617170
    C         -3.023340         1.741140         0.616900
    C          3.026200         1.741030         0.617140
    C         -0.714980         3.418190         0.616970
    C          0.717910         3.418160         0.617030
    C         -0.708880        -2.944690         1.835940
    C          0.711460        -2.944720         1.836000
    C         -2.584540        -1.581890         1.835850
    C          2.587170        -1.581990         1.836050
    C         -3.023420        -0.231050         1.835820
    C          3.026110        -0.231170         1.836060
    C         -2.306930         1.973930         1.835820
    C          2.309710         1.973840         1.836000
    C         -1.157830         2.808760         1.835850
    C          1.160650         2.808720         1.835940
    C         -1.425340        -1.958560         2.589220
    C          1.427900        -1.958610         2.589330
    C         -2.306980         0.755050         2.589150
    C          2.309660         0.754960         2.589340
    C          0.001370         2.432090         2.589220
    C         -0.715180        -0.981100         3.335950
    C          0.717720        -0.981130         3.336010
    C         -1.157940         0.381680         3.335920
    C          1.160540         0.381630         3.336010
    C          0.001320         1.223890         3.335950
   |---------------------------------------------------------------------------
   Molecular polarizability (Angstrom**3):
    75.07245027     0.00003655     0.00007911
     0.00003655    75.07245069     0.00002885
     0.00007911     0.00002885    75.07238398
   Mean:     75.07   (Angstrom**3)
             11.12   (a.u.)
   | ---------------------------------------------------------------------------
```

The parameters used in this case is taken from P+Q ISO [R] model from Ref. [4].

# References

1. Applequist, J., Carl, J. R., & Fung, K. K. (1972). An Atom Dipole Interaction Model for Molecular Polarizability. Application to Polyatomic Molecules and Determination of Atom Polarizabilities. Journal of the American Chemical Society, 94(9), 2952–2960. https://doi.org/10.1021/ja00764a010
2. Olson, M. L., & Sundberg, K. R. (1978). An atom monopole-dipole interaction model with charge transfer for the treatment of polarizabilities of π-bonded molecules. The Journal of Chemical Physics, 69(12), 5400–5404. https://doi.org/10.1063/1.436570
3. Applequist, J. (1993). Atom charge transfer in molecular polarizabilities. Application of the Olson-Sundberg model to aliphatic and aromatic hydrocarbons. Journal of Physical Chemistry, 97(22), 6016–6023. https://doi.org/10.1021/j100124a039
4. Mayer, A. (2007). Formulation in terms of normalized propagators of a charge-dipole model enabling the calculation of the polarization properties of fullerenes and carbon nanotubes. Physical Review B, 75(4), 045407. https://doi.org/10.1103/PhysRevB.75.045407
5. Jensen, L. L., & Jensen, L. (2008). Electrostatic Interaction Model for the Calculation of the Polarizability of Large Noble Metal Nanoclusters. The Journal of Physical Chemistry C, 112(40), 15697–15703. https://doi.org/10.1021/jp804116z
6. Jensen, L. L., & Jensen, L. (2009). Atomistic Electrodynamics Model for Optical Properties of Silver Nanoclusters. The Journal of Physical Chemistry C, 113(34), 15182–15190. https://doi.org/10.1021/jp904956f
