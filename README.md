This fortran demo program implement charge-dipole interaction model with Gaussian-type density functions, i.e., Q+P ISO [R] in Ref. [4].
See [this link](https://link.aps.org/doi/10.1103/PhysRevB.75.045407) for more details.

# Installation

One needs to install this program first by running

```console
make
```
and it is assumed that `gfortron` compiler and libraries `llapack` and `lblas` have been installed.

# Usage

Run `./main.x {xyz_file}` to compute the molecular dipole polarizability and `xyz_file` is molecular geometry file with `xyz` format, e.g. `C60.xyz`:

```console
./main.x C60.xyz
```

One should obtain the output as follows:

```console
|---------------------------Input Geometry-----------------------------------------
All coordinates are in unit of Angstrom.
Number of atoms   60

 C          1.253859         0.000000         0.000000
 C          0.387464        -1.192491         0.000000
 C          0.387464         1.192491         0.000000
 C         -1.014393        -0.737000         0.000000
 C         -1.014393         0.737000         0.000000
 C          2.444770         0.000000         0.736023
 C          0.755476        -2.325115         0.736023
 C          0.755476         2.325115         0.736023
 C         -1.977861        -1.437000         0.736023
 C         -1.977861         1.437000         0.736023
 C          2.832234        -1.192491         1.510951
 C          2.832234         1.192491         1.510951
 C          2.009335        -2.325115         1.510951
 C          2.009335         2.325115         1.510951
 C         -0.258918        -3.062115         1.510951
 C         -0.258918         3.062115         1.510951
 C         -1.590397        -2.629491         1.510951
 C         -1.590397         2.629491         1.510951
 C         -2.992254        -0.700000         1.510951
 C         -2.992254         0.700000         1.510951
 C          3.459164        -0.737000         2.764810
 C          3.459164         0.737000         2.764810
 C          1.769869        -3.062115         2.764810
 C          1.769869         3.062115         2.764810
 C          0.368012        -3.517606         2.764810
 C          0.368012         3.517606         2.764810
 C         -2.365325        -2.629491         2.764810
 C         -2.365325         2.629491         2.764810
 C         -3.231720        -1.437000         2.764810
 C         -3.231720         1.437000         2.764810
 C          3.231720        -1.437000         3.955722
 C          3.231720         1.437000         3.955722
 C          2.365325        -2.629491         3.955722
 C          2.365325         2.629491         3.955722
 C         -0.368012        -3.517606         3.955722
 C         -0.368012         3.517606         3.955722
 C         -1.769869        -3.062115         3.955722
 C         -1.769869         3.062115         3.955722
 C         -3.459164        -0.737000         3.955722
 C         -3.459164         0.737000         3.955722
 C          2.992254        -0.700000         5.209581
 C          2.992254         0.700000         5.209581
 C          1.590397        -2.629491         5.209581
 C          1.590397         2.629491         5.209581
 C          0.258918        -3.062115         5.209581
 C          0.258918         3.062115         5.209581
 C         -2.009335        -2.325115         5.209581
 C         -2.009335         2.325115         5.209581
 C         -2.832234        -1.192491         5.209581
 C         -2.832234         1.192491         5.209581
 C          1.977861        -1.437000         5.984509
 C          1.977861         1.437000         5.984509
 C         -0.755476        -2.325115         5.984509
 C         -0.755476         2.325115         5.984509
 C         -2.444770         0.000000         5.984509
 C          1.014393        -0.737000         6.720532
 C          1.014393         0.737000         6.720532
 C         -0.387464        -1.192491         6.720532
 C         -0.387464         1.192491         6.720532
 C         -1.253859         0.000000         6.720532
|---------------------------------------------------------------------------
Molecular polarizability:
 78.43810985     0.00000000     0.00000032
 -0.00000000    78.43811711     0.00000000
  0.00000032    -0.00000000    78.43811037
| ---------------------------------------------------------------------------
 Note: The following floating-point exceptions are signalling: IEEE_DENORMAL
```

# References

1. Applequist, J., Carl, J. R., & Fung, K. K. (1972). An Atom Dipole Interaction Model for Molecular Polarizability. Application to Polyatomic Molecules and Determination of Atom Polarizabilities. Journal of the American Chemical Society, 94(9), 2952–2960. https://doi.org/10.1021/ja00764a010
2. Olson, M. L., & Sundberg, K. R. (1978). An atom monopole-dipole interaction model with charge transfer for the treatment of polarizabilities of π-bonded molecules. The Journal of Chemical Physics, 69(12), 5400–5404. https://doi.org/10.1063/1.436570
3. Applequist, J. (1993). Atom charge transfer in molecular polarizabilities. Application of the Olson-Sundberg model to aliphatic and aromatic hydrocarbons. Journal of Physical Chemistry, 97(22), 6016–6023. https://doi.org/10.1021/j100124a039
4. Mayer, A. (2007). Formulation in terms of normalized propagators of a charge-dipole model enabling the calculation of the polarization properties of fullerenes and carbon nanotubes. Physical Review B, 75(4), 045407. https://doi.org/10.1103/PhysRevB.75.045407
5. Jensen, L. L., & Jensen, L. (2008). Electrostatic Interaction Model for the Calculation of the Polarizability of Large Noble Metal Nanoclusters. The Journal of Physical Chemistry C, 112(40), 15697–15703. https://doi.org/10.1021/jp804116z
6. Jensen, L. L., & Jensen, L. (2009). Atomistic Electrodynamics Model for Optical Properties of Silver Nanoclusters. The Journal of Physical Chemistry C, 113(34), 15182–15190. https://doi.org/10.1021/jp904956f
