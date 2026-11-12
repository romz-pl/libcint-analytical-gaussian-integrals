# libcint: An Analytical Gaussian Integral Library

**libcint** is an open-source C library for evaluating analytical integrals over Gaussian-type orbitals (GTOs), the basis functions used in most molecular electronic-structure calculations. It was written by Qiming Sun and collaborators, originally to serve the PySCF quantum chemistry framework, and it is distributed on GitHub under a permissive BSD-style license. Its main purpose is to compute integral batches quickly and accurately, leaving higher-level algorithms such as Hartree-Fock, density functional theory, and coupled cluster to the calling code.

## Scope of the Integrals

libcint covers the integral types needed for molecular calculations:

- **One-electron integrals**: overlap, kinetic energy, nuclear attraction (including point-charge potentials), and multipole and dipole-type operators.
- **Two-electron repulsion integrals (ERIs)**: the four-center integrals (ij|kl) that dominate the cost of Hartree-Fock and many correlated methods. These are the main performance target of the library.
- **Auxiliary integrals**: two-center (2c2e) and three-center (3c2e) Coulomb integrals, which are central to density fitting and resolution-of-the-identity methods.
- **Range-separated and attenuated operators**: Coulomb interactions modified by an error-function or similar attenuation, used in screened hybrid functionals and long-range corrected methods.
- **Derivative integrals**: first and higher derivatives with respect to nuclear coordinates, needed for analytic gradients and Hessians.
- **Spinor integrals**: versions in a two-component spinor representation, for relativistic calculations.

The integrals are named systematically. The name combines an integral class (such as `int1e`, `int2e`, `int2c2e`, or `int3c2e`) with a suffix that encodes a derivative (for example `_ip1` for the first derivative with respect to the first center) or a representation (`_sph` for spherical harmonics, `_cart` for Cartesian functions, `_spinor` for spinors). Calling code builds the name it needs and passes it to a generic dispatch function.

## Mathematical Approach

The core difficulty is evaluating integrals over products of Gaussians with arbitrary angular momentum efficiently. libcint relies on well-established recurrence methods:

- **Obara-Saika and Head-Gordon-Pople recurrences** reduce high-angular-momentum integrals to simpler base integrals through horizontal and vertical recurrence relations. This lets the library handle arbitrary contraction and high angular momentum without deriving formulas by hand for each case.
- **Rys quadrature** is the main approach for the electron repulsion integrals. It rewrites the integrals as weighted sums over a small number of quadrature roots, which is numerically stable and efficient for the four-center case.
- **Boys function evaluation** supplies the auxiliary function values that underlie many of these expressions, with tabulation and interpolation used to keep evaluation fast.

Contracted basis functions are handled by summing over primitive Gaussians, and the contraction coefficients are applied in the final accumulation.

## Code Generation and Performance

Because the number of angular momentum combinations is large, libcint's kernels are produced by a code generator rather than written manually. The generator emits specialized C code for each combination of angular momenta, which is then compiled into the library. This approach gives:

- **Fast inner loops** with the recurrence coefficients precomputed, avoiding branches on angular momentum.
- **Vectorization-friendly code** that compilers can optimize well.
- **Multi-threading support** through OpenMP, so integral batches can be evaluated in parallel.
- **Screening support** through optimizer objects (such as those created by the two-electron optimizer initialization), which estimate integral magnitudes and let the caller skip negligible contributions.

## The Interface

libcint's C interface is deliberately low-level. A molecule is described by three flat arrays, a convention shared with PySCF:

- `atm`, an integer array of shape (natm, 6) describing each atom, including its charge and the offset of its coordinates in `env`.
- `bas`, an integer array of shape (nbas, 8) describing each shell: its atom, angular momentum, number of primitives, number of contractions, and offsets into `env` for exponents and coefficients.
- `env`, a flat double array holding coordinates, exponents, contraction coefficients, and other numerical data.

An integral is requested by passing a shell index list (`shls`), the output buffer, the dimensions, these three arrays, an optimizer, and a scratch cache. Helper functions such as `CINTcgto_spheric` and `CINTcgto_cart` return the number of basis functions in a shell, which the caller needs to size buffers correctly.

This design keeps the library independent of any particular molecule or program. It does not own the geometry or basis set; the calling code does, and libcint simply evaluates what it is asked for.

## Relationship to PySCF

In practice, most users encounter libcint through PySCF, which wraps it in Python. Calls such as `mol.intor('int1e_ovlp')` or `mol.intor('int2e')` dispatch to libcint behind the scenes. PySCF also uses it for gradients, auxiliary-basis integrals, and many other operations. This makes libcint a central dependency of the PySCF ecosystem, and its development has been driven largely by the needs of that framework.

## Scope Limits

A few boundaries are worth knowing:

- libcint handles Gaussian-type orbitals. Effective core potential (ECP) integrals are not part of its core functionality; they are typically provided by a separate library, libecpint.
- It is an integral library, not a complete electronic-structure code. Orbital optimization, Fock matrix construction, and correlation methods must be implemented by the caller.
- The high-level API is C and requires the caller to manage the atomic and basis-set arrays, which is more involved than a Python interface but gives full control over memory and performance.

## Building and Use

The library is built with CMake and requires only a C compiler, with OpenMP as an optional dependency for parallel evaluation. It has no mandatory external numerical libraries, which keeps it easy to embed in other projects.

## Summary

libcint is a focused, high-performance engine for analytical Gaussian integrals. Its strengths are broad integral coverage, support for arbitrary angular momentum through recurrence relations and generated kernels, screening and parallelism, and a lean C interface that fits naturally into larger quantum chemistry codes. If you are working on a method that needs many integrals quickly, or you want to understand how PySCF computes its integrals under the hood, libcint is the component to study. If you tell me which integral type or part of the interface you're working with, I can go into more detail on that.
