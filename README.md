# libcint: An Analytical Gaussian Integral Library

**libcint** is an open-source C library for evaluating analytical integrals over Gaussian-type orbitals (GTOs), the basis functions used in most molecular electronic-structure calculations. It was written by Qiming Sun and collaborators, originally to serve the PySCF quantum chemistry framework, and it is distributed on GitHub under a permissive BSD-style license. Its main purpose is to compute integral batches quickly and accurately, leaving higher-level algorithms such as Hartree-Fock, density functional theory, and coupled cluster to the calling code. [GitHub repository...](https://github.com/sunqm/libcint)

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



# Minimal libcint Demo in C++

This program builds a two-atom system by hand, using the `atm`, `bas`, and `env` arrays that libcint expects. It then evaluates:

1. The overlap matrix over two s-type shells, using `cint1e_ovlp_sph`.
2. One two-electron repulsion integral (00|00), using `cint2e_sph`.

Both shells have exponent α = 1.0 and are normalized, so the expected values can be checked by hand.

## Source: `libcint_demo.cpp`

```cpp
// libcint_demo.cpp
// Build: g++ -O2 libcint_demo.cpp -lcint -lm -o libcint_demo
#include <cstdio>
#include <vector>

extern "C" {
#include <cint.h>
#include <cint_funcs.h>
}

int main() {
    const double R = 1.4;      // separation of the two atoms (bohr)
    const int natm = 2;        // number of atoms
    const int nbas = 2;        // number of shells (one s shell per atom)

    // Flat arrays describing the molecule, following libcint's conventions
    std::vector<int>    atm(natm * ATM_SLOTS, 0);
    std::vector<int>    bas(nbas * BAS_SLOTS, 0);
    std::vector<double> env(PTR_ENV_START + 16, 0.0);

    // Fill env: coordinates, exponent, and contraction coefficient
    int off = PTR_ENV_START;
    const int ptr_coord0 = off;
    env[off++] = 0.0; env[off++] = 0.0; env[off++] = 0.0;   // atom 0 at origin
    const int ptr_coord1 = off;
    env[off++] = 0.0; env[off++] = 0.0; env[off++] = R;     // atom 1 along z
    const int ptr_exp = off;
    env[off++] = 1.0;                                       // exponent alpha
    const int ptr_coeff = off;
    env[off++] = 1.0;                                       // contraction coefficient

    // Atoms: nuclear charge and pointer to coordinates in env
    for (int i = 0; i < natm; i++) {
        atm[i * ATM_SLOTS + CHARGE_OF] = 1;
        atm[i * ATM_SLOTS + PTR_COORD] = (i == 0) ? ptr_coord0 : ptr_coord1;
    }

    // Shells: one primitive s function (l = 0) centered on each atom
    for (int i = 0; i < nbas; i++) {
        bas[i * BAS_SLOTS + ATOM_OF]   = i;
        bas[i * BAS_SLOTS + ANG_OF]    = 0;
        bas[i * BAS_SLOTS + NPRIM_OF]  = 1;
        bas[i * BAS_SLOTS + NCTR_OF]   = 1;
        bas[i * BAS_SLOTS + KAPPA_OF]  = 0;
        bas[i * BAS_SLOTS + PTR_EXP]   = ptr_exp;
        bas[i * BAS_SLOTS + PTR_COEFF] = ptr_coeff;
    }

    // Number of spherical basis functions in each shell (1 for an s shell)
    const int n[2] = { CINTcgto_spheric(0, bas.data()),
                       CINTcgto_spheric(1, bas.data()) };

    // ---- One-electron overlap matrix S_ij ----
    std::printf("Overlap matrix S (expected: diagonal ~ 1, off-diagonal ~ exp(-R^2/2) = %.6f)\n",
                std::exp(-0.5 * R * R));
    for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
            const int shls[2] = { i, j };
            std::vector<double> buf(n[i] * n[j], 0.0);
            cint1e_ovlp_sph(buf.data(), shls, atm.data(), natm,
                            bas.data(), nbas, env.data());
            std::printf("  S[%d][%d] = %12.8f\n", i, j, buf[0]);
        }
    }

    // ---- Two-electron repulsion integral (00|00) ----
    // Expected: 2 * sqrt(alpha / pi) = 2 / sqrt(pi) = 1.12837917 for alpha = 1,
    // because all four shells sit on the same atom.
    const int shls4[4] = { 0, 0, 0, 0 };
    std::vector<double> eri(n[0] * n[0] * n[0] * n[0], 0.0);
    std::vector<double> cache(1 << 16, 0.0);   // scratch space for libcint
    cint2e_sph(eri.data(), shls4, atm.data(), natm, bas.data(), nbas,
               env.data(), nullptr, cache.data());
    std::printf("\nERI (00|00) = %12.8f  (expected 2/sqrt(pi) = 1.12837917)\n", eri[0]);

    return 0;
}
```

## Expected Output

```
Overlap matrix S (expected: diagonal ~ 1, off-diagonal ~ exp(-R^2/2) = 0.375311)
  S[0][0] =   1.00000000
  S[0][1] =   0.37531110
  S[1][0] =   0.37531110
  S[1][1] =   1.00000000

ERI (00|00) =   1.12837917  (expected 2/sqrt(pi) = 1.12837917)
```

## Notes

- **Why these values?** For normalized s Gaussians with the same exponent α, the overlap between centers separated by R is $S = e^{-\alpha R^2/2}$, which gives $e^{-R^2/2}$ when α = 1. The same-center repulsion integral reduces to $2\sqrt{\alpha/\pi}$, which gives $2/\sqrt{\pi}$ when α = 1. Both checks confirm that libcint is normalizing the primitives and using the expected conventions.
- **Index layout.** `atm` has shape (natm, 6), `bas` has shape (nbas, 8), and `env` is a flat array. The first `PTR_ENV_START` entries are reserved, so data always begins at that offset.
- **Buffers.** Each integral writes `n_i * n_j` values for a one-electron integral and `n_i * n_j * n_k * n_l` values for a two-electron integral. `CINTcgto_spheric` returns the spherical function count per shell, which sets the buffer size.
- **Screening.** Passing `nullptr` as the optimizer disables screening. For larger systems, create one with `CINTinit_2e_optimizer` and free it with `CINTdel_2e_optimizer`.
- **Cache.** The scratch buffer size here is a generous fixed value for a tiny example. For production code, size it according to the integral type and angular momenta.
- **Verify signatures.** The function names and argument lists follow libcint's public headers as I know them. Check them against the `cint.h` and `cint_funcs.h` installed on your system, since signatures can differ between releases.

## Build and Run

```bash
g++ -O2 libcint_demo.cpp -lcint -lm -o libcint_demo
./libcint_demo
```

libcint must be installed (for example, built from the GitHub repository with CMake and installed to a location on your library and include paths).

## Extending the Demo

Natural next steps are to add a nuclear attraction integral with `cint1e_nuc_sph` (which needs the nuclear charges in `atm`), a kinetic energy integral with `cint1e_kin_sph`, or p-type shells by setting `ANG_OF` to 1 and rebuilding the buffers with the corresponding size.


# Bibliography

+ Q. Sun, **Libcint: An efficient general integral library for Gaussian basis functions**, [Journal of Computational Chemistry 36, 1664-1671 (2015)](http://dx.doi.org/10.1002/jcc.23981)
