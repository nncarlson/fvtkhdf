# GNU polymorphic array bounds reproducer

`gnu_polymorphic_bounds.f90` is a standalone 13-line program. Its expected
output is `42`: the assumed-shape dummy `x(:)` has bounds `1:2`, so `x(1)`
denotes the caller's `a(0)`. GNU instead prints `99`, the value of `a(1)`.
There is no inheritance, allocation, MPI, HDF5, or library dependency.

Results on 2026-10-27, tested at both `-O0` and `-O2`:

| Compiler | Additional flags | Output at both optimization levels |
| --- | --- | --- |
| GNU Fortran 15.2.0 | `-fcheck=all` | `99` (incorrect) |
| GNU Fortran 16.1.1 20260517 | `-fcheck=all` | `99` (incorrect) |
| NAG 7.2, build 7245 | `-C=all` | `42` |
| Intel ifx 2026.0.0 20260331 | `-check all` | `42` |
| LLVM Flang 22.1.5 | none | `42` |

GNU 16.1.1 also prints `99` without runtime checking at both optimization
levels. The two-element case deliberately keeps the incorrectly indexed
element within the actual array's bounds, exposing the wrong value without
relying on a crash.

Two independent control variants print `42` with all five compiler versions
at both optimization levels, using the same checking flags:

- Change the actual declaration from `a(0:1)` to `a(1:2)`.
- Change the dummy declaration from `class(t)` to `type(t)`.

To reproduce, load the appropriate compiler module, choose one compilation
command below, and run `./repro`. Run from this directory; replace `-O0` with
`-O2` to check optimized compilation.

```sh
gfortran -O0 -fcheck=all gnu_polymorphic_bounds.f90 -o repro
nagfor -O0 -C=all gnu_polymorphic_bounds.f90 -o repro
ifx -O0 -check all gnu_polymorphic_bounds.f90 -o repro
flang -O0 gnu_polymorphic_bounds.f90 -o repro
./repro
```

The compiler modules used were `gcc/15.2.0`, `gcc/16.1.1`, `nag/7.2`,
`intel/2026.0`, and `llvm/22.1.5`. Executables were built and run in temporary
directories. NAG required execution outside the sandbox for licence validation.
