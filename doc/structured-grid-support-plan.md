# RectilinearGrid and StructuredGrid integration

Implementation handoff updated 2026-11-05, branch `structured`.
Work in progress: implementation and documentation review continue. The initial
new-format implementation was committed as `b1b314e`; this checkpoint includes
the subsequent file/core merger and runtime-check policy changes.

## Settled scope and architecture

The first implementation adds standalone RectilinearGrid and StructuredGrid
writers with static or fixed-geometry temporal data. Moving geometry and
public MultiBlock support for the new types remain deferred.

The implementation now merges each dataset core into its public file type:

```text
vtkhdf_simple_dataset                 abstract, vtkhdf_simple_dataset_class
├── vtkhdf_ug_file                    vtkhdf_ug_file_type
└── vtkhdf_structured_base            abstract, vtkhdf_structured_base_class
    ├── vtkhdf_image_file             vtkhdf_image_file_type
    ├── vtkhdf_rectilinear_file       vtkhdf_rectilinear_file_type
    └── vtkhdf_structured_file        vtkhdf_structured_file_type
```

The simple ancestor specifies only the shared lifecycle contract: close, root,
start/finalize time step, temporal-enabled and step-open queries. The concrete
branches retain their implementation and state. Structured geometry interfaces
differ and do not belong to the ancestor. Numeric overload generation remains
Fypp-based, but the shared file-wrapper include and forwarding layer are gone.

The structured base declares deferred INIT and SUPPORTS_MODE bindings. Each
concrete file type supplies its VTK type string to the shared INIT_GRID method
and defines its permitted modes. All three currently permit static and fixed
mesh; future moving-geometry support can have a different policy per type.

CREATE owns a standalone file and duplicated communicator through a
vtkhdf_file_context component. Each public file type also has an undocumented
internal INIT entry point taking a parent group and borrowed context for MB.
Close releases dataset resources, closes any owned file and communicator, then
resets the object. Implicit finalization performs no collective cleanup.

Handle creation and owner validation reside in the implementation that owns the
registry. MB sets/checks the existing block-index component and forwards the
handle to its leaf. MB stores leaves as class(vtkhdf_simple_dataset), allocatable,
using polymorphic lifecycle calls and type selection for format-specific calls.
Its public block creation remains UG-only in this increment. The nested-leaf
regression exercises all four concrete types through the ancestor independently
of that public API.

Future moving geometry will need separate geometry snapshots/offsets and
lifecycle hooks; do not encode it in unrelated concrete extensions now.

## Agreed public geometry interfaces

```fortran
call image%write_mesh(point_bounds, origin, spacing [, direction])
call rectilinear%write_mesh(point_bounds, x, y, z)
call structured%write_mesh(point_bounds)
call structured%write_mesh_geometry(points, array_origin [, write_bounds])
```

Rectilinear calls remain collective. Rank zero supplies complete coordinate
vectors; their values and lengths are ignored on other ranks, where empty
arrays suffice. All arrays/ranks use the same real kind, real32 or real64.
Singleton axes still have one coordinate. Distributed axis segments are deferred.

Structured coordinates have shape `(3,nx,ny,nz)` and support the same rectangular
region selection, halo exclusion, repeated transfers, and empty contributions
as point data. The first geometry write establishes the stored real kind;
subsequent calls must match. This choice avoids adding a kind argument to
`write_mesh`. Geometry must be supplied before time stepping.

## Logical indices and file representation

All three writers retain arbitrary user logical indices in their APIs:

- Inclusive `point_bounds(3,2)`: lower bounds in column 1, upper bounds in column 2.
- `array_origin(3)`: global logical index of the first spatial array element,
  independent of Fortran subscripts.
- Optional inclusive `write_bounds(3,2)` in the same logical space; it must fit
  within global association bounds and the array's mapped footprint.
- Scalar data `(nx,ny,nz)`, vector data `(ncomp,nx,ny,nz)`, including singleton axes.
- Cells indexed by their lower-index corner; singleton directions keep their index.
- Empty contributions use a zero spatial extent or empty write bounds.
  `array_origin` remains required but need not identify an in-grid index then.

The storage mapping subtracts global lower bounds. ImageData stores WholeExtent
and its Origin/Spacing/Direction attributes. The two new types store Dimensions;
readers expose zero-based extents while physical coordinates/data are preserved.
Rectilinear coordinates are 1D XCoordinates/YCoordinates/ZCoordinates datasets.
Structured Points has HDF5 shape `(nz,ny,nx,3)`, retaining singleton axes. Point
and cell datasets squeeze inactive spatial axes as before.

Fixed geometry has no temporal dimension or geometry offset datasets, even when
associated data is temporal. Existing data registration, snapshots, and reuse
apply unchanged. Real data and coordinates use quiet NaN fill; integer data uses
zero. Hidden points/cells use vtkGhostType flags 2/32 respectively.

## Validation and related cautions

Regression coverage includes both real kinds, all active-axis combinations,
nonzero logical bounds, coordinate ordering, noncontiguous halo arrays, repeated
regions, empty contributions, temporal reuse, lifecycle errors, and installed
consumers. Check layout with h5py and physical geometry/data independently with
VTK. Readers before 9.7.2 have a shared structured-data bug for nontrailing
singleton axes; reader tests skip only those cases while layout checks remain
active. ImageData additionally has the positive-global-lower-bound reader bug.
The user confirmed both fixes were backported for VTK 9.7.2.

The prior CI problems are resolved. Preserve GNU_BUG20261103 (GNU before 16.1)
and HDF5_NULL_MEMSPACE_BUG (parallel HDF5 before 1.13). Existing CI images remain
the baseline. The shared exported vtkhdf_version parameter is [2,8], matching
the targeted specification in README.

The [runtime checking policy](runtime-checks.md) records the error-handling
philosophy and how to choose checks. Agreement checks now use debug-only ASSERT;
argument validity checks use local INSIST. Collective library-success checks and
algorithmic reductions remain unconditional. Method-identification tags and
empty-rank padding collectives were removed. The three MPI disagreement tests
are Debug-only.

Latest GNU 16.1.1 results with HDF5 1.14.6 and MPICH 5.0.1:

| Configuration | Serial passed | MPI passed | Reader skips in each suite |
| --- | ---: | ---: | ---: |
| Debug | 85 | 95 | 9 |
| Release | 85 | 92 | 9 |

Release object inspection confirmed no references to the agreement helpers from
the structured writers and HDF5 support. The merger also passed GNU 14.2 and
LLVM 22 serial tests, installed serial/MPI consumers, a clean build, and the
Sphinx build. NAG 7.2 serial and MPI tests passed before the subsequent checking
refinements. NAG compilation needs escalated access to its licence. For a NAG
parallel build, configure a fresh directory with FC=mpifort and CC=mpicc after
loading `ml nag mpich hdf5`; selecting nagfor directly causes linker problems.
fpm source generation was checked, but fpm itself was unavailable locally.

## References

- [VTKHDF specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html).
- [VTK 9.7.1 reader](https://github.com/Kitware/VTK/blob/v9.7.1/IO/HDF/vtkHDFReader.cxx).
- [VTK reader implementation](https://github.com/Kitware/VTK/blob/v9.7.1/IO/HDF/vtkHDFReaderImplementation.cxx).
- [GNU packing reproducer](https://github.com/nncarlson/fortran-compiler-tests/blob/master/gfortran-bugs/gfortran-20261103.f90).
