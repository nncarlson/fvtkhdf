# RectilinearGrid and StructuredGrid integration

Implementation handoff updated 2026-11-04, branch `structured`.

## Settled scope and architecture

The first implementation adds standalone RectilinearGrid and StructuredGrid
writers with static or fixed-geometry temporal data. Moving geometry and
MultiBlock integration are deferred. Preserve an architecture that can later
place each concrete dataset type in a composite without using file wrappers.

The internal hierarchy is:

```text
vtkhdf_structured_base                abstract, vtkhdf_structured_base_class
├── vtkhdf_image                      vtkhdf_image_type
├── vtkhdf_rectilinear                vtkhdf_rectilinear_type
└── vtkhdf_structured                 vtkhdf_structured_type
```

The base owns the dataset group, logical grid bounds and mapping, point/cell
registries and storage, field data, and time stepping. Geometry belongs to the
concrete types. The file wrappers own the file and communicator; datasets retain
an unowned communicator copy. The three wrappers share forwarding code through
`vtkhdf_structured_file_template.fypp`. No deferred geometry interface is needed because
the concrete geometry signatures differ.

Future MB support may introduce an abstract `vtkhdf_dataset` ancestor shared
with UG and an owning `class(vtkhdf_dataset), allocatable` component in each
block. That hierarchy and MB dispatch are not part of this increment. Future
moving geometry will need separate geometry snapshots/offsets and lifecycle
hooks; do not encode it in unrelated concrete extensions now.

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

## References

- [VTKHDF specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html).
- [VTK 9.7.1 reader](https://github.com/Kitware/VTK/blob/v9.7.1/IO/HDF/vtkHDFReader.cxx).
- [VTK reader implementation](https://github.com/Kitware/VTK/blob/v9.7.1/IO/HDF/vtkHDFReaderImplementation.cxx).
- [GNU packing reproducer](https://github.com/nncarlson/fortran-compiler-tests/blob/master/gfortran-bugs/gfortran-20261103.f90).
