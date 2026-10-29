# VTK 9.7.1 VTKHDF ImageData reader failures

Verified 22 October 2026.

## Findings and upstream status

Two geometry/array-readback defects are reproducible with VTK 9.7.1:

1. **Positive nonzero ImageData extents cause out-of-range HDF5 reads.**
2. **Singleton axes other than a trailing z axis are handled inconsistently.**
   Some cases fail with a dataset-rank error; a singleton y axis can silently
   truncate PointData and return empty CellData.

The standalone reproducer below writes small files directly with h5py. It
requires neither fVTKHDF nor Fortran, MPI, or a simulation application. PointData
and CellData are tested in separate files, so a CellData failure is observed
independently of any preceding PointData failure.

All 18 geometry/association checks pass using the same reproducer with a VTK
development build reporting `9.7.20261022`, built from
[commit 9f0588e2fe0d10573555834a66ac3dab95be4114](https://github.com/Kitware/VTK/commit/9f0588e2fe0d10573555834a66ac3dab95be4114).
The built reader source was checked against that pinned commit. Both geometry
defects are fixed in that development snapshot. Both fixes have since been
backported to the VTK 9.7 release branch and are queued for VTK 9.7.2. The
recorded tests used VTK 9.7.1 and the development snapshot above; no ParaView
release was tested.

**Correction concerning temporal FieldData:** the earlier investigation also
found truncation of fVTKHDF's former flattened multicomponent temporal FieldData stream.
Additional controls show that VTK 9.7.1 correctly reads the corresponding
tuple-by-component representation. That observation should be treated as a
fVTKHDF layout interoperability issue, rather than included as a third confirmed
VTK reader defect. Details are provided separately below.

## Current user constraint and project policy

For VTK 9.7.0/9.7.1 or ParaView builds using the affected reader, use zero-based
global point bounds on all axes and trailing singleton axes: `(nx,ny,1)` for 2D
or `(nx,1,1)` for 1D. In particular, one-based whole extents are affected by the
positive-lower-bound defect. Positive local block indices and nonzero physical
Origin values remain usable with zero-based global bounds.

fVTKHDF preserves the supplied global bounds and physical origin. It will not
add a compatibility option or silently rebase the stored extent. The geometry
constraints can be lifted when using a VTK or ParaView build containing the
verified upstream fixes. The public API continues to support general logical
bounds and singleton-axis choices.

## Environment

| Component | Measured value |
| --- | --- |
| Python interpreter | `/usr/bin/python3` |
| Python | 3.14.8 |
| Operating system | Linux x86-64, kernel `7.2.5-100.fc43.x86_64`, glibc 2.42 |
| Installed VTK | `vtkVersion.GetVTKVersion() == "9.7.1"` |
| Installed VTK source version | `vtk version 9.7.1` |
| NumPy | 2.3.5 |
| h5py | 3.13.0 |
| HDF5 used by h5py to write fixtures | 1.14.6 |
| HDF5 version in VTK's read-error banner | 1.14.6 |
| Development comparison | VTK `9.7.20261022`, commit linked above |

All fixtures use VTKHDF `Version=[2,8]`, `Type="ImageData"`, identity Direction,
zero physical Origin, unit Spacing, and scalar numeric arrays. The failing cases
do not involve direction matrices, vector components, temporal mesh arrays,
partitioning, or compression.

## Defect 1: positive nonzero extents become absolute HDF5 offsets

Suggested issue title: **vtkHDFReader: ImageData arrays fail with positive
nonzero WholeExtent lower bounds**.

### Minimal case

| Property | Value |
| --- | --- |
| Logical point dimensions `(nx,ny,nz)` | `(4,3,2)` |
| `WholeExtent` | `[2,5, 3,5, 4,5]` |
| PointData dataset shape, in HDF5 order | `(2,3,4)` |
| CellData dataset shape, in HDF5 order | `(1,2,3)` |
| Expected scalar point tuples | 24 |
| Expected scalar cell tuples | 6 |

The point values are `0..23`, stored in z,y,x order. The cell values are `0..5`.
Moving the logical lower bounds must not change either dataset's storage size:
the lower corner of each entire array is at relative HDF5 offset `(0,0,0)`.

### Actual VTK 9.7.1 behavior

Both independent reads fail and the requested array is absent from the output.
Representative diagnostics, with timestamps omitted:

```text
selection + offset not within extent for file dataspace
Error H5Dread start: 4, 3, 2,  count: 2 3 4
```

For CellData, the corresponding diagnostic has the same start `(4,3,2)` and
count `(1,2,3)`. The start is the reversed **logical** lower bound, rather than
the required dataset-relative zero offset.

Two controls with the same dimensions and values pass:

| WholeExtent lower bounds | Point tuples | Cell tuples |
| --- | --- | --- |
| `(0,0,0)` | 24 | 6 |
| `(-2,-3,-4)` | 24 | 6 |

The negative-bound control is a full-domain read only; it establishes no general
claim about streamed subextents with negative global bounds.

An additional native-writer control reproduces the positive-bound failure:
VTK 9.7.1's own `vtkHDFWriter` successfully writes an image with this extent
(`Write()` returns 1), producing PointData shape `(2,3,4)` and CellData shape
`(1,2,3)`. Its `vtkHDFReader` then fails to read those arrays. Thus this case also
fails a native VTK writer/reader round trip.

### Source evidence and comparison

In [VTK 9.7.1's ReadStructuredData](https://github.com/Kitware/VTK/blob/v9.7.1/IO/HDF/vtkHDFReader.cxx#L299-L358),
the extent buffer retains logical indices. A negative selected lower bound is
shifted to zero, but a positive lower bound is left unchanged. The global lower
bound is not generally subtracted. This explains both the failing positive
case and the passing negative full-domain control.

The [tested development implementation](https://github.com/Kitware/VTK/blob/9f0588e2fe0d10573555834a66ac3dab95be4114/IO/HDF/vtkHDFReader.cxx#L272-L320)
computes each spatial selection relative to the corresponding global lower
bound. The positive-bound fixtures then read correctly, including every value.

## Defect 2: nontrailing singleton axes cause rank errors or silent truncation

Suggested issue title: **vtkHDFReader: ImageData singleton x/y axes produce
incorrect spatial selections and array lengths**.

### Effective restriction versus the ImageData model

For the compact layouts tested here, VTK 9.7.1 behaves as though singleton
axes must be trailing: `(nx,ny,nz)` for 3D, `(nx,ny,1)` for 2D, and `(nx,1,1)`
for 1D, with each active dimension greater than one. Thus 1D is supported along
logical x; prohibiting every singleton axis except z would overstate the
restriction. This is an implicit, unchecked assumption in this reader path.
The VTKHDF specification does not explicitly establish it as a format rule.

The underlying ImageData model supports a broader set of orientations:
[VTK 9.7.1's structured-data definitions](https://github.com/Kitware/VTK/blob/v9.7.1/Common/DataModel/vtkStructuredData.h#L45-L56)
include x-, y-, and z-directed lines and XY, YZ, and XZ planes. Direct calls to
`vtkImageData.GetDataDimension()` in 9.7.1 also return 1 for all three line
orientations and 2 for all three plane orientations. The report therefore
distinguishes the old VTKHDF reader's effective limitation from a general
ImageData dimensionality constraint; maintainers can clarify whether that
limitation was ever intentional.

### On-disk convention used by these fixtures

These fixtures omit globally singleton spatial axes from the attribute dataset.
This is the compact convention accepted by VTK 9.7.1 for ordinary XY images and
x-directed lines, and by the tested development reader for every axis choice.
Active cell axes are retained even if they contain only one cell.

The [published ImageData specification](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html#image-data)
describes spatial order as z,y,x and does not fully explain singleton-axis
storage. An upstream report should include the exact shapes below, making this
convention explicit and allowing maintainers to clarify the specification.

### Measured matrix

Every row uses zero logical lower bounds. Each association is read from its own
file. A missing array accompanies a reader error. Counts in the final two
columns are scalar-array tuple counts, rather than geometry counts.

| Case | Point dimensions `(x,y,z)` | PointData shape | CellData shape | Point tuples: expected → actual | Cell tuples: expected → actual |
| --- | --- | --- | --- | --- | --- |
| 3D control | `(4,3,2)` | `(2,3,4)` | `(1,2,3)` | 24 → 24 | 6 → 6 |
| XY plane; singleton z | `(4,3,1)` | `(3,4)` | `(2,3)` | 12 → 12 | 6 → 6 |
| YZ plane; singleton x | `(1,4,3)` | `(3,4)` | `(2,3)` | 12 → missing | 6 → missing |
| XZ plane; singleton y | `(4,1,3)` | `(3,4)` | `(2,3)` | 12 → 4 | 6 → 0 |
| x-directed line | `(4,1,1)` | `(4,)` | `(3,)` | 4 → 4 | 3 → 3 |
| y-directed line | `(1,4,1)` | `(4,)` | `(3,)` | 4 → missing | 3 → missing |
| z-directed line | `(1,1,4)` | `(4,)` | `(3,)` | 4 → missing | 3 → missing |

The singleton-x plane emits:

```text
Dataset: Expecting ndims >= 3, got: 2
```

The y- and z-directed lines emit:

```text
Dataset: Expecting ndims >= 2, got: 1
```

### Silent failure with singleton y

This case is especially significant because neither independent read emits a
reader `ErrorEvent`:

```text
WholeExtent = [0,3, 0,0, 0,2]
Geometry: 12 points, 6 cells
PointData expected: [0,1,2,3,4,5,6,7,8,9,10,11]
PointData actual:   [0,1,2,3]
CellData expected:  [0,1,2,3,4,5]
CellData actual:    []
```

Checking only the output class, geometry counts, or absence of reader errors
would miss the incorrect array lengths.

### Source evidence and comparison

[GetNDims and ReduceDimension in VTK 9.7.1](https://github.com/Kitware/VTK/blob/v9.7.1/IO/HDF/vtkHDFReader.cxx#L59-L86)
reduce dimensionality by inspecting z and y, then retain the first coordinate
pairs. They do not select the actual active axes. For an XZ plane, the retained
pairs describe x and singleton y, omitting active z. The subsequent point/cell
extent adjustment explains the four point values and zero cell values.

The development implementation linked above skips each globally singleton axis
while constructing its reversed spatial selection. All plane orientations and
all three line orientations then pass both tuple-count and value checks.

## Standalone reproducer

Save the following as `reproduce_vtkhdf_reader.py` and run it with the Python
interpreter whose VTK installation is being tested:

```bash
python3 -m pip install vtk==9.7.1 numpy h5py
python3 reproduce_vtkhdf_reader.py > results.txt 2> errors.txt
```

The script leaves small `.vtkhdf` files in `vtkhdf-reader-repro/`, suitable for
attaching to an issue together with `results.txt` and `errors.txt`. `FAIL` lines
are reported findings; the script deliberately completes the full matrix
without stopping at the first failed case.

```python
#!/usr/bin/env python3
"""Generate standalone VTKHDF fixtures and report vtkHDFReader results."""
import platform
from pathlib import Path

import h5py
import numpy as np
from vtkmodules.vtkCommonCore import vtkVersion
from vtkmodules.vtkIOHDF import vtkHDFReader
from vtkmodules.util.numpy_support import vtk_to_numpy

OUT = Path('vtkhdf-reader-repro')
OUT.mkdir(exist_ok=True)
print('Python:', platform.python_version(), 'Platform:', platform.platform())
print('VTK:', vtkVersion.GetVTKVersion())
print('NumPy:', np.__version__, 'h5py:', h5py.__version__,
      'h5py HDF5:', h5py.version.hdf5_version)


def header(f, lower, dims):
    g = f.create_group('VTKHDF')
    g.attrs['Version'] = np.array([2, 8], dtype=np.int32)
    g.attrs['Type'] = np.bytes_('ImageData')
    g.attrs['WholeExtent'] = np.array(
        [[lo, lo + n - 1] for lo, n in zip(lower, dims)],
        dtype=np.int32).ravel()
    g.attrs['Origin'] = np.zeros(3, dtype=np.float64)
    g.attrs['Spacing'] = np.ones(3, dtype=np.float64)
    g.attrs['Direction'] = np.eye(3, dtype=np.float64).ravel()
    for name in ['PointData', 'CellData', 'FieldData']:
        g.create_group(name)
    return g


def read(path, association, name, step=None):
    reader = vtkHDFReader()
    reader.SetFileName(str(path))
    errors = []
    reader.AddObserver('ErrorEvent', lambda *_: errors.append('ErrorEvent'))
    if step is not None:
        reader.UpdateInformation()
        reader.SetStep(step)
    reader.Update()
    image = reader.GetOutputDataObject(0)
    if image is None:
        return errors, None, None
    arrays = {'PointData': image.GetPointData,
              'CellData': image.GetCellData,
              'FieldData': image.GetFieldData}
    array = arrays[association]().GetArray(name)
    if array is None:
        return errors, image, None
    # Copy values while the reader and its output are still alive.
    return errors, image, (array.GetNumberOfComponents(),
                          array.GetNumberOfTuples(),
                          vtk_to_numpy(array).copy())


cases = [
    ('xyz_zero', [0, 0, 0], [4, 3, 2]),
    ('xyz_positive', [2, 3, 4], [4, 3, 2]),
    ('xyz_negative', [-2, -3, -4], [4, 3, 2]),
    ('xy_plane', [0, 0, 0], [4, 3, 1]),
    ('yz_plane', [0, 0, 0], [1, 4, 3]),
    ('xz_plane', [0, 0, 0], [4, 1, 3]),
    ('x_line', [0, 0, 0], [4, 1, 1]),
    ('y_line', [0, 0, 0], [1, 4, 1]),
    ('z_line', [0, 0, 0], [1, 1, 4]),
]
print('\nGEOMETRY: expected tuples -> actual tuples; status')
for case, lower, dims in cases:
    for association in ['PointData', 'CellData']:
        # Omit globally singleton point axes, retaining active cell axes even
        # when they have just one cell. Spatial order on disk is z,y,x.
        active = np.array(dims) > 1
        sizes = np.array(dims)
        if association == 'CellData':
            sizes = np.maximum(1, sizes - 1)
        shape = tuple(sizes[::-1][active[::-1]])
        count = int(np.prod(sizes))
        expected = np.arange(count, dtype=np.float64)
        path = OUT / f'{case}_{association}.vtkhdf'
        with h5py.File(path, 'w') as f:
            g = header(f, lower, dims)
            g[association]['values'] = expected.reshape(shape)
        errors, image, actual = read(path, association, 'values')
        correct = (actual is not None and actual[0] == 1 and
                   actual[1] == count and
                   np.array_equal(actual[2], expected) and not errors)
        observed = 'missing' if actual is None else str(actual[1])
        print(f'{case:14} {association:9} {count:2} -> {observed:7} '
              f'{"PASS" if correct else "FAIL"}')
        if case == 'xz_plane' and actual is not None:
            print('  values:', actual[2].tolist(), 'ErrorEvents:', len(errors))

print('\nFIELDDATA CONTROL: same logical values, two storage layouts')
for layout in ['flat', 'tuples']:
    path = OUT / f'fielddata_{layout}.vtkhdf'
    with h5py.File(path, 'w') as f:
        g = header(f, [0, 0, 0], [2, 2, 2])
        data = np.array([[10, 20], [30, 40]], dtype=np.int32)
        g['FieldData']['v'] = data.ravel() if layout == 'flat' else data
        steps = g.create_group('Steps')
        steps.attrs['NSteps'] = np.int32(2)
        steps['Values'] = np.array([0., 1.])
        offsets = [0, 2] if layout == 'flat' else [0, 1]
        steps.create_group('FieldDataOffsets')['v'] = np.array(offsets, dtype=np.int64)
        steps.create_group('FieldDataSizes')['v'] = np.array([[2, 1], [2, 1]], dtype=np.int32)
    for step in [0, 1]:
        errors, image, actual = read(path, 'FieldData', 'v', step)
        print(layout, 'step', step, 'components', actual[0],
              'tuples', actual[1], 'values', actual[2].tolist())
```

With VTK 9.7.1, the geometry section prints:

```text
xyz_zero       PointData 24 -> 24      PASS
xyz_zero       CellData   6 -> 6       PASS
xyz_positive   PointData 24 -> missing FAIL
xyz_positive   CellData   6 -> missing FAIL
xyz_negative   PointData 24 -> 24      PASS
xyz_negative   CellData   6 -> 6       PASS
xy_plane       PointData 12 -> 12      PASS
xy_plane       CellData   6 -> 6       PASS
yz_plane       PointData 12 -> missing FAIL
yz_plane       CellData   6 -> missing FAIL
xz_plane       PointData 12 -> 4       FAIL
  values: [0.0, 1.0, 2.0, 3.0] ErrorEvents: 0
xz_plane       CellData   6 -> 0       FAIL
  values: [] ErrorEvents: 0
x_line         PointData  4 -> 4       PASS
x_line         CellData   3 -> 3       PASS
y_line         PointData  4 -> missing FAIL
y_line         CellData   3 -> missing FAIL
z_line         PointData  4 -> missing FAIL
z_line         CellData   3 -> missing FAIL
```

The tested development build prints `PASS` for all 18 geometry checks.

## Separate finding: fVTKHDF temporal FieldData representation

This is an interoperability finding affecting the existing UG layout and its
ImageData counterpart. The controls do **not** support attributing it to a
general inability of VTK to read multicomponent temporal FieldData.

The reproducer compares two representations of one two-component tuple at each
of two time steps. Both have `FieldDataSizes/v = [[2,1],[2,1]]`:

| Representation | `FieldData/v` values and shape | `FieldDataOffsets/v` | VTK 9.7.1 result at steps 0 and 1 |
| --- | --- | --- | --- |
| Flattened scalar stream, matching fVTKHDF's former approach | `[10,20,30,40]`, shape `(4,)` | `[0,2]` | One component, one tuple: `[10]`, then `[30]` |
| Tuple-by-component storage | `[[10,20],[30,40]]`, shape `(2,2)` | `[0,1]` | Two components, one tuple: `[[10,20]]`, then `[[30,40]]` |

The tested development build produces exactly the same FieldData results as
9.7.1. The working control demonstrates that a component dimension in the
dataset and offsets measured along the tuple axis preserve both components.

This matches the native implementation: VTK 9.7.1's
[FieldData writer](https://github.com/Kitware/VTK/blob/v9.7.1/IO/HDF/vtkHDFWriter.cxx#L2340-L2357)
passes the component count to
[InitDynamicDataset](https://github.com/Kitware/VTK/blob/v9.7.1/IO/HDF/vtkHDFWriterImplementation.cxx#L1703-L1713),
which retains a component dimension for multicomponent arrays.
[AppendDataArrayOffset](https://github.com/Kitware/VTK/blob/v9.7.1/IO/HDF/vtkHDFWriter.cxx#L2640-L2690)
advances ordinary attribute offsets by tuple counts. The
[reader's NewFieldArray](https://github.com/Kitware/VTK/blob/v9.7.1/IO/HDF/vtkHDFUtilities.cxx#L1209-L1264)
uses tuple selections and the stored dataset dimensions; `FieldDataSizes`
does not turn a one-dimensional scalar stream into component tuples.

Matching ImageData output to the former UG layout established behavioral
parity, but did not establish format conformance. Both implementations now
write tuple-by-component datasets with two extendible dimensions and offsets
measured in tuples. The public API still allows component and tuple counts to
vary between steps; per-step sizes exclude padding in earlier tuple rows.
The upstream [FieldDataSizes proposal](https://discourse.vtk.org/t/improve-fielddata-in-vtkhdf-temporal-format/13669)
explicitly describes changing component and tuple counts and illustrates the
padded representation.

Strict regression tests replace the former multicomponent reader skip. They
cover all four FieldData numeric kinds, changing component counts, changing
tuple counts, noncontiguous input, initial omission, later reuse, empty writes,
and fields never written. These tests also read both UG and ImageData files.
The standalone flattened-layout control above is retained to document the
original failure; current fVTKHDF no longer writes that representation.
