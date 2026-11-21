"""Independent HDF5 layout and VTK component readback for paired interfaces."""
import sys
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError as exc:
    print(exc)
    raise SystemExit(77)

operation, directory = sys.argv[1:]
directory = Path(directory)
if operation == 'reader':
    try:
        from vtkmodules.vtkIOHDF import vtkHDFReader
        from vtkmodules.vtkCommonCore import vtkVersion
        from vtkmodules.vtkCommonDataModel import vtkHyperTreeGrid, vtkHyperTreeGridNonOrientedGeometryCursor
        from vtkmodules.util.numpy_support import vtk_to_numpy
    except ImportError as exc:
        print(exc)
        raise SystemExit(77)
    if (vtkVersion.GetVTKMajorVersion(), vtkVersion.GetVTKMinorVersion()) < (9, 7):
        raise SystemExit(77)


def ownership(nproc, phase):
    return [[t for t in range(1, 4) if (t + phase) % nproc == rank] for rank in range(nproc)]


def payload(ids, phase, step, intercepts=False):
    values = np.array([10*t+j+100*phase for t in ids
                       for j in range(1, (3 if (t+phase) % 2 == 0 else 1)+1)])
    if intercepts:
        return -2*values[:, None]+10*step-np.arange(1, 4)[None, :]*.25
    return values[:, None]+10*step+np.arange(1, 4)[None, :]*.125


def pieces(output):
    if isinstance(output, vtkHyperTreeGrid):
        return [output]
    return [output.GetPartitionAsDataObject(i) for i in range(output.GetNumberOfPartitions())]


def indices(grid, ids, phase):
    cursor = vtkHyperTreeGridNonOrientedGeometryCursor()
    result = []
    for t in ids:
        grid.InitializeNonOrientedGeometryCursor(cursor, t-1, False)
        assert cursor.HasTree()
        result.append(cursor.GetGlobalNodeIndex())
        if (t+phase) % 2 == 0:
            assert not cursor.IsLeaf()
            for child in range(2):
                cursor.ToChild(child)
                assert cursor.IsLeaf()
                result.append(cursor.GetGlobalNodeIndex())
                cursor.ToParent()
        else:
            assert cursor.IsLeaf()
    assert len(result) == grid.GetNumberOfCells()
    return np.array(result, dtype=int)


def check(path, kind, mode, temporal, reset=False):
    names = ('vtkInterfaceNormals', 'vtkInterfaceIntercepts') if reset else ('vtkInterfaceNormals_1', 'vtkInterfaceIntercepts_2')
    dtype = np.dtype('f4' if kind == 'real32' else 'f8')
    with h5py.File(path) as file:
        root = file['VTKHDF']
        nproc = int(root['Steps/NumberOfParts'][0]) if mode else len(root['NumberOfTrees'])
        steps = 4 if temporal else (2 if mode else 1)
        all_names = set(names) if reset else {
            'vtkInterfaceNormals', 'vtkInterfaceNormals_1', 'vtkInterfaceNormals_2',
            'vtkInterfaceIntercepts', 'vtkInterfaceIntercepts_1', 'vtkInterfaceIntercepts_2', 'vtkInterfaceIntercepts_3'}
        assert set(root['CellData']) == all_names, path
        for component, name in enumerate(names):
            base = ('vtkInterfaceNormals', 'vtkInterfaceIntercepts')[component]
            assert root.attrs[base.removeprefix('vtk')+'Name'] == name.encode(), path
            ds = root['CellData'][name]
            assert ds.dtype == dtype, (path, name, ds.dtype)
            assert ds.attrs['Name'] == base.encode()
            assert set(ds.attrs) == {'Name'}, (path, name, dict(ds.attrs))
            phases = [0, 1 if mode == 3 else 0] if temporal else [0]
            updates = [0, 2] if temporal else [0]
            want = np.concatenate([payload(sum(ownership(nproc, phase), []), phase, step, component == 1)
                                   for phase, step in zip(phases, updates)])
            np.testing.assert_array_equal(ds[:], want, err_msg=f'{path}:{name}')
            if temporal:
                np.testing.assert_array_equal(root['Steps/CellDataOffsets'][name], [0, 0, 5, 5])
                assert ds.chunks == (max(1, 64//(3*dtype.itemsize)), 3)
                assert root['Steps/CellDataOffsets'][name].chunks == (512,)
            else:
                assert ds.chunks is None
                if mode:
                    assert name not in root['Steps/CellDataOffsets']
        if mode:
            np.testing.assert_array_equal(root['Steps/Values'], np.arange(steps))
        # No active Normals or Vectors role is assigned by this API.
        assert not root['CellData'].attrs

    if operation != 'reader':
        return
    reader = vtkHDFReader()
    errors = []
    reader.AddObserver('ErrorEvent', lambda *args: errors.append(args))
    reader.SetFileName(str(path))
    reader.UpdateInformation()
    assert reader.GetNumberOfSteps() == steps
    # Read forwards and backwards to catch stale arrays and name changes.
    for step in list(range(steps)) + list(reversed(range(steps))):
        reader.SetStep(step)
        reader.Update()
        output = pieces(reader.GetOutputDataObject(0))
        assert len(output) == nproc
        phase = 1 if temporal and mode == 3 and step >= 2 else 0
        update = (2 if step >= 2 else 0) if temporal else 0
        for grid, ids in zip(output, ownership(nproc, phase)):
            assert grid.GetHasInterface()
            assert grid.GetInterfaceNormalsName() == names[0]
            assert grid.GetInterfaceInterceptsName() == names[1]
            order = indices(grid, ids, phase)
            for component, name in enumerate(names):
                array = grid.GetCellData().GetArray(name)
                assert array is not None and array.GetNumberOfComponents() == 3
                actual = vtk_to_numpy(array)
                assert actual.dtype == dtype
                np.testing.assert_array_equal(actual[order], payload(ids, phase, update, component == 1),
                                              err_msg=f'{path}:{name} step={step}')
        assert not errors, (path, errors)


for kind in ('real32', 'real64'):
    for temporal, modes in ((False, range(3)), (True, range(1, 4))):
        for mode in modes:
            path = directory / f'hypertree_interface_{kind}_{mode}_{"T" if temporal else "F"}.vtkhdf'
            check(path, kind, mode, temporal)
            check(Path(str(path)+'.reset'), kind, 0, False, reset=True)
print(f'HyperTreeGrid interface {operation} checks passed')
