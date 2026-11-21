"""Check the moving-circle example through VTK's full-grid and interface geometry filters."""
import sys

try:
    import numpy as np
    from vtkmodules.vtkCommonCore import vtkVersion
    from vtkmodules.vtkCommonDataModel import vtkHyperTreeGrid, vtkHyperTreeGridNonOrientedGeometryCursor
    from vtkmodules.vtkFiltersHyperTree import vtkHyperTreeGridGeometry, vtkHyperTreeGridToUnstructuredGrid
    from vtkmodules.vtkIOHDF import vtkHDFReader
    from vtkmodules.util.numpy_support import vtk_to_numpy
except ImportError as exc:
    print(exc)
    raise SystemExit(77)

if (vtkVersion.GetVTKMajorVersion(), vtkVersion.GetVTKMinorVersion()) < (9, 7):
    raise SystemExit(77)

reader = vtkHDFReader()
reader.SetFileName(sys.argv[1])
errors = []
reader.AddObserver('ErrorEvent', lambda *args: errors.append(args))
reader.UpdateInformation()
assert reader.GetNumberOfSteps() == 9
radius = 1.2
cell_counts = []


def check_tree(grid, cursor, center):
    """Check values in every cell, including refined ancestors."""
    index = cursor.GetGlobalNodeIndex()
    cell_data = grid.GetCellData()
    normal = np.array(cell_data.GetArray('vtkInterfaceNormals').GetTuple(index))
    intercept = np.array(cell_data.GetArray('vtkInterfaceIntercepts').GetTuple(index))
    bounds = [0.]*6
    cursor.GetBounds(bounds)
    origin = np.array([bounds[0], bounds[2]])
    size = np.array([bounds[1]-bounds[0], bounds[3]-bounds[2]])
    corners = origin + size*np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    phi = np.sum((corners-center)**2, axis=1)-radius**2
    mixed = phi.min() < 0 < phi.max()
    assert intercept[2] == (1 if mixed else 2)
    assert not grid.GetMask().GetTuple1(index)  # both sides remain in the full grid
    assert intercept[0] == 0 and normal[2] == 0
    distance = np.linalg.norm(origin+size/2-center)-radius
    np.testing.assert_allclose(cell_data.GetArray('distance').GetTuple1(index), distance, atol=2e-6)
    assert cell_data.GetArray('level').GetTuple1(index) == cursor.GetLevel()
    if mixed:
        np.testing.assert_allclose(np.linalg.norm(normal), 1, atol=2e-6)
        assert np.dot(normal[:2], origin+size/2-center) > 0  # outward
        # Independently locate the interpolated circle crossings on cell edges.
        for a in range(4):
            b = (a+1) % 4
            if (phi[a] <= 0) != (phi[b] <= 0):
                crossing = corners[a]+phi[a]/(phi[a]-phi[b])*(corners[b]-corners[a])
                np.testing.assert_allclose(np.dot(normal[:2], crossing)+intercept[1], 0, atol=3e-6)
        if cursor.IsLeaf():
            assert cursor.GetLevel() == 4
    else:
        assert np.all(normal == 0)
    count = 1
    if not cursor.IsLeaf():
        for child in range(4):
            cursor.ToChild(child)
            count += check_tree(grid, cursor, center)
            cursor.ToParent()
    return count


for step in range(9):
    reader.SetStep(step)
    reader.Update()
    assert not errors, errors
    center = np.array([.5+5*step/8, 1+2*step/8])
    output = reader.GetOutputDataObject(0)
    pieces = ([output] if isinstance(output, vtkHyperTreeGrid) else
              [output.GetPartitionAsDataObject(i) for i in range(output.GetNumberOfPartitions())])
    area = 0.
    boundary_points = 0
    cells = 0
    for rank, grid in enumerate(pieces):
        assert grid.GetHasInterface()
        assert grid.GetInterfaceNormalsName() == 'vtkInterfaceNormals'
        assert grid.GetInterfaceInterceptsName() == 'vtkInterfaceIntercepts'
        cursor = vtkHyperTreeGridNonOrientedGeometryCursor()
        count = 0
        for tree in range(rank, 24, len(pieces)):
            grid.InitializeNonOrientedGeometryCursor(cursor, tree, False)
            assert cursor.HasTree()
            count += check_tree(grid, cursor, center)
        assert count == grid.GetNumberOfCells()
        cells += count

        # This conversion preserves full leaf cells on BOTH sides of the interface.
        full_grid = vtkHyperTreeGridToUnstructuredGrid()
        full_grid.AddObserver('ErrorEvent', lambda *args: errors.append(args))
        full_grid.SetInputData(grid)
        full_grid.Update()
        mesh = full_grid.GetOutput()
        assert mesh.GetNumberOfCells() == grid.GetNumberOfLeaves()
        for i in range(mesh.GetNumberOfCells()):
            polygon = mesh.GetCell(i)
            assert polygon.GetNumberOfPoints() == 4
            bounds = polygon.GetBounds()
            area += (bounds[1]-bounds[0])*(bounds[3]-bounds[2])

        # Interface-only output traces the circle on top, without internal edges.
        geometry = vtkHyperTreeGridGeometry()
        geometry.AddObserver('ErrorEvent', lambda *args: errors.append(args))
        geometry.SetInputData(grid)
        geometry.FillMaterialOff()
        geometry.Update()
        boundary = geometry.GetOutput()
        if boundary.GetNumberOfPoints():
            xy = vtk_to_numpy(boundary.GetPoints().GetData())[:, :2]
            np.testing.assert_allclose(np.linalg.norm(xy-center, axis=1), radius, atol=.002)
            assert boundary.GetNumberOfLines() > 0
            boundary_points += len(xy)
    assert not errors, errors
    assert boundary_points > 0
    np.testing.assert_allclose(area, 6*4, atol=1e-6)
    cell_counts.append(cells)
    print(f'step {step}: {cells} cells, full-grid area {area:g}, circle interface verified')
assert len(set(cell_counts)) > 1  # refinement follows the moving circle
print('Moving-circle material-interface geometry checks passed')
