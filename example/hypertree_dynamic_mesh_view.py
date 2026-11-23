"""Save a ParaView state showing the full moving-circle grid and interface overlay.

Run with pvpython --force-offscreen-rendering to write a .pvsm file and exit
without opening a window. Load that file in ParaView through File > Load State.
"""
import argparse
from pathlib import Path


def create_view(filename):
    from paraview.simple import (
        AssignViewToLayout, ColorBy, CreateLayout, CreateView, GetAnimationScene,
        GetColorTransferFunction, GetScalarBar,
        HyperTreeGridGeometryFilter, HyperTreeGridToUnstructuredGrid,
        OpenDataFile, Show, _DisableFirstRenderCameraReset,
    )

    _DisableFirstRenderCameraReset()
    reader = OpenDataFile(str(Path(filename).resolve()))
    reader.UpdatePipeline()
    mesh = HyperTreeGridToUnstructuredGrid(registrationName='Full adaptive grid', Input=reader)
    interface = HyperTreeGridGeometryFilter(registrationName='Material interface', Input=reader)
    interface.FillMaterial = 0

    view = CreateView('RenderView')
    layout = CreateLayout(name='Moving circle')
    AssignViewToLayout(view=view, layout=layout, hint=0)
    view.ViewSize = [1200, 800]
    view.Background = [1., 1., 1.]
    view.UseColorPaletteForBackground = 0
    mesh_display = Show(mesh, view)
    mesh_display.Representation = 'Surface With Edges'
    mesh_display.EdgeColor = [.25, .25, .25]
    ColorBy(mesh_display, ('CELLS', 'level'))
    colors = GetColorTransferFunction('level')
    colors.ApplyPreset('Viridis', True)
    colors.RescaleTransferFunction(0, 4)
    mesh_display.SetScalarBarVisibility(view, True)
    legend = GetScalarBar(colors, view)
    legend.Title = 'Refinement level'
    legend.ComponentTitle = ''
    legend.TitleColor = [.1, .1, .1]
    legend.LabelColor = [.1, .1, .1]
    legend.WindowLocation = 'Upper Right Corner'

    interface_display = Show(interface, view)
    ColorBy(interface_display, None)
    interface_display.DiffuseColor = [.9, .1, .1]
    interface_display.LineWidth = 3.
    # Lift the line display slightly above the 2D surface to avoid z-fighting.
    interface_display.Translation = [0., 0., .001]

    view.CameraPosition = [3., 2., 20.]
    view.CameraFocalPoint = [3., 2., 0.]
    view.CameraViewUp = [0., 1., 0.]
    view.CameraParallelProjection = 1
    view.CameraParallelScale = 2.6
    GetAnimationScene().UpdateAnimationUsingDataTimeSteps()
    return reader, mesh, interface, view


def main():
    from paraview.simple import SaveScreenshot, SaveState

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('filename', nargs='?', default='hypertree_dynamic_mesh.vtkhdf')
    parser.add_argument('--state', help='output .pvsm path; defaults to the data file with .pvsm suffix')
    parser.add_argument('--screenshot', help='optional output image path')
    args = parser.parse_args()
    _, _, _, view = create_view(args.filename)
    state = args.state or str(Path(args.filename).with_suffix('.pvsm'))
    SaveState(state)
    if args.screenshot:
        SaveScreenshot(args.screenshot, view)
    print(f'Load this state in ParaView: {Path(state).resolve()}')


if __name__ == '__main__':
    main()
