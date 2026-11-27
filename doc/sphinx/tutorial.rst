An fVTKHDF Tutorial
===================

In this tutorial, we'll write several small Fortran programs step by step
to create VTKHDF files with different dataset types. After each exercise,
we'll compile and run the program, then open the resulting file in
ParaView to inspect what we wrote. The exercises cover a basic selection of
dataset types and writer features, and each includes a complete Fortran
program.
The :repo:`example directory <example/README.md>` shows additional features
in complete programs, and the :doc:`reference manual <reference>`
describes the available dataset types and writer interfaces in detail.

To build the programs, you’ll need CMake 3.28 or newer, an installed
MPI-enabled fVTKHDF library, and the compatible C, Fortran, MPI, and
parallel HDF5 environment used to build it. You’ll also need a version of
ParaView that supports the VTKHDF datasets in the exercises. See the
:doc:`installation` for how to build and install the library and
:ref:`using-fvtkhdf` for how to build your tutorial programs.

Temperature on a regular grid (ImageData)
-----------------------------------------

We’ll write a temperature field on a regular two-dimensional grid and view
it in ParaView. We’ll start with one MPI process to get the file-writing
sequence working for a single snapshot. Then we’ll make the hot spot move
and distribute the same grid across two processes. We’ll finish by writing
temperature at the grid points and assigning the shared points to one rank.

.. figure:: _static/tutorial/temperature-static.png
   :alt: Rectangular grid colored by temperature, with a hot spot above and left of center.
   :width: 65%

   The first exercise places a hot spot above and left of the center of the grid.

ImageData is a natural fit for this regular grid because its geometry
is defined by point bounds, an origin, and uniform spacing along each
axis. We’ll use the ``vtkhdf_image_file_type`` module and declare a
``vtkhdf_image_file`` object for the output.

.. code-block:: fortran

   use vtkhdf_image_file_type
   type(vtkhdf_image_file) :: file

Write a static temperature field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We'll start by writing a single snapshot of the temperature at the cell
centers using one MPI process. Before creating the file, we initialize MPI
and check that ``MPI_COMM_WORLD`` has the expected size.

.. code-block:: fortran

   use mpi

   ! Declarations omitted.

   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
   if (nproc /= 1) then
     call MPI_Finalize(ierr)
     error stop 'Run image_static with exactly 1 MPI process.'
   end if

Our grid covers a 1 by 0.75 rectangle in the z = 0 plane. We use
``nx = 24`` cells along x and ``ny = 18`` cells along y. The physical
origin is ``(0.0, 0.0, 0.0)``, with cell widths ``dx = 1.0/nx`` and
``dy = 0.75/ny``.

We set the inclusive logical point bounds to ``[0,0,0]`` at the lower
corner and ``[nx,ny,0]`` at the upper corner. That gives ``nx+1`` points
along x and ``ny+1`` along y. The lower and upper z bounds are both
0, making this a two-dimensional grid. We retain a singleton third
dimension in the Fortran cell temperature array, which has shape
``(nx,ny,1)``.

We set ``array_origin`` to the global logical index of the first array
element. A cell takes the logical index of its lower corner, the point
with the smaller index along both x and y.

.. figure:: _static/tutorial/cell-indexing.svg
   :alt: The lower-left cell stores temperature(1,1,1). Its lower corner has
         logical point index (0, 0, 0), and its physical center is at
         (dx/2, dy/2, 0.0).
   :width: 100%

   A close-up of the lower-left corner of the grid.

The half-cell offsets in the loop give the physical coordinates of
each cell center, where we evaluate a Gaussian centered at ``(0.35,0.45)``.

.. literalinclude:: ../../example/tutorial/image_static.F90
   :language: fortran
   :start-at: point_bounds(:,1)
   :end-before: #ifdef USE_MPI

We create the file in static mode with ``MPI_COMM_WORLD``. The call
returns ``stat``, which we check before writing. If creation fails, we
report ``errmsg``.

.. code-block:: fortran

   call file%create('image_static.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_STATIC)
   if (stat /= 0) then
     call MPI_Finalize(ierr)
     error stop errmsg
   end if

We write the mesh once and then write temperature as cell data.
``write_mesh`` takes a three-component spacing vector. We pass ``1.0``
for z as a placeholder because the only logical z index is 0, so the
z spacing has no effect on point coordinates. The writer places the
array using ``array_origin`` and its shape.

.. literalinclude:: ../../example/tutorial/image_static.F90
   :language: fortran
   :start-at: call file%write_mesh
   :end-at: call file%close

The program calls ``MPI_Finalize`` after closing the file. The
:doc:`ImageData <image>` and :doc:`structured writer <structured>`
reference pages describe the mesh and cell-data arguments in detail.

.. raw:: html

   <details>
   <summary>View complete source <code>image_static.F90</code></summary>

.. literalinclude:: ../../example/tutorial/image_static.F90
   :language: fortran

:download:`Download image_static.F90 <../../example/tutorial/image_static.F90>`

.. raw:: html

   </details>

Follow :ref:`using-fvtkhdf`, replacing ``my_app`` with ``image_static``
and using the downloaded
:download:`image_static.F90 <../../example/tutorial/image_static.F90>`
as the source. Include the ``USE_MPI`` target definition so the
preprocessor selects the MPI code. Then run the compiled program as shown
below, assuming the executable is in ``build/`` relative to the working
directory.

.. code-block:: bash

   mpiexec -n 1 ./build/image_static

The run produces ``image_static.vtkhdf`` in the current directory.
Open it in ParaView. The default **Slice** representation can appear
blank. Switch **Representation** to **Surface** or **Surface With Edges**
to show the grid. Color by ``temperature`` and look for the off-center
hot spot shown above.

Write a time-dependent temperature field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now we’ll make the hot spot move by writing a new temperature field at
each time step. We’ll keep the same mesh, array shape, and one-process
MPI setup.

We create ``image_temporal.vtkhdf`` with ``mode=VTKHDF_FIXED_MESH`` and
write the mesh once. Before the first time step, we register ``temperature``.
The real value ``mold=0.0`` specifies the type and kind of the scalar
field. We keep the returned handle in ``temperature_var`` so each step
can write to the registered field.

.. code-block:: fortran

   type(vtkhdf_cell_data_handle) :: temperature_var

.. literalinclude:: ../../example/tutorial/image_temporal.F90
   :language: fortran
   :start-at: temperature_var = file%register
   :end-at: call file%close

The Gaussian starts at the same position as in the static example and
moves along x as time goes from 0 to 1. We start and finalize each time
step around its field write. After the loop, we close the file and
finalize MPI.

.. raw:: html

   <details>
   <summary>View complete source <code>image_temporal.F90</code></summary>

.. literalinclude:: ../../example/tutorial/image_temporal.F90
   :language: fortran

:download:`Download image_temporal.F90 <../../example/tutorial/image_temporal.F90>`

.. raw:: html

   </details>

Build ``image_temporal`` from ``image_temporal.F90``, again
defining ``USE_MPI``.

.. code-block:: bash

   mpiexec -n 1 ./build/image_temporal

Open ``image_temporal.vtkhdf``, choose **Surface** or **Surface With Edges**,
and color by ``temperature``. Play the time series to see the hot spot
move to the right on the fixed grid.

.. figure:: _static/tutorial/temperature-temporal.gif
   :alt: The hot spot moves from left to right across the fixed rectangular grid.
   :width: 65%

   Watch the hot spot move from x = 0.35 to x = 0.70 while the cell edges stay fixed.

.. note::

   To build the first two programs against a serial installation, follow
   the serial instructions in :ref:`using-fvtkhdf` and use
   ``build-serial`` as the build directory. Run the programs directly
   from the project directory.

   .. code-block:: bash

      ./build-serial/image_static
      ./build-serial/image_temporal

   The two-process programs require an MPI-enabled installation.

Distribute the same field across two MPI processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next we’ll write the same time series to one file using two MPI
processes. The writer calls are collective over ``MPI_COMM_WORLD``, so
both ranks make the same calls in the same order. They supply the same
global mesh description, and each rank supplies its own local field values.

We split the grid into two equal slabs along x and use each process’s
rank in ``MPI_COMM_WORLD`` to decide which slab it writes.
Each rank owns 12 by 18 cells and allocates space for its own slab.

.. figure:: _static/tutorial/cell-slabs.svg
   :alt: Rank 0 owns x cell indices 0 through 11 on the left. Rank 1 owns
         x cell indices 12 through 23 on the right. The slabs meet at x = 0.5.
   :width: 56.25%

Both slabs span the full y range. Each rank sets ``array_origin``
from its rank number and the slab width.

.. literalinclude:: ../../example/tutorial/image_parallel.F90
   :language: fortran
   :start-at: array_origin =
   :end-at: allocate(temperature

After writing the mesh, we fill an integer array of shape
``(local_nx,ny,1)`` with the local MPI rank. Both ranks
write their slab to the static cell array ``mpi_rank`` using the same
``array_origin`` as for temperature. The ownership stays fixed for every
time step, so we write it once before the time loop.

.. literalinclude:: ../../example/tutorial/image_parallel.F90
   :language: fortran
   :start-at: cell_rank = rank
   :end-at: call file%write_cell_data('mpi_rank'

We use the slab’s global starting index when calculating x so both ranks
evaluate the Gaussian at the same physical positions as the one-process
program. For example, rank 1’s ``temperature(1,j,1)`` is centered at
``x = (12 + 0.5)*dx``.

.. literalinclude:: ../../example/tutorial/image_parallel.F90
   :language: fortran
   :start-at: do step = 0, 10
   :end-before: call file%close

Each rank supplies exactly its owned cells, so the array shape and
``array_origin`` are enough to place both slabs without explicit
``write_bounds``. After the last time step, both ranks close the file
and finalize MPI. See the :doc:`common writer conventions <common>`
for collective call requirements.

.. raw:: html

   <details>
   <summary>View complete source <code>image_parallel.F90</code></summary>

.. literalinclude:: ../../example/tutorial/image_parallel.F90
   :language: fortran

:download:`Download image_parallel.F90 <../../example/tutorial/image_parallel.F90>`

.. raw:: html

   </details>

Build ``image_parallel`` from ``image_parallel.F90`` and
launch it with two processes.

.. code-block:: bash

   mpiexec -n 2 ./build/image_parallel

Open ``image_parallel.vtkhdf`` in ParaView, set **Representation** to
**Surface** or **Surface With Edges**, and color by ``mpi_rank``.
Look for rank 0 on the left and rank 1 on the right, meeting at
``x = 0.5``. Switch to ``temperature`` and play the time series.
Check that the grid still spans the 1 by 0.75 rectangle with 24 by 18 cells
and that the hot spot follows the same path across the slab boundary.
Both processes have now written their parts of the same temperature field
into one file.

Write point temperature on two MPI processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This time we’ll write temperature at the grid points, including points
shared by the two slabs. The mesh, the slabs, and the path of the hot
spot stay the same. Each slab is 12 cells wide, so its temperature array
has 13 point columns and 19 point rows, with a singleton z dimension.
The calculation uses point coordinates, so there is no half-cell offset.

The two local arrays both include the point column at logical index 12
(``x = 0.5``). We assign shared points to the lower rank so each global
point is written exactly once.

.. list-table::
   :header-rows: 1

   * - Rank
     - Point columns in local array
     - Point columns written
   * - 0
     - 0–12
     - 0–12
   * - 1
     - 12–24
     - 13–24

We keep ``array_origin`` at the first point stored in each local array.
We start with ``write_bounds`` covering the whole array, then advance
rank 1’s lower x write bound by one to skip its first column. Both
``array_origin`` and ``write_bounds`` use global logical indices.

.. literalinclude:: ../../example/tutorial/image_parallel_points.F90
   :language: fortran
   :start-at: array_origin =
   :end-at: if (rank == 1) write_bounds

To see who wrote each point, we fill a point array with the local MPI
rank and write it once as ``mpi_rank``. We use the same ``write_bounds``
for ``mpi_rank`` and ``temperature``.

.. literalinclude:: ../../example/tutorial/image_parallel_points.F90
   :language: fortran
   :start-at: point_rank = rank
   :end-at: call file%write_point_data('mpi_rank'

We register ``temperature`` as temporal point data, so this version
declares ``temperature_var`` as a point-data handle.

.. code-block:: fortran

   type(vtkhdf_point_data_handle) :: temperature_var

Both ranks follow the same time loop and call ``write_point_data`` with
their own bounds at every step.

.. literalinclude:: ../../example/tutorial/image_parallel_points.F90
   :language: fortran
   :start-at: temperature_var = file%register
   :end-before: call file%close

.. raw:: html

   <details>
   <summary>View complete source <code>image_parallel_points.F90</code></summary>

.. literalinclude:: ../../example/tutorial/image_parallel_points.F90
   :language: fortran

:download:`Download image_parallel_points.F90 <../../example/tutorial/image_parallel_points.F90>`

.. raw:: html

   </details>

Build ``image_parallel_points`` from
``image_parallel_points.F90`` and run it with two processes too.

.. code-block:: bash

   mpiexec -n 2 ./build/image_parallel_points

Open ``image_parallel_points.vtkhdf`` in ParaView. Choose **Points** and
color by ``mpi_rank`` to check that the center column has the same color
as the points on its left. Then switch to **Surface With Edges**, color
by ``temperature``, and watch the hot spot move through the shared
column. Surface coloring interpolates between point values within each
cell.

We’ve written the complete point field at every time step, with rank 0
writing the shared column once per step.

Temperature on a nonuniform grid (RectilinearGrid)
--------------------------------------------------

We’ll start with the two-rank, cell-centered ImageData program and vary
the spacing along y to put more rows near the hot spot’s path.
RectilinearGrid takes a vector of point coordinates for each axis, so
neighboring coordinates can be different distances apart. To inspect
the files in this chapter, use ParaView 6.2 or later.

.. code-block:: fortran

   use vtkhdf_rectilinear_file_type
   type(vtkhdf_rectilinear_file) :: file

Write temperature on a fixed nonuniform grid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We’ll keep the 1 by 0.75 domain, the 24 by 18 cells, and the same moving
Gaussian. The mesh stays fixed throughout the run. We space the points
uniformly along x and concentrate them near y = 0.45, where the hot spot
travels.

To set the y coordinates, we start with evenly spaced values of ``s``
from −1 to 1.
The factor ``0.35 + 0.65*s**2`` draws the interior points toward the
focus while keeping the endpoints fixed. We scale each side by its
distance from the focus to the corresponding end of the domain.

.. literalinclude:: ../../example/tutorial/rectilinear_parallel.F90
   :language: fortran
   :start-at: function clustered_axis
   :end-at: end function

Both ranks compute the same vectors because both need the coordinates
to evaluate temperature. Both call ``write_mesh`` with the same global
point bounds. The writer takes the coordinate vectors from rank 0.
The singleton z axis supplies one coordinate, ``0.0``.

.. literalinclude:: ../../example/tutorial/rectilinear_parallel.F90
   :language: fortran
   :start-at: x_points =
   :end-at: call file%write_mesh

The logical grid and cell slabs work as before. The coordinate arrays
start at Fortran subscript 0, so we can use global logical indices
directly to look up the points. To find a cell center, we average its
two bounding coordinates along each axis. We use ``array_origin`` to
find the global x index for each local cell, and both slabs start at
logical index 0 along y.

.. literalinclude:: ../../example/tutorial/rectilinear_parallel.F90
   :language: fortran
   :start-at: do j = 1, ny
   :end-at: call file%write_cell_data(temperature_var

The collective time loop and the writes for ``mpi_rank`` work as in the
ImageData version. See the :doc:`RectilinearGrid reference <rectilinear>`
for the coordinate-vector arguments.

.. raw:: html

   <details>
   <summary>View complete source <code>rectilinear_parallel.F90</code></summary>

.. literalinclude:: ../../example/tutorial/rectilinear_parallel.F90
   :language: fortran

:download:`Download rectilinear_parallel.F90 <../../example/tutorial/rectilinear_parallel.F90>`

.. raw:: html

   </details>

Build ``rectilinear_parallel`` from ``rectilinear_parallel.F90``
and launch it with two processes.

.. code-block:: bash

   mpiexec -n 2 ./build/rectilinear_parallel

Open ``rectilinear_parallel.vtkhdf`` in ParaView and color by
``temperature`` with **Surface With Edges**. Play the time series and
look for the narrow rows near y = 0.45.
The temperature changes from step to step while the cell edges stay put.

.. figure:: _static/tutorial/temperature-rectilinear.gif
   :alt: The hot spot moves along a band of closely spaced rows on a fixed
         grid with evenly spaced columns.
   :width: 65%

   The hot spot moves along the band of closely spaced rows.

Move the mesh with the hot spot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We’ll now cluster the x coordinates as well, moving their concentration
with the hot spot. We keep the y coordinates from the fixed-mesh
exercise and use the same ``clustered_axis`` function for x, with the
current Gaussian center as its focus. The ends of the domain stay fixed.
Each rank still owns 12 columns, and the physical boundary between the
two slabs now follows the hot spot.

We create the file with ``mode=VTKHDF_MOVING_MESH``. We pass only the
point bounds to ``write_mesh`` to define the logical grid once. Each
logical cell keeps its owner, so we write the static ``mpi_rank`` array
before supplying the coordinates. We register ``temperature`` as temporal
cell data before the first time step.

.. literalinclude:: ../../example/tutorial/rectilinear_moving.F90
   :language: fortran
   :start-at: y_points = clustered_axis
   :end-at: temperature_var = file%register_temporal_cell_data

At each step, both ranks calculate the new x coordinates and call
``write_mesh_geometry`` before writing temperature. This call takes
all three coordinate vectors, including the same y and z coordinates
each time.

.. literalinclude:: ../../example/tutorial/rectilinear_moving.F90
   :language: fortran
   :start-at: do step = 0, 10
   :end-at: call file%write_mesh_geometry

The temperature loop still averages neighboring point coordinates to
find each cell center. Because it uses the current x vector, the centers
move with the mesh. We write the resulting temperatures and finalize
the step as before.

.. raw:: html

   <details>
   <summary>View complete source <code>rectilinear_moving.F90</code></summary>

.. literalinclude:: ../../example/tutorial/rectilinear_moving.F90
   :language: fortran

:download:`Download rectilinear_moving.F90 <../../example/tutorial/rectilinear_moving.F90>`

.. raw:: html

   </details>

Build ``rectilinear_moving`` from ``rectilinear_moving.F90`` and run
it with two processes.

.. code-block:: bash

   mpiexec -n 2 ./build/rectilinear_moving

Open ``rectilinear_moving.vtkhdf`` in ParaView and color by
``temperature`` with **Surface With Edges**. As you play the time series,
watch the narrow columns follow the hot spot while the row boundaries
stay in place.

.. figure:: _static/tutorial/temperature-rectilinear-moving.gif
   :alt: Closely spaced vertical grid lines follow the hot spot from left
         to right while the horizontal grid lines stay fixed.
   :width: 65%

   The mesh concentrates cells around the hot spot at each time step.

Temperature on a curvilinear grid (StructuredGrid)
--------------------------------------------------

We’ll take the two-rank point-data ImageData program and bend its logical
rectangle into a 60° annular sector. The cells still connect to their logical
neighbors, and each point now has its own physical x, y, and z coordinates.
We’ll write one snapshot of the curved grid, then rotate the grid through
a stationary hot spot. To inspect the files in this chapter, use ParaView
6.2 or later.

.. code-block:: fortran

   use vtkhdf_structured_file_type
   type(vtkhdf_structured_file) :: file

Write a point field on an annular sector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We’ll write a static temperature field on an annular sector in the z = 0
plane. We keep the ImageData grid’s 24 by 18 cells. The first logical
index steps clockwise from 100° to 40°, and the second steps outward
from radius 0.65 to 1.25. The Gaussian is centered at (0.0, 0.95),
halfway across the radial span and near the 100° end of the sector.

We split the grid along the first logical axis as before. Each rank’s
array holds 13 angular points and 19 radial points. Rank 0 owns shared
column 12, which now forms a radial line at 70°. The ``points`` array has
a leading dimension of length 3 for the x, y, and z coordinates of each
point.

.. literalinclude:: ../../example/tutorial/structured_parallel.F90
   :language: fortran
   :start-at: real :: points
   :end-at: real :: temperature

We use the global first index to get the angle and the second index to
get the radius. Then we convert from polar to Cartesian coordinates and
evaluate the temperature at the same point.

.. literalinclude:: ../../example/tutorial/structured_parallel.F90
   :language: fortran
   :start-at: do j = 1, ny+1
   :end-before: call file%create

We create the file with ``mode=VTKHDF_STATIC``. Both ranks pass the same
global point bounds to ``write_mesh``, then call ``write_mesh_geometry``
with their local coordinates. We use the same ``array_origin`` and
``write_bounds`` for geometry, ``mpi_rank``, and ``temperature``, so
rank 0 writes all three arrays along the shared radial line. See the
:doc:`StructuredGrid reference <structured_grid>` for the coordinate
array and bounds arguments.

.. literalinclude:: ../../example/tutorial/structured_parallel.F90
   :language: fortran
   :start-at: call file%write_mesh(point_bounds)
   :end-at: call file%write_point_data('temperature'

.. raw:: html

   <details>
   <summary>View complete source <code>structured_parallel.F90</code></summary>

.. literalinclude:: ../../example/tutorial/structured_parallel.F90
   :language: fortran

:download:`Download structured_parallel.F90 <../../example/tutorial/structured_parallel.F90>`

.. raw:: html

   </details>

Build ``structured_parallel`` from ``structured_parallel.F90`` and run
it with two processes.

.. code-block:: bash

   mpiexec -n 2 ./build/structured_parallel

Open ``structured_parallel.vtkhdf`` in ParaView, choose **Surface With
Edges**, and color by point ``temperature``. Look for the hot spot halfway
between the inner and outer edges, near the left end of the sector.
To check ownership, switch to **Points** and color by ``mpi_rank``.
The shared radial line belongs to rank 0.

.. figure:: _static/tutorial/temperature-structured.png
   :alt: A 60° annular sector with a hot spot near its left end.
   :width: 65%

   The hot spot sits near the 100° end of the sector.

Rotate the grid through a fixed hot spot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now we’ll rotate the sector 40° counterclockwise about the origin while
keeping the Gaussian centered at (0.0, 0.95) in physical space. The
sector finishes at 80° to 140°, so the hot spot is near the opposite
end of the grid.

We create the file with ``mode=VTKHDF_MOVING_MESH`` and define the logical
grid once. Each logical point keeps its owner as the grid moves, so we can
write the static ``mpi_rank`` array during setup, before writing any
coordinates. The grid points move through the Gaussian, so we register
``temperature`` as temporal point data.

.. literalinclude:: ../../example/tutorial/structured_moving.F90
   :language: fortran
   :start-at: call file%write_mesh(point_bounds)
   :end-at: temperature_var = file%register_temporal_point_data

At each step, we add the current rotation to each point’s angle and
evaluate the same Gaussian at the new coordinates. We write the geometry
first, then the temperatures.

.. literalinclude:: ../../example/tutorial/structured_moving.F90
   :language: fortran
   :start-at: do step = 0, 10
   :end-before: call file%close

.. raw:: html

   <details>
   <summary>View complete source <code>structured_moving.F90</code></summary>

.. literalinclude:: ../../example/tutorial/structured_moving.F90
   :language: fortran

:download:`Download structured_moving.F90 <../../example/tutorial/structured_moving.F90>`

.. raw:: html

   </details>

Build ``structured_moving`` from ``structured_moving.F90`` and run it
with two processes.

.. code-block:: bash

   mpiexec -n 2 ./build/structured_moving

Open ``structured_moving.vtkhdf`` in ParaView, choose **Surface With
Edges**, and color by point ``temperature``. Play the time series and
watch the hot spot stay in place as the grid rotates through it. The
hot spot starts near one end of the sector and finishes near the other.

.. figure:: _static/tutorial/temperature-structured-moving.gif
   :alt: The annular sector rotates counterclockwise while the hot spot
         stays at the same x and y coordinates.
   :width: 65%

   The sector rotates through a Gaussian whose center stays at (0.0, 0.95).
