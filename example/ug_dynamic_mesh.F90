program ug_dynamic_mesh_example

  use,intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_ug_file_type
  use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none

  real, allocatable :: x(:,:), points(:,:), pressure(:), velocity(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  character(:), allocatable :: errmsg
  integer :: stat, nproc, rank, ncells, npoints, step, j
  real :: dt, t

  type(vtkhdf_ug_file) :: file
  type(vtkhdf_cell_data_handle)  :: prs_var
  type(vtkhdf_point_data_handle) :: vel_var
  type(vtkhdf_field_data_handle) :: cpu_var

#ifdef USE_MPI
  call MPI_Init(stat)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, stat)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)
#else
  nproc = 1
  rank = 0
#endif

#ifdef USE_MPI
  call file%create('example.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=UG_DYNAMIC_MESH)
  if (stat /= 0) then
    call MPI_Finalize(stat)
    error stop errmsg
  end if
#else
  call file%create('example.vtkhdf', stat, errmsg, mode=UG_DYNAMIC_MESH)
  if (stat /= 0) error stop errmsg
#endif

  ! Initial 2D mesh. The local mesh on each rank is a 2-triangle discretization
  ! of a unit square, suitably x-shifted so they only overlap on a boundary.
  ! The mesh topologies are identical; point coordinates differ. 3D coordinates
  ! are required; we set z=0.

  points = reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0], shape=[3,4])
  points(1,:) = points(1,:) + rank
  cnode = [1,2,3,  1,3,4]
  xcnode = [1,4,7]
  types = spread(VTK_TRIANGLE, dim=1, ncopies=2)

  ! local mesh sizes
  ncells = size(types)
  npoints = size(points,dim=2)

  ! In UG_MOVING_MESH mode all datasets are necessarily temporal and need
  ! to be registered. For this example we have cell-centered pressure and
  ! point-centered velocity, and a scalar "cpu time" field data.

  associate (scalar_mold => 0.0, vector_mold => [0.0, 0.0, 0.0])
    prs_var = file%register_temporal_cell_data('pressure', scalar_mold)
    vel_var = file%register_temporal_point_data('velocity', vector_mold)
    cpu_var = file%register_temporal_field_data('cpu time', scalar_mold)
  end associate

  ! Start time stepping.
  dt = 4*atan(1.0)/10
  x = points ! initial point coordinates as reference
  do step = 0, 5
    t = step * dt
    call file%start_time_step(t)

    ! Field data can be written at any time during the step.
    call file%write_field_data(cpu_var, 5*t)

    ! Write mesh before mesh-centered data.
    points = (1 + 0.1*sin(t)**2) * x
    if (step == 0) then ! topology and geometry
      call file%write_mesh(points, cnode, xcnode, types)
    else ! new geometry; topology unchanged.
      call file%write_mesh_geometry(points)
    end if

    ! Write cell-centered pressure on the first step;
    ! subsequent steps reuse the most recently written value (this one).
    if (step == 0) then
      pressure = [(rank + j*0.2, j=0, ncells-1)]
      call file%write_cell_data(prs_var, pressure)
    end if

    ! Write point-centered velocity for this step.
    velocity = 0.1*sin(2*t) * x
    call file%write_point_data(vel_var, velocity)

    ! Return file to a complete, readable state until start of next step.
    call file%finalize_time_step
    call file%flush
  end do

  ! New mesh. Each unit square is subdivided into 4 triangles by adding
  ! a new point at the center.

  x = reshape([x, [rank+0.5, 0.5, 0.0]], shape=[3,5])
  cnode = [1,2,5,  2,3,5,  3,4,5,  4,1,5]
  xcnode = [1,4,7,10,13]
  types = spread(VTK_TRIANGLE, dim=1, ncopies=4)

  ! local mesh sizes
  ncells = size(types)
  npoints = size(x,dim=2)

  ! Continue time stepping.
  do step = 6, 10
    t = step * dt
    call file%start_time_step(t)

    ! Field data can be written at any time during the step.
    call file%write_field_data(cpu_var, 5*t)

    ! Write the new mesh before mesh-centered data.
    points = (1 + 0.1*sin(t)**2) * x
    if (step == 6) then ! topology and geometry
      call file%write_mesh(points, cnode, xcnode, types)
    else ! new geometry; topology unchanged.
      call file%write_mesh_geometry(points)
    end if

    ! Write cell-centered pressure on the new mesh;
    ! subsequent steps reuse the most recently written value (this one).
    if (step == 6) then
      pressure = [(rank + j*0.2, j=0, ncells-1)]
      call file%write_cell_data(prs_var, pressure)
    end if

    ! Write point-centered velocity for this step.
    velocity = 0.1*sin(2*t) * x
    call file%write_point_data(vel_var, velocity)

    ! Return file to a complete, readable state until start of next step.
    call file%finalize_time_step
    call file%flush
  end do
  
  call file%close

#ifdef USE_MPI
  call MPI_Finalize(stat)
#endif

end program
