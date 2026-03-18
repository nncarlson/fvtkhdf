program ug_fixed_mesh_example

  use,intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_ug_file_type
  use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE
#ifdef USE_MPI
  use mpi_f08
#endif
  implicit none

  real, allocatable :: points(:,:)
  real, allocatable :: pressure(:), temperature(:), velocity(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  character(:), allocatable :: errmsg
  integer :: stat, nproc, rank, npoints, ncells, step, j
  real :: t

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

  ! Create file with support for time-dependent data but fixed mesh.
#ifdef USE_MPI
  call file%create('example.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=UG_FIXED_MESH)
  if (stat /= 0) then
    call MPI_Finalize(stat)
    error stop errmsg
  end if
#else
  call file%create('example.vtkhdf', stat, errmsg, mode=UG_FIXED_MESH)
  if (stat /= 0) error stop errmsg
#endif

  ! Fixed 2D mesh. The local mesh on each rank is a 2-triangle discretization
  ! of a unit square, suitably x-shifted so they only overlap on a boundary.
  ! The mesh topologies are identical; point coordinates differ. 3D coordinates
  ! are required; we set z=0.

  points = reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0], shape=[3,4])
  points(1,:) = points(1,:) + rank
  cnode = [1,2,3,  1,3,4]
  xcnode = [1,4,7]
  types = spread(VTK_TRIANGLE, dim=1, ncopies=2)

  ! local mesh sizes
  npoints = size(points,dim=2)
  ncells  = size(xcnode) - 1

  ! Each rank writes its part of the mesh collectively.
  call file%write_mesh(points, cnode, xcnode, types)

  ! We have time-dependent cell-centered pressure and point-centered velocity,
  ! and a scalar "cpu time" field data. These need to be registered before
  ! starting time stepping.

  allocate(pressure(ncells), velocity(3,npoints))
  associate (scalar_mold => pressure(1), vector_mold => velocity(:,1))
    prs_var = file%register_temporal_cell_data('pressure', scalar_mold)
    vel_var = file%register_temporal_point_data('velocity', vector_mold)
    cpu_var = file%register_temporal_field_data('cpu time', scalar_mold)
  end associate

  ! Start time stepping.
  do step = 0, 10
    t = step * 0.1
    call file%start_time_step(t)

    ! Field data can be written at any time during the step.
    call file%write_field_data(cpu_var, 5*t)

    ! Write cell-centered pressure on the first step and every third step
    ! thereafter; unwritten steps reuse the most recently written value.
    if (modulo(step,3) == 0) then
      pressure = [(rank + j*t, j=1, ncells)]
      call file%write_cell_data(prs_var, pressure)
    end if

    !! Write point-centered velocity for this step.
    velocity = exp(t)*points
    call file%write_point_data(vel_var, velocity)

    ! Return file to a complete, readable state until start of next step.
    call file%finalize_time_step
    call file%flush
  end do

  ! We have time-independent (static) cell-centered temperature.
  ! Static data can be written at any time after the mesh, but its
  ! name must be unique among data of its mesh entity type.
  allocate(temperature(ncells), source=real(rank))
  call file%write_cell_data('temperature', temperature)

  call file%close

#ifdef USE_MPI
  call MPI_Finalize(stat)
#endif

end program
