program vtkhdf_ug_test

  use,intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_ug_file_type
  use vtkhdf_vtk_cell_types
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  integer :: rank, stat
  character(:), allocatable :: errmsg
  real, allocatable :: points(:,:), p(:)
  integer,  allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  real :: time

  type(vtkhdf_ug_file) :: vizfile
  type(vtkhdf_point_data_handle) :: pvar

#ifdef USE_MPI
  call MPI_Init(stat)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)
  call vizfile%create('dyn_ug_test.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=UG_DYNAMIC_MESH)
#else
  rank = 0
  call vizfile%create('dyn_ug_test.vtkhdf', stat, errmsg, mode=UG_DYNAMIC_MESH)
#endif
  if (stat /= 0) error stop errmsg

  ! Starting mesh. The data is for a basic mesh unit. The global mesh is a
  ! collection of non-overlapping shifts of this basic unit, one for each
  ! rank, right-shifted proportional to the rank.
  call get_mesh1_data(points, cnode, xcnode, types)
  points(1,:) = points(1,:) + rank ! shift right
  call vizfile%write_mesh_topology(cnode, xcnode, types)

  ! Register time-dependent, point-centered pressure dataset
  pvar = vizfile%register_temporal_point_data('pressure', mold=0.0)

  !!!! Step 0: topology already written, write geometry and pressure

  time = 0.0
  call vizfile%start_time_step(time)
  call vizfile%write_mesh_geometry(points)
  call get_pressure_data(points, time, p)
  call vizfile%write_point_data(pvar, p)

  !!!! Step 1: same mesh, write pressure

  time = 0.1
  call vizfile%start_time_step(time)
  call get_pressure_data(points, time, p)
  call vizfile%write_point_data(pvar, p)

  !!!! Step 2: new mesh, write pressure

  time = 0.2
  call vizfile%start_time_step(time)
  call get_mesh2_data(points, cnode, xcnode, types)
  points(1,:) = points(1,:) + rank ! shift right
  call vizfile%write_mesh(points, cnode, xcnode, types)
  call get_pressure_data(points, time, p)
  call vizfile%write_point_data(pvar, p)

  !!!! Step 3: move point, write pressure

  time = 0.3
  call vizfile%start_time_step(time)
  points(2,5) = points(2,5) + 0.1
  call vizfile%write_mesh_geometry(points)
  call get_pressure_data(points, time, p)
  call vizfile%write_point_data(pvar, p)

  !!!! Step 4: move point, don't write pressure

  time = 0.4
  call vizfile%start_time_step(time)
  points(2,5) = points(2,5) + 0.1
  call vizfile%write_mesh_geometry(points)

  call vizfile%close
#ifdef USE_MPI
  call MPI_Finalize(stat)
#endif

contains

  ! Single quad cell mesh.
  subroutine get_mesh1_data(points, cnode, xcnode, types)
    real, allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    points = 0.8 * reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0], shape=[3,4])
    cnode = [1,2,3,4]
    xcnode = [1,5]
    types = [VTK_QUAD]
  end subroutine

  ! Four tri remeshing of mesh1.
  subroutine get_mesh2_data(points, cnode, xcnode, types)
    real, allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    points = 0.4 * reshape([0,0,0, 2,0,0, 2,2,0, 0,2,0, 1,1,0], shape=[3,5])
    cnode = [1,2,5, 2,3,5, 3,4,5, 4,1,5]
    xcnode = [1,4,7,10,13]
    types = spread(VTK_TRIANGLE, dim=1, ncopies=4)
  end subroutine

  subroutine get_pressure_data(points, time, pdata)
    real, intent(in) :: points(:,:), time
    real, allocatable, intent(out) :: pdata(:)
    integer :: j
    allocate(pdata(size(points,dim=2)))
    do j = 1, size(points,dim=2)
      pdata(j) = norm2(points(:,j)) + time
    end do
  end subroutine

end program
