program vtkhdf_mb_test

  use,intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_mb_file_type
  use vtkhdf_vtk_cell_types, only: VTK_TETRA
  implicit none

  real(r8), allocatable :: points(:,:)
  real(r8), allocatable :: pressure(:), temperature(:), velocity(:,:)
  integer, allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  real(r8) :: time
  character(:), allocatable :: errmsg
  integer :: stat, j, npoints, ncells

  type(vtkhdf_mb_file) :: vizfile

  !! Create the file
  call vizfile%create('mb_demo.vtkhdf', stat, errmsg)
  if (stat /= 0) error stop errmsg

  !! We have 2 UnstructuredGrid blocks: "liquid" and "solid"
  !! Both blocks have a time-independent temperature.
  !! The liquid block also has time-dependent pressure and velocity.
  !! The solid block is static; the liquid block is temporal.

  !! Generate the mesh for the liquid block; the mesh for the
  !! solid block with be a non-overlapping shift of the liquid block.
  call get_mesh_data(points, cnode, xcnode, types)

  npoints = size(points,dim=2)
  ncells  = size(xcnode) - 1

  !! Add the blocks to the file and write their meshes.
  call vizfile%add_block('liquid', stat, errmsg, is_temporal=.true.)
  if (stat /= 0) error stop errmsg
  call vizfile%write_mesh('liquid', points, cnode, xcnode, types)

  !! Shift the points right to get the solid block mesh
  points(1,:) = points(1,:) + 1.2_r8
  !NB: A bug in the current reader that requires it to be temporal.
  !call vizfile%add_block('solid', stat, errmsg)
  call vizfile%add_block('solid', stat, errmsg, is_temporal=.true.)
  if (stat /= 0) error stop errmsg
  call vizfile%write_mesh('solid', points, cnode, xcnode, types)

  !! Register the time-dependent cell-centered pressure and point-centered
  !! velocity for the liquid block.
  allocate(pressure(ncells), velocity(3,npoints))
  associate (scalar_mold => pressure(1), vector_mold => velocity(:,1))
    call vizfile%register_temporal_cell_data('liquid', 'pressure', scalar_mold)
    call vizfile%register_temporal_point_data('liquid', 'velocity', vector_mold)
  end associate

  !! Start simulation time stepping
  do j = 0, 10
    time = j*0.1_r8

    !! Start the time step
    call vizfile%write_time_step(time)

    !! Generate some arbitrary time-dependent data and write it.
    pressure = cos(time)
    velocity = spread([cos(time),sin(time),1.0_r8],dim=2,ncopies=npoints)
    call vizfile%write_temporal_cell_data('liquid', 'pressure', pressure)
    call vizfile%write_temporal_point_data('liquid', 'velocity', velocity)
  end do

  !! We have time-independent (static) point-centered temperature in
  !! both blocks. Static data can be written at any time after the mesh,
  !! but its name must be unique among data of its mesh entity type.
  allocate(temperature(npoints), source=1.0_r8)
  call vizfile%write_point_data('liquid', 'temperature', temperature)
  call vizfile%write_point_data('solid',  'temperature', 1+temperature)

  call vizfile%close

contains

  ! A 5-tet subdivision of a unit cube.
  subroutine get_mesh_data(points, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    points = reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1], shape=[3,8])
    cnode = [1,2,4,5, 2,3,4,7, 2,5,6,7, 4,5,7,8, 2,4,5,7]
    xcnode = [1,5,9,13,17,21]
    types = spread(VTK_TETRA, dim=1, ncopies=5)
  end subroutine

end program
