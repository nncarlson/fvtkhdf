program vtkhdf_mb_test

  use,intrinsic :: iso_fortran_env, only: r8 => real64, int8
  use vtkhdf_mb_file_type
  use vtkhdf_vtk_cell_types, only: VTK_HEXAHEDRON
  use mpi_f08
  implicit none

  real(r8), allocatable :: points(:,:)
  real(r8), allocatable :: pressure(:), temperature(:), velocity(:,:)
  integer,  allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  real(r8) :: time
  character(:), allocatable :: errmsg
  integer :: stat, j, nproc, rank
  integer :: npoints_liquid, ncells_liquid, npoints_solid, ncells_solid

  type(vtkhdf_mb_file) :: vizfile
  type(vtkhdf_block_handle) :: bliq, bsol

  call MPI_Init(stat)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, stat)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)

  if (nproc < 2) then
    call MPI_Finalize(stat)
    error stop 'must be run with at least 2 processes'
  end if

  !! Create the file
  call vizfile%create('mb_demo.vtkhdf', MPI_COMM_WORLD, stat, errmsg)
  if (stat /= 0) error stop errmsg

  !! We have 2 UnstructuredGrid blocks: "liquid" and "solid".
  !! Every rank has a piece of the liquid block mesh.
  !! Every rank has a piece of the solid block mesh, EXCEPT rank 0.
  !! Both blocks have a time-independent temperature.
  !! The liquid block also has time-dependent pressure and velocity.
  !! The solid block is static; the liquid block is temporal.

  !! Get the local piece of the liquid block mesh.
  call get_liquid_mesh_data(points, cnode, xcnode, types)

  !! Add the liquid block and write the local mesh piece (COLLECTIVE!)
  bliq = vizfile%add_block('liquid', is_temporal=.true.)
  call vizfile%write_mesh(bliq, points, cnode, xcnode, types)

  !! local mesh sizes
  npoints_liquid = size(points,dim=2)
  ncells_liquid  = size(xcnode) - 1

  !! Get the local piece of the solid block mesh.
  !! NB: rank 0 gets a 0-sized empty mesh!
  call get_solid_mesh_data(points, cnode, xcnode, types)

  !! Add the solid block and write the local mesh piece (COLLECTIVE!)
  !! Rank 0 must participate with its 0-sized mesh!
  !NB: A bug in the current reader requires it to be temporal.
  !bsol = vizfile%add_block('solid')
  bsol = vizfile%add_block('solid', is_temporal=.true.)
  call vizfile%write_mesh(bsol, points, cnode, xcnode, types)

  !! local mesh sizes
  npoints_solid = size(points,dim=2)
  ncells_solid  = size(xcnode) - 1

  !! Register the time-dependent cell-centered pressure and point-centered
  !! velocity for the liquid block. (COLLECTIVE!)
  allocate(pressure(ncells_liquid), velocity(3,npoints_liquid))
  associate (scalar_mold => pressure(1), vector_mold => velocity(:,1))
    call vizfile%register_temporal_cell_data(bliq, 'pressure', scalar_mold)
    call vizfile%register_temporal_point_data(bliq, 'velocity', vector_mold)
  end associate

  !! Start simulation time stepping
  do j = 0, 10
    time = j*0.1_r8

    !! Start the time step
    call vizfile%write_time_step(time) ! COLLECTIVE!

    !! Generate some arbitrary time-dependent data and write it. (COLLECTIVE!)
    pressure = cos(time) + rank
    velocity = spread([cos(time+rank),sin(time+rank),1.0_r8],dim=2,ncopies=npoints_liquid)
    call vizfile%write_temporal_cell_data(bliq, 'pressure', pressure)
    call vizfile%write_temporal_point_data(bliq, 'velocity', velocity)
  end do

  !! Write the time-independent point-centered temperature for both blocks
  !! (COLLECTIVE!) Rank 0 must participate in the call for the solid block
  !! with its 0-sized temperature data. Static data can be written at any
  !! time after the mesh, but its name must be unique among data of its mesh
  !! entity type.
  temperature = spread(rank, dim=1, ncopies=npoints_liquid)
  call vizfile%write_point_data(bliq, 'temperature', temperature)

  temperature = spread(rank, dim=1, ncopies=npoints_solid)
  call vizfile%write_point_data(bsol, 'temperature', temperature)

  call vizfile%close
  call MPI_Finalize(stat)

contains

  ! A single hex cell mesh
  subroutine get_unit_mesh(points, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    points = 0.8*reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1], [3,8])
    cnode = [1,2,3,4,5,6,7,8]
    xcnode = [1,9]
    types = [VTK_HEXAHEDRON]
  end subroutine

  ! Each rank gets one hex cell
  subroutine get_liquid_mesh_data(points, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    call get_unit_mesh(points, cnode, xcnode, types)
    points(1,:) = points(1,:) + rank ! spread in x
  end subroutine

  ! Each rank gets one hex cell, EXCEPT rank 0
  subroutine get_solid_mesh_data(points, cnode, xcnode, types)
    real(r8), allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    call get_unit_mesh(points, cnode, xcnode, types)
    points(1,:) = points(1,:) + rank ! spread in x
    points(2,:) = points(2,:) + 1    ! common shift in y
    if (rank == 0) then ! 0-sized mesh
      points = reshape([real(r8)::], [3,0])
      cnode = [integer::]
      xcnode = [1]
      types = [integer(int8)::]
    end if
  end subroutine

end program
