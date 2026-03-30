program mb_example1

  use,intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_mb_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  real, allocatable :: points(:,:)
  real, allocatable :: pressure(:), temperature(:), velocity(:,:)
  integer,  allocatable :: cnode(:), xcnode(:)
  integer(int8), allocatable :: types(:)
  character(:), allocatable :: errmsg
  integer :: stat, nproc, rank, step, j
  integer :: npoints_liquid, ncells_liquid, npoints_solid, ncells_solid
  real :: t

  type(vtkhdf_mb_file) :: file
  type(vtkhdf_block_handle) :: liq_blk, sol_blk
  type(vtkhdf_cell_data_handle) :: prs_var
  type(vtkhdf_point_data_handle) :: vel_var

#ifdef USE_MPI
  call MPI_Init(stat)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, stat)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)
#else
  nproc = 1
  rank = 0
#endif

  ! Create the file
#ifdef USE_MPI
  call file%create('example.vtkhdf', MPI_COMM_WORLD, stat, errmsg)
  if (stat /= 0) then
    call MPI_Finalize(stat)
    error stop errmsg
  end if
#else
  call file%create('example.vtkhdf', stat, errmsg)
  if (stat /= 0) error stop errmsg
#endif

  ! We have 2 UnstructuredGrid blocks: "liquid" and "solid":
  ! - Every rank has a piece of the liquid block mesh.
  ! - Every rank has a piece of the solid block mesh, except the last
  !   when nproc > 1.
  ! - Both blocks have a time-independent temperature.
  ! - The liquid block also has time-dependent pressure and velocity.
  ! - The solid block is UG_STATIC; the liquid block is UG_FIXED_MESH.

  ! Get the local piece of the liquid block mesh.
  call get_liquid_mesh_data(points, cnode, xcnode, types)

  ! Add the liquid block and write the local mesh piece.
  liq_blk = file%add_block('liquid', mode=UG_FIXED_MESH)
  call file%write_mesh(liq_blk, points, cnode, xcnode, types)

  ! local mesh sizes
  npoints_liquid = size(points,dim=2)
  ncells_liquid  = size(xcnode) - 1

  ! Get the local piece of the solid block mesh.
  ! NB: last rank gets a 0-sized empty mesh when nproc > 1
  call get_solid_mesh_data(points, cnode, xcnode, types)

  ! Add the solid block and write the local mesh piece.
  ! Collective: rank with 0-sized mesh must participate!
  ! NB: Paraview reader bug requires it to be temporal; see issue #27.
  !sol_blk = file%add_block('solid', mode=UG_STATIC)
  sol_blk = file%add_block('solid', mode=UG_FIXED_MESH)
  call file%write_mesh(sol_blk, points, cnode, xcnode, types)

  ! local mesh sizes
  npoints_solid = size(points,dim=2)
  ncells_solid  = size(xcnode) - 1

  ! Register the time-dependent cell-centered pressure and point-centered
  ! velocity for the liquid block.
  allocate(pressure(ncells_liquid), velocity(3,npoints_liquid))
  associate (scalar_mold => pressure(1), vector_mold => velocity(:,1))
    prs_var = file%register_temporal_cell_data(liq_blk, 'pressure', scalar_mold)
    vel_var = file%register_temporal_point_data(liq_blk, 'velocity', vector_mold)
  end associate

  ! Start time stepping.
  do step = 0, 10
    t = step*0.1
    call file%start_time_step(t)

    ! Write cell-centered pressure for the liquid block.
    pressure = [(rank + j*t, j=1, ncells_liquid)]
    call file%write_cell_data(liq_blk, prs_var, pressure)

    ! Write point-centered velocity for the liquid block.
    velocity = spread([cos(t+rank),sin(t+rank),0.0],dim=2,ncopies=npoints_liquid)
    call file%write_point_data(liq_blk, vel_var, velocity)

    ! Return file to a complete, readable state until start of next step.
    call file%finalize_time_step
    call file%flush
  end do

  !! Write the time-independent point-centered temperature for both blocks
  !! Rank with 0-sized mesh must participate in the call for the solid block
  !! with its 0-sized temperature data.
  temperature = spread(rank, dim=1, ncopies=npoints_liquid)
  call file%write_point_data(liq_blk, 'temperature', temperature)

  temperature = spread(rank, dim=1, ncopies=npoints_solid)
  call file%write_point_data(sol_blk, 'temperature', temperature)

  call file%close

#ifdef USE_MPI
  call MPI_Finalize(stat)
#endif

contains

  ! A single hex cell mesh
  subroutine get_unit_mesh(points, cnode, xcnode, types)
    use vtkhdf_vtk_cell_types, only: VTK_HEXAHEDRON
    real, allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    points = reshape([0,0,0, 1,0,0, 1,1,0, 0,1,0, &
                      0,0,1, 1,0,1, 1,1,1, 0,1,1], [3,8])
    cnode = [1,2,3,4,5,6,7,8]
    xcnode = [1,9]
    types = [VTK_HEXAHEDRON]
  end subroutine

  ! Liquid mesh: each rank gets one hex cell
  subroutine get_liquid_mesh_data(points, cnode, xcnode, types)
    real, allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    call get_unit_mesh(points, cnode, xcnode, types)
    points(1,:) = points(1,:) + rank ! spread in x
  end subroutine

  ! Solid mesh: each rank gets one hex cell, except last when nproc > 1.
  subroutine get_solid_mesh_data(points, cnode, xcnode, types)
    real, allocatable, intent(out) :: points(:,:)
    integer, allocatable, intent(out) :: cnode(:), xcnode(:)
    integer(int8), allocatable :: types(:)
    call get_unit_mesh(points, cnode, xcnode, types)
    points(1,:) = points(1,:) + rank ! spread in x
    points(2,:) = points(2,:) + 1    ! common shift in y
    if (rank >= max(1,nproc-1)) then ! 0-sized mesh
      points = reshape([real::], [3,0])
      cnode = [integer::]
      xcnode = [1]
      types = [integer(int8)::]
    end if
  end subroutine

end program
