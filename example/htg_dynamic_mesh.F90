program htg_dynamic_mesh_example

  use,intrinsic :: iso_fortran_env, only: int8
  use vtkhdf_htg_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  integer, parameter :: nx = 6, ny = 4      ! root cells in x and y
  integer, parameter :: max_level = 4       ! deepest refinement level
  real, parameter :: radius = 1.2

  type(vtkhdf_htg_file) :: file
  type(vtkhdf_cell_data_handle) :: level_var, dist_var
  type(vtkhdf_field_data_handle) :: center_var

  integer, allocatable :: tree_ids(:), xtree(:)
  logical, allocatable :: refined(:)
  integer(int8), allocatable :: level(:)
  real, allocatable :: dist(:), x(:), y(:), z(:)
  character(:), allocatable :: errmsg
  integer :: stat, nproc, rank, step, i, j, tree
  real :: center(2), t

#ifdef USE_MPI
  call MPI_Init(stat)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, stat)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)
#else
  nproc = 1
  rank = 0
#endif

  ! Create the file with support for a time-dependent tree refinement.
#ifdef USE_MPI
  call file%create('htg_dynamic_mesh.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=VTKHDF_DYNAMIC_MESH)
  if (stat /= 0) then
    call MPI_Finalize(stat)
    error stop errmsg
  end if
#else
  call file%create('htg_dynamic_mesh.vtkhdf', stat, errmsg, mode=VTKHDF_DYNAMIC_MESH)
  if (stat /= 0) error stop errmsg
#endif

  ! Define the 2D grid of root cells by its point dimensions; z is a singleton
  ! axis, so every refined cell has 4 children (a quadtree).
  call file%write_mesh([nx+1, ny+1, 1])

  ! Rank zero supplies the coordinate vectors; other ranks pass empty arrays.
  if (rank == 0) then
    x = [(real(i), i = 0, nx)]
    y = [(real(j), j = 0, ny)]
    z = [0.0]
  else
    allocate(x(0), y(0), z(0))
  end if
  call file%write_mesh_geometry(x, y, z)

  ! With a changing topology all cell data is temporal and must be registered.
  level_var = file%register_temporal_cell_data('level', 0_int8)
  dist_var  = file%register_temporal_cell_data('distance', 0.0, attribute='Scalars')
  center_var = file%register_temporal_field_data('center', 0.0)

  do step = 0, 8
    t = step / 8.0
    call file%start_time_step(t)

    ! A circle moves across the grid; cells crossed by it are refined.
    center = [0.5 + 5.0*t, 1.0 + 2.0*t]
    call file%write_field_data(center_var, center, as_vector=.true.)

    ! Each rank describes its own trees. Root cells are numbered from 1 with
    ! x varying fastest; here they are dealt to the ranks in turn. The cells
    ! of each tree are listed level by level, starting with its root cell.
    allocate(tree_ids(0), refined(0), level(0), dist(0))
    xtree = [1]
    do j = 1, ny
      do i = 1, nx
        tree = i + (j-1)*nx
        if (modulo(tree-1, nproc) /= rank) cycle
        call append_tree(real(i-1), real(j-1))
        tree_ids = [tree_ids, tree]
        xtree = [xtree, size(refined) + 1]
      end do
    end do

    ! Write the new topology before the cell data for this step.
    call file%write_mesh_topology(tree_ids, xtree, refined)
    call file%write_cell_data(level_var, level)
    call file%write_cell_data(dist_var, dist)

    ! Return file to a complete, readable state until start of next step.
    call file%finalize_time_step
    call file%flush
    deallocate(tree_ids, xtree, refined, level, dist)
  end do

  call file%close

#ifdef USE_MPI
  call MPI_Finalize(stat)
#endif

contains

  !! Append the cells of the tree whose unit root cell has lower corner
  !! (X0,Y0), in breadth-first order: all cells of one level precede those of
  !! the next, and the 4 children of a refined cell are contiguous, ordered
  !! with x varying fastest, and follow the order of their parents.
  subroutine append_tree(x0, y0)
    real, intent(in) :: x0, y0
    real, allocatable :: xlo(:), ylo(:), xnext(:), ynext(:)
    logical, allocatable :: split(:)
    real :: h, d
    integer :: l, n, c
    xlo = [x0]
    ylo = [y0]
    h = 1.0
    do l = 0, max_level
      allocate(split(size(xlo)), xnext(0), ynext(0))
      do n = 1, size(xlo)
        d = norm2([xlo(n), ylo(n)] + 0.5*h - center) - radius ! signed distance at the cell center
        split(n) = (l < max_level) .and. (abs(d) < 0.75*h)
        dist = [dist, d]
        if (split(n)) then
          xnext = [xnext, (xlo(n) + 0.5*h*modulo(c,2), c = 0, 3)]
          ynext = [ynext, (ylo(n) + 0.5*h*(c/2), c = 0, 3)]
        end if
      end do
      refined = [refined, split]
      level = [level, spread(int(l,int8), 1, size(xlo))]
      deallocate(split)
      if (size(xnext) == 0) exit
      call move_alloc(xnext, xlo)
      call move_alloc(ynext, ylo)
      h = 0.5*h
    end do
  end subroutine

end program
