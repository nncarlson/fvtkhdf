!! Fixture producer for the HyperTreeGrid writer. The trees, ownership, and
!! data follow simple rules that vtkhdf_htg_check.py reproduces independently:
!!
!! * Tree g (0-based) in topology phase p has pattern mod(g+p,3): a single
!!   leaf; a refined root; or a refined root whose first and last children
!!   are also refined.
!! * Tree g belongs to rank mod(g+p,nactive), where nactive is the number of
!!   ranks except in phase 2, which leaves the last of several ranks empty.
!! * A leaf with (0-based) breadth-first index l in its tree is masked when
!!   mod(l,5) == 4.
!! * Cell data values encode the tree, the cell's index in its tree, and the
!!   step; coordinates encode the axis, the point index, and a geometry phase.

program vtkhdf_htg_test

  use,intrinsic :: iso_fortran_env, only: int8, int32, int64, real32, real64
  use vtkhdf_test_env
  use vtkhdf_htg_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  integer :: rank, nproc

  call test_env_init(rank, nproc)

  call write_static('htg_static.vtkhdf', [4,3,2], 2, with_mask=.false.)
  call write_static('htg_static_mask.vtkhdf', [3,1,3], 3, with_mask=.true.)
  call write_coarse
  call write_fixed
  call write_moving
  call write_dynamic

  call test_env_finalize()

contains

  subroutine open_file(file, name, mode)
    type(vtkhdf_htg_file), intent(out) :: file
    character(*), intent(in) :: name
    integer, intent(in) :: mode
    integer :: stat
    character(:), allocatable :: errmsg
#ifdef USE_MPI
    call file%create(name, MPI_COMM_WORLD, stat, errmsg, mode=mode)
#else
    call file%create(name, stat, errmsg, mode=mode)
#endif
    if (stat /= 0) error stop errmsg
  end subroutine

  !! Rank-zero coordinate vector for an axis; other ranks get an empty vector.
  function axis_coordinates(axis, n, phase) result(c)
    integer, intent(in) :: axis, n, phase
    real(real64), allocatable :: c(:)
    integer :: i
    if (rank == 0) then
      c = [(axis*i + 0.25_real64*i*i + 0.5_real64*phase*axis, i = 0, n-1)]
    else
      allocate(c(0))
    end if
  end function

  integer function child_count(dims, branch_factor) result(nchild)
    integer, intent(in) :: dims(3), branch_factor
    nchild = branch_factor**count(dims > 1)
  end function

  !! This rank's trees for a topology phase, with the tree (GID) and index
  !! within its tree (LID) of each cell for generating data.
  subroutine local_trees(dims, branch_factor, phase, tree_ids, xtree, refined, mask, gid, lid)
    integer, intent(in) :: dims(3), branch_factor, phase
    integer, allocatable, intent(out) :: tree_ids(:), xtree(:), gid(:), lid(:)
    logical, allocatable, intent(out) :: refined(:), mask(:)
    integer :: g, l, n, nchild, nactive
    logical, allocatable :: r(:)
    nchild = child_count(dims, branch_factor)
    nactive = nproc
    if (phase == 2) nactive = max(nproc-1, 1)
    allocate(tree_ids(0), gid(0), lid(0), refined(0))
    xtree = [1]
    do g = 0, product(max(dims-1,1)) - 1
      if (modulo(g+phase, nactive) /= rank) cycle
      select case (modulo(g+phase, 3))
      case (0)
        r = [.false.]
      case (1)
        r = [.true., (.false., l = 1, nchild)]
      case (2)
        allocate(r(1+3*nchild), source=.false.)
        r([1, 2, nchild+1]) = .true.
      end select
      n = size(r)
      tree_ids = [tree_ids, g+1]
      xtree = [xtree, xtree(size(xtree)) + n]
      refined = [refined, r]
      gid = [gid, (g, l = 0, n-1)]
      lid = [lid, (l, l = 0, n-1)]
      deallocate(r)
    end do
    mask = (.not.refined) .and. (modulo(lid, 5) == 4)
  end subroutine

  subroutine write_static(name, dims, branch_factor, with_mask)
    character(*), intent(in) :: name
    integer, intent(in) :: dims(3), branch_factor
    logical, intent(in) :: with_mask
    type(vtkhdf_htg_file) :: file
    integer, allocatable :: tree_ids(:), xtree(:), gid(:), lid(:)
    logical, allocatable :: refined(:), mask(:)
    real(real32), allocatable :: vec(:,:)
    integer :: j
    call open_file(file, name, VTKHDF_STATIC)
    call file%write_mesh(dims, branch_factor=branch_factor)
    call local_trees(dims, branch_factor, 0, tree_ids, xtree, refined, mask, gid, lid)
    if (with_mask) then ! topology first, int64 indexing, real32 coordinates
      call file%write_mesh_topology(int(tree_ids,int64), int(xtree,int64), refined, mask)
      call file%write_mesh_geometry(real(axis_coordinates(1,dims(1),0),real32), &
          real(axis_coordinates(2,dims(2),0),real32), real(axis_coordinates(3,dims(3),0),real32))
    else ! geometry first
      call file%write_mesh_geometry(axis_coordinates(1,dims(1),0), &
          axis_coordinates(2,dims(2),0), axis_coordinates(3,dims(3),0))
      call file%write_mesh_topology(tree_ids, xtree, refined)
    end if
    call file%write_cell_data('value', 1000.0_real64*gid + lid, attribute='Scalars')
    allocate(vec(3,size(gid)))
    do j = 1, size(gid)
      vec(:,j) = [real(gid(j)), real(lid(j)), 0.0]
    end do
    call file%write_cell_data('vec', vec)
    call file%write_cell_data('ids', int(1000*gid + lid, int32))
    call file%write_cell_data('big ids', int(1000*gid + lid, int64))
    call file%write_cell_data('flag', int(modulo(lid,2), int8))
    call file%write_field_data('field-value', 42.0_real64)
    call file%close
  end subroutine

  !! No tree is refined: there are no descriptor bits at all.
  subroutine write_coarse
    type(vtkhdf_htg_file) :: file
    integer, parameter :: dims(3) = [5,1,1]
    integer, allocatable :: tree_ids(:), xtree(:)
    logical, allocatable :: refined(:)
    integer :: g
    call open_file(file, 'htg_coarse.vtkhdf', VTKHDF_STATIC)
    call file%write_mesh(dims)
    call file%write_mesh_geometry(axis_coordinates(1,dims(1),0), &
        axis_coordinates(2,dims(2),0), axis_coordinates(3,dims(3),0))
    tree_ids = pack([(g+1, g = 0, 3)], [(modulo(g,nproc) == rank, g = 0, 3)])
    xtree = [(g, g = 1, size(tree_ids)+1)]
    allocate(refined(size(tree_ids)), source=.false.)
    call file%write_mesh_topology(tree_ids, xtree, refined)
    call file%write_cell_data('value', 1000.0_real64*(tree_ids-1))
    call file%close
  end subroutine

  !! Fixed mesh: static and temporal cell data, with omission and reuse.
  subroutine write_fixed
    type(vtkhdf_htg_file) :: file
    type(vtkhdf_cell_data_handle) :: value_var, vec_var
    type(vtkhdf_field_data_handle) :: time_var
    integer, parameter :: dims(3) = [3,4,1]
    integer, allocatable :: tree_ids(:), xtree(:), gid(:), lid(:)
    logical, allocatable :: refined(:), mask(:)
    real(real32), allocatable :: vec(:,:)
    integer :: step, j
    call open_file(file, 'htg_fixed.vtkhdf', VTKHDF_FIXED_MESH)
    call file%write_mesh(dims)
    call local_trees(dims, 2, 0, tree_ids, xtree, refined, mask, gid, lid)
    call file%write_mesh_topology(tree_ids, xtree, refined, mask)
    call file%write_mesh_geometry(axis_coordinates(1,dims(1),0), &
        axis_coordinates(2,dims(2),0), axis_coordinates(3,dims(3),0))
    call file%write_cell_data('ids', int(1000*gid + lid, int32))
    value_var = file%register_temporal_cell_data('value', 0.0_real64, attribute='Scalars')
    vec_var = file%register_temporal_cell_data('vec', [0.0_real32, 0.0_real32, 0.0_real32])
    time_var = file%register_temporal_field_data('time', 0.0_real64)
    allocate(vec(3,size(gid)))
    do step = 0, 3
      call file%start_time_step(real(step,real64))
      call file%write_field_data(time_var, real(step,real64))
      if (step /= 2) call file%write_cell_data(value_var, 1000.0_real64*gid + lid + 0.5_real64*step)
      if (step == 0) then
        do j = 1, size(gid)
          vec(:,j) = [real(gid(j)), real(lid(j)), real(step)]
        end do
        call file%write_cell_data(vec_var, vec)
      end if
      call file%finalize_time_step
      call file%flush
    end do
    call file%close
  end subroutine

  !! Moving mesh: new coordinates at steps 0 and 2.
  subroutine write_moving
    type(vtkhdf_htg_file) :: file
    type(vtkhdf_cell_data_handle) :: value_var
    integer, parameter :: dims(3) = [3,3,2]
    integer, allocatable :: tree_ids(:), xtree(:), gid(:), lid(:)
    logical, allocatable :: refined(:), mask(:)
    integer :: step
    call open_file(file, 'htg_moving.vtkhdf', VTKHDF_MOVING_MESH)
    call file%write_mesh(dims)
    call local_trees(dims, 2, 0, tree_ids, xtree, refined, mask, gid, lid)
    value_var = file%register_temporal_cell_data('value', 0.0_real64)
    do step = 0, 3
      call file%start_time_step(real(step,real32))
      if (step == 0) call file%write_mesh_topology(tree_ids, xtree, refined)
      if (step == 0 .or. step == 2) call file%write_mesh_geometry(axis_coordinates(1,dims(1),step/2), &
          axis_coordinates(2,dims(2),step/2), axis_coordinates(3,dims(3),step/2))
      if (step == 1) call file%write_cell_data('ids', int(1000*gid + lid, int32))
      call file%write_cell_data(value_var, 1000.0_real64*gid + lid + 0.5_real64*step)
    end do
    call file%close
  end subroutine

  !! Dynamic mesh: topology phases 0,0,1,2 and geometry phases 0,0,0,1. Only
  !! the phase 1 topology supplies a mask. Changing ownership changes the
  !! number of parts, and hence of stored coordinate copies, in MPI runs.
  subroutine write_dynamic
    type(vtkhdf_htg_file) :: file
    type(vtkhdf_cell_data_handle) :: value_var, vec_var
    integer, parameter :: dims(3) = [4,3,2]
    integer, parameter :: topology(0:3) = [0,0,1,2], geometry(0:3) = [0,0,0,1]
    integer, allocatable :: tree_ids(:), xtree(:), gid(:), lid(:)
    logical, allocatable :: refined(:), mask(:)
    real(real32), allocatable :: vec(:,:)
    integer :: step, j
    call open_file(file, 'htg_dynamic.vtkhdf', VTKHDF_DYNAMIC_MESH)
    call file%write_mesh(dims)
    value_var = file%register_temporal_cell_data('value', 0.0_real64)
    vec_var = file%register_temporal_cell_data('vec', [0.0_real32, 0.0_real32, 0.0_real32])
    call file%write_mesh_geometry(axis_coordinates(1,dims(1),0), &
        axis_coordinates(2,dims(2),0), axis_coordinates(3,dims(3),0))
    do step = 0, 3
      call file%start_time_step(0.25_real64*step)
      if (step == 0 .or. topology(step) /= topology(max(step-1,0))) then
        call local_trees(dims, 2, topology(step), tree_ids, xtree, refined, mask, gid, lid)
        if (topology(step) == 1) then
          call file%write_mesh_topology(tree_ids, xtree, refined, mask)
        else
          call file%write_mesh_topology(tree_ids, xtree, refined)
        end if
      end if
      if (step > 0 .and. geometry(step) /= geometry(max(step-1,0))) then
        call file%write_mesh_geometry(axis_coordinates(1,dims(1),geometry(step)), &
            axis_coordinates(2,dims(2),geometry(step)), axis_coordinates(3,dims(3),geometry(step)))
      end if
      call file%write_cell_data(value_var, 1000.0_real64*gid + lid + 0.5_real64*step)
      if (step /= 1) then ! step 1 reuses the step 0 values
        if (allocated(vec)) deallocate(vec)
        allocate(vec(3,size(gid)))
        do j = 1, size(gid)
          vec(:,j) = [real(gid(j)), real(lid(j)), real(step)]
        end do
        call file%write_cell_data(vec_var, vec)
      end if
    end do
    call file%close
  end subroutine

end program vtkhdf_htg_test
