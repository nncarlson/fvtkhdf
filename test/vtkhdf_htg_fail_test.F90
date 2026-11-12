!! Expected-failure contract checks for the HyperTreeGrid writer. Each case
!! must terminate through the assertion machinery; reaching the end of the
!! program is the test failure.

program vtkhdf_htg_fail_test

  use,intrinsic :: iso_fortran_env, only: real64
  use vtkhdf_test_env
  use vtkhdf_htg_file_type
#ifdef USE_MPI
  use mpi
#endif
  implicit none

  type(vtkhdf_htg_file) :: file
  type(vtkhdf_cell_data_handle) :: var
  character(64) :: case_name
  character(:), allocatable :: errmsg
  integer :: rank, nproc, stat, mode
  real(real64), parameter :: x(3) = [0.0_real64, 1.0_real64, 2.0_real64], y(2) = [0.0_real64, 1.0_real64]
  real(real64), parameter :: z(1) = [0.0_real64]

  !! Two quadtrees: a single leaf and a refined root (6 cells).
  integer, parameter :: tree_ids(2) = [1,2], xtree(3) = [1,2,7]
  logical, parameter :: refined(6) = [.false., .true., .false., .false., .false., .false.]

  call test_env_init(rank, nproc)
  call get_command_argument(1, case_name)

  select case (case_name)
  case ('geometry_change', 'late_register')
    mode = VTKHDF_FIXED_MESH
  case ('mesh_after_data')
    mode = VTKHDF_MOVING_MESH
  case ('missing_data', 'static_data')
    mode = VTKHDF_DYNAMIC_MESH
  case default
    mode = VTKHDF_STATIC
  end select

#ifdef USE_MPI
  call file%create('htg_fail_'//trim(case_name)//'.vtkhdf', MPI_COMM_WORLD, stat, errmsg, mode=mode)
#else
  call file%create('htg_fail_'//trim(case_name)//'.vtkhdf', stat, errmsg, mode=mode)
#endif
  if (stat /= 0) error stop errmsg

  select case (case_name)
  case ('no_grid') ! root grid must be defined first
    call file%write_mesh_topology(tree_ids, xtree, refined)
  case ('bad_dims')
    call file%write_mesh([1,1,1])
  case ('bad_branch_factor')
    call file%write_mesh([3,2,1], branch_factor=4)
  case ('grid_twice')
    call file%write_mesh([3,2,1])
    call file%write_mesh([3,2,1])
  case ('axis_length')
    call file%write_mesh([3,2,1])
    call file%write_mesh_geometry(x, x, z)
  case ('tree_id_range')
    call file%write_mesh([3,2,1])
    call file%write_mesh_topology([1,3], xtree, refined)
  case ('missing_cells') ! refined root without all of its children
    call file%write_mesh([3,2,1])
    call file%write_mesh_topology(tree_ids, [1,2,6], refined(:5))
  case ('extra_cells') ! cells beyond the leaves of the last level
    call file%write_mesh([3,2,1])
    call file%write_mesh_topology(tree_ids, [1,3,8], [.false., refined])
  case ('masked_refined')
    call file%write_mesh([3,2,1])
    call file%write_mesh_topology(tree_ids, xtree, refined, mask=refined)
  case ('topology_twice')
    call file%write_mesh([3,2,1])
    call file%write_mesh_topology(tree_ids, xtree, refined)
    call file%write_mesh_topology(tree_ids, xtree, refined)
  case ('data_before_topology')
    call file%write_mesh([3,2,1])
    call file%write_mesh_geometry(x, y, z)
    call file%write_cell_data('value', [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
  case ('data_before_geometry')
    call file%write_mesh([3,2,1])
    call file%write_mesh_topology(tree_ids, xtree, refined)
    call file%write_cell_data('value', [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
  case ('wrong_cell_count')
    call write_complete_mesh
    call file%write_cell_data('value', [1.0, 2.0, 3.0, 4.0, 5.0])
  case ('geometry_change') ! fixed mesh geometry cannot be replaced
    call write_complete_mesh
    call file%start_time_step(0.0)
    call file%start_time_step(1.0)
    call file%write_mesh_geometry(x, y, z)
  case ('late_register') ! temporal registration precedes time stepping
    call write_complete_mesh
    call file%start_time_step(0.0)
    var = file%register_temporal_cell_data('value', 0.0)
  case ('mesh_after_data')
    call write_complete_mesh
    var = file%register_temporal_cell_data('value', 0.0)
    call file%start_time_step(0.0)
    call file%write_cell_data(var, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    call file%write_mesh_geometry(x, y, z)
  case ('static_data') ! a changing topology precludes static cell data
    call write_complete_mesh
    call file%write_cell_data('value', [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
  case ('missing_data') ! reuse requires an unchanged cell count
    call write_complete_mesh
    var = file%register_temporal_cell_data('value', 0.0)
    call file%start_time_step(0.0)
    call file%write_cell_data(var, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    call file%start_time_step(1.0)
    call file%write_mesh_topology(tree_ids, [1,2,3], refined([1,1]))
    call file%finalize_time_step
  case default
    error stop 'unknown case'
  end select

  call file%close
  call test_env_finalize()

contains

  subroutine write_complete_mesh
    call file%write_mesh([3,2,1])
    call file%write_mesh_geometry(x, y, z)
    call file%write_mesh_topology(tree_ids, xtree, refined)
  end subroutine

end program vtkhdf_htg_fail_test
