#include "vtkhdf_assert.inc"
program vtkhdf_assert_test
  use, intrinsic :: iso_fortran_env, only: output_unit
  implicit none
  character(4) :: scenario
  logical :: condition
  integer :: unit

  interface
    subroutine vtkhdf_assert(file, line)
      character(*), intent(in) :: file
      integer, intent(in) :: line
    end subroutine
  end interface

  call get_command_argument(1, scenario)
  if (scenario /= 'pass' .and. scenario /= 'fail') error stop 'invalid scenario'
  condition = scenario == 'pass'
  write(output_unit,'(a)') 'ASSERT_BEFORE'
  flush(output_unit)
  INSIST(condition)
  open(newunit=unit, file='after-assertion', status='new', action='write')
  write(unit,'(a)') 'ASSERT_AFTER'
  close(unit)
  write(output_unit,'(a)') 'ASSERT_AFTER'
  flush(output_unit)
end program
