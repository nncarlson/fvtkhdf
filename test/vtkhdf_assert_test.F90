#include "vtkhdf_assert.inc"
program vtkhdf_assert_test
  implicit none
  character(4) :: scenario
  logical :: condition

  interface
    subroutine vtkhdf_assert(file, line)
      character(*), intent(in) :: file
      integer, intent(in) :: line
    end subroutine
  end interface

  call get_command_argument(1, scenario)
  if (scenario /= 'pass' .and. scenario /= 'fail') error stop 'invalid scenario'
  condition = scenario == 'pass'
  INSIST(condition)
  ! End normally so the checker detects a failed assertion that returns.
end program
