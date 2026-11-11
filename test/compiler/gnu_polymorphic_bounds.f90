program repro
  implicit none
  type :: t
    integer :: value
  end type
  type(t) :: a(0:1) = [t(42), t(99)]
  call check(a)
contains
  subroutine check(x)
    class(t), intent(in) :: x(:)
    print *, x(1)%value ! Expected: 42
  end subroutine
end program
