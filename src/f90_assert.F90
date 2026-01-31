subroutine f90_assert(file, line)
  
  character(len=*), intent(in) :: file
  integer,          intent(in) :: line

  character(256) :: string
  
  write(string,fmt='(a,i4.4)') "Assertion failed at " // file // ":", line
  error stop trim(string)

end subroutine
