program camb_p
  use planck_config
  use ls_misc_utils
  use camb_main

  implicit none

  integer(i4b) :: retval

  retval = camb_module2(getcmdline())
  if (retval/=0) call exit_with_status(retval)
end program camb_p
