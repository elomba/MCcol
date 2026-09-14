!===============================================================================
! Module: Util
!
! Purpose:
!   Provides high-resolution wall-clock / CPU timing utilities for performance
!   profiling and runtime monitoring in gpMC.
!===============================================================================
Module Util
  Use set_precision
  Implicit None
Contains

  !-----------------------------------------------------------------------------
  ! Function: cputime
  !
  ! Purpose:
  !   Returns the elapsed wall/system time in seconds based on the Fortran
  !   intrinsic `system_clock`.
  !
  ! Returns:
  !   cputime (Real(wp)) : Elapsed time in seconds.
  !-----------------------------------------------------------------------------
  Function cputime()
    Implicit None
    Integer :: count, count_rate
    Real (wp) :: cputime
    Call System_clock(count, count_rate)
    cputime = Real(count, wp) / Real(count_rate, wp)
  End Function cputime

End Module Util
