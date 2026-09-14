!===============================================================================
! Module: set_precision
!
! Purpose:
!   Centralizes floating-point kind parameters for the gpMC package.
!   Allows switching between single, double, and extended precision across
!   all modules by modifying the working precision parameter `wp`.
!
! Precision Constants:
!   skind : Kind for standard 32-bit single-precision IEEE floating point.
!   dkind : Kind for standard 64-bit double-precision IEEE floating point.
!   wp    : Working precision kind used throughout the simulation code.
!           Set to `skind` or `dkind` as desired.
!
! Note on Quadruple Precision:
!   Quadruple precision (`qkind`) is non-standard and varies by compiler.
!   Guidance for Intel, IBM, and NAG compilers is provided below.
!===============================================================================
MODULE set_precision
  IMPLICIT NONE
  INTRINSIC KIND

  ! Standard IEEE floating-point precisions
  INTEGER, PARAMETER :: skind = KIND(0.0E0)  ! Single precision (32-bit)
  INTEGER, PARAMETER :: dkind = KIND(0.0D0)  ! Double precision (64-bit)

  ! Working precision for gpMC simulation calculations
  INTEGER, PARAMETER :: wp = skind

  !-----------------------------------------------------------------------------
  ! Optional / compiler-dependent quad precision:
  ! Intel/IBM:
  !   INTEGER, PARAMETER :: qkind = KIND(0.0Q0)
  ! NAG compiler:
  !   USE, INTRINSIC :: f90_kind
  !   INTEGER, PARAMETER :: qkind = quad
  !-----------------------------------------------------------------------------

END MODULE set_precision
