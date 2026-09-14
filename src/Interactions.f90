!===============================================================================
! File: Interactions.f90
!
! Purpose:
!   Evaluates pairwise interatomic potential energies for Lennard-Jones and
!   Morse interactions, with optional real-space Ewald electrostatic contributions
!   and energy shifting at the cutoff distance.
!
! Potentials implemented:
!   1. fpot_LJ:
!      Standard Lennard-Jones (12-6) potential:
!        u_LJ(r) = 4 * epsilon * [ (sigma/r)^12 - (sigma/r)^6 ] - u_cut
!      In code: al = 4*epsilon/kT, bl2 = sigma^2, r6 = (bl2/r2)^3.
!
!   2. fpot_elecLJ:
!      Lennard-Jones (12-6) + Real-space screened Coulomb potential:
!        u(r) = u_LJ(r) + [ q_i * q_j / (4*pi*epsilon_0) ] * erfc(kappa*r) / r
!
!   3. fpot_Morse:
!      Morse potential for covalent / ionic pair interactions:
!        u_Morse(r) = D_e * [ exp(-2*alpha*(r - r_0)) - 2*exp(-alpha*(r - r_0)) ] - u_cut
!      In code: al = D_e/kT, bl = alpha, cl = r_0, e12 = exp(-bl*(r - cl)).
!
!   4. fpot_elecMorse:
!      Morse + Real-space screened Coulomb potential:
!        u(r) = u_Morse(r) + [ q_i * q_j / (4*pi*epsilon_0) ] * erfc(kappa*r) / r
!
! Arguments for all functions:
!   r2  (Input) : Squared Euclidean interatomic distance r_ij^2 (Angstrom^2).
!   nit (Input) : Pair interaction index from matrix itp(iti, itj).
!
! Returns:
!   Pair potential energy in reduced units (dimensionless, scaled by 1/kT).
!===============================================================================

!-------------------------------------------------------------------------------
! Function: fpot_LJ
! Purpose: Pure Lennard-Jones (12-6) pair potential.
!-------------------------------------------------------------------------------
function fpot_LJ(r2,nit)
  use set_precision
  Use potential, Only : al, bl, cl, bl2, kappa, qprod, ctr, elect, keyp, ucut
  implicit none
  Real(wp) :: fpot_LJ, rr, r2, fpel, r6
  integer :: nit

  r6 = (bl2(nit)/r2)**3
  fpot_LJ = al(nit)*r6*(r6-1.0d0) - ucut(nit)
end function fpot_LJ


!-------------------------------------------------------------------------------
! Function: fpot_elecLJ
! Purpose: Lennard-Jones (12-6) + real-space screened Coulomb (Ewald).
!-------------------------------------------------------------------------------
function fpot_elecLJ(r2,nit)
  use set_precision
  Use potential, Only : al, bl, cl, bl2, kappa, qprod, ctr, elect, keyp, ucut
  implicit none
  Real(wp) :: fpot_elecLJ, rr, r2, fpel, r6
  integer :: nit

  rr = Sqrt(r2)
  fpel = ctr*qprod(nit)*Erfc(kappa*rr)/rr
  r6 = (bl2(nit)/r2)**3
  fpot_elecLJ = al(nit)*r6*(r6-1.0d0) + fpel - ucut(nit)
end function fpot_elecLJ


!-------------------------------------------------------------------------------
! Function: fpot_Morse
! Purpose: Pure Morse pair potential.
!-------------------------------------------------------------------------------
function fpot_Morse(r2,nit)
  use set_precision
  Use potential, Only : al, bl, cl, bl2, kappa, qprod, ctr, elect, keyp, ucut
  implicit none
  Real(wp) :: fpot_Morse, rr, r2, fpel, e12
  integer :: nit

  rr = Sqrt(r2)
  e12 = Exp(-bl(nit)*(rr-cl(nit)))
  fpot_Morse = al(nit)*e12*(e12-2.0d0) - ucut(nit)
end function fpot_Morse


!-------------------------------------------------------------------------------
! Function: fpot_elecMorse
! Purpose: Morse + real-space screened Coulomb (Ewald).
!-------------------------------------------------------------------------------
function fpot_elecMorse(r2,nit)
  use set_precision
  Use potential, Only : al, bl, cl, bl2, kappa, qprod, ctr, elect, keyp, ucut
  implicit none
  Real(wp) :: fpot_elecMorse, rr, r2, fpel, r6, e12
  integer :: nit

  rr = Sqrt(r2)
  fpel = ctr*qprod(nit)*Erfc(kappa*rr)/rr
  e12 = Exp(-bl(nit)*(rr-cl(nit)))
  fpot_elecMorse = al(nit)*e12*(e12-2.0d0) + fpel - ucut(nit)
end function fpot_elecMorse

