function fpot_LJ(r2,nit)
  use set_precision
  Use potential, Only : al, bl, cl, bl2, kappa, qprod, ctr, elect, keyp, ucut
  implicit none
  !
  !
  Real(wp) :: fpot_LJ, rr, r2, fpel, r6
  integer :: nit
  !
  !
  ! Lennard-Jones interaction
  r6 = (bl2(nit)/r2)**3
  fpot_LJ = al(nit)*r6*(r6-1.0d0)-ucut(nit)

end function fpot_LJ

function fpot_elecLJ(r2,nit)
  use set_precision
  Use potential, Only : al, bl, cl, bl2, kappa, qprod, ctr, elect, keyp, ucut
  implicit none
  !
  !
  Real(wp) :: fpot_elecLJ, rr, r2, fpel, r6
  integer :: nit
  !
  !
  ! Real space  electrostatics.
  rr = Sqrt(r2)
  fpel = ctr*qprod(nit)*Erfc(kappa*rr)/rr
  ! Lennard-Jones interaction+real space contribution of Ewald term
  r6 = (bl2(nit)/r2)**3
  fpot_elecLJ = al(nit)*r6*(r6-1.0d0)+fpel-ucut(nit)

end function fpot_elecLJ

function fpot_Morse(r2,nit)
  use set_precision
  Use potential, Only : al, bl, cl, bl2, kappa, qprod, ctr, elect, keyp, ucut
  implicit none
  !
  !
  Real(wp) :: fpot_Morse, rr, r2, fpel, e12
  integer :: nit
  !
  !
  rr = Sqrt(r2)
  ! Morse potential and real space contribution of Ewald term
  e12 = Exp(-bl(nit)*(rr-cl(nit)))
  ! Morse potential and real space contribution of Ewald term
  fpot_Morse = al(nit)*e12*(e12-2)-ucut(nit)

end function fpot_Morse

function fpot_elecMorse(r2,nit)
  use set_precision
  Use potential, Only : al, bl, cl, bl2, kappa, qprod, ctr, elect, keyp, ucut
  implicit none
  !
  !
  Real(wp) :: fpot_elecMorse, rr, r2, fpel, r6, e12
  integer :: nit
  !
  !
  rr = Sqrt(r2)
  fpel = ctr*qprod(nit)*Erfc(kappa*rr)/rr
  e12 = Exp(-bl(nit)*(rr-cl(nit)))
  ! Morse potential and real space contribution of Ewald term
  fpot_elecMorse = al(nit)*e12*(e12-2)+fpel-ucut(nit)
  
end function fpot_elecMorse

