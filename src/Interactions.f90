function fpot(r2,nit)
  use set_precision
  Use potential, Only : al, bl, cl, bl2, kappa, qprod, ctr, elect, keyp
  implicit none
  !
  !
  Real(wp) :: fpot, rr, r2, fpel, r6
  integer :: nit
  !
  !
  if (elect) then
     ! Compute only when electrostatics are needed.
     rr = Sqrt(r2)
     fpel = ctr*qprod(nit)*Erfc(kappa*rr)/rr
  else
     fpel = 0.0d0
  endif
  select case (keyp(nit))
  case (1)
     If (.Not. elect) rr = Sqrt(r2)
  ! Morse potential and real space contribution of Ewald term
     fpot = al(nit)*((1-Exp(-bl(nit)*(rr-cl(nit))))&
          &**2-1)+fpel
  case(2)
  ! Lennard-Jones interaction+real space contribution of Ewald term
     r6 = (bl2(nit)/r2)**3
     fpot = 4*al(nit)*r6*(r6-1.0d0)+fpel
  case default
     print *, " *** Error: interaction not implemented "
  end select
     
end function fpot

