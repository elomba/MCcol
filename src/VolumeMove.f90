Module VolumeChange

  Use Set_precision
  Use Interfaces
  Use Linkcell, Only : use_cell
  Use Thermo
  Use Energy, Only : Energ, Energ_cell
  Use Init, Only : Init_fourier

Contains 

 Subroutine move_volume
  implicit none

  if (scaling == 'cubic') then

     call move_volume_cubic 

  elseif (scaling == 'ortho'.or.scaling == 'isotr') then 

     call move_volume_ortho 

  Endif

 End subroutine move_volume


 Subroutine move_volume_cubic
!
!
!    Beware that to speed up the code the kappa is scaled together with the simulation box
!    This will cause problems if you want to evaluate the critical point
!    This part is not working properly: the dispersion and real parts should be separated
!    As kappa is also scaled, it is not possible to use the tabulated potential
!
  Use linkcell
  Use interp   
  implicit none
  Integer    ::  i, j, ix, jy, kz, icell, cell, iti,itj, nit, ir
  Real (wp)  ::  rd2, rr, xmu, mu, mu2, y0, y1, y2, y3, a0, a1, a2, a3, eng, rdd(ndim)
  Real (wp)  ::  lnvn, xi, deltaEt, v_new, factor, E_sr_new, E_fourier_new, selfe_new
  Real (wp)  ::  E_total_new, side_old, factor_inv

  !Call Random_number(xi)
  !lnvn = Log(v0)+(2*xi-1)*vdmax
  !v_new = Exp(lnvn)
  side_old = side(1)
  Call Random_number(xi)
  side(1) = side(1) + (2*xi-1)*vdmax
  factor = side(1) / side_old
  side(2) = side(2) * factor
  side(3) = side(3) * factor
  v_new = side(1)*side(2)*side(3)
  rc2 = rc2 * factor**2.


  If(use_cell) Then
    !
    eng = 0.
    Do i = 1, natoms
       !
       ! Locate cell
       !
       ix=int((r(i,1)+0.5d0)/cellx)
       jy=int((r(i,2)+0.5d0)/celly)
       kz=int((r(i,3)+0.5d0)/cellz)
       cell=(ix*maxj+jy)*maxk+kz
       !
       ! Loop over neighbouring cells
       !
       do icell = 1, nn

          j = head(neigh(cell,icell))
          do while (j .ne. 0)
             if (j.ne.i) then
                iti = iatype(i)
                itj = iatype(j)
                nit = itp(iti,itj)
                rdd(:) = r(i,:)-r(j,:)
                rd2 = dist2(rdd)
                If (rd2 < rc2(nit)) Then
                   rr =Sqrt(rd2)
!
! Use smooth cubic interpolation with coefficients as suggested by
! Paul Breeuwsma
!
                   xmu = rr/dr
                   ir = int(xmu)
                   mu = xmu-ir
                   mu2 = mu*mu
                   y0 = utab(ir-1,nit)
                   y1 = utab(ir,nit)
                   y2 = utab(ir+1,nit)
                   y3 = utab(ir+2,nit)
                   a0 = -0.5d0*y0 + 1.5d0*y1 - 1.5d0*y2 + 0.5d0*y3
                   a1 = y0 - 2.5d0*y1 + 2*y2 - 0.5d0*y3
                   a2 = -0.5d0*y0 + 0.5d0*y2
                   a3 = y1
                   eng = eng+a0*mu*mu2+a1*mu2+a2*mu+a3
                End If
             end if
             j = list(j)
          End Do
       End Do
    end Do
    eng = eng/2
    E_sr_new = eng

  Else

    eng = 0.
    Do i = 1, natoms-1
       Do j = i+1, natoms
          iti = iatype(i)
          itj = iatype(j)
          nit = itp(iti,itj)
          rdd(:) = r(i,:)-r(j,:)
          rd2 = dist2(rdd)
          If (rd2 < rc2(nit)) Then
             iti = iatype(i)
             itj = iatype(j)
             eng = eng + fpot(rd2,nit)
          End If
       End Do
    End Do
    E_sr_new = eng

  EndIf

  factor_inv =  1.d0 /factor
  E_Fourier_new = E_Fourier * factor_inv
  selfe_new = selfe * factor_inv
  E_total_new= E_sr_new + E_Fourier_new + selfe_new

  deltaEt =  E_total_new - (E_sr + E_Fourier + selfe)
  
  Call Random_number(xi)
  !If (Exp(-(deltaEt+pres*(v_new-v0))/kT+dble(natoms+1)*log(v_new/v0)) .Gt. xi) Then
  If (Exp(-(deltaEt+pres*(v_new-v0))/kT+dble(natoms)*log(v_new/v0)) .Gt. xi) Then
         nvaccept = nvaccept + 1
         E_sr = E_sr_new
         E_coulomb = E_coulomb * factor_inv
         E_fourier = E_fourier * factor_inv
         selfe = selfe * factor_inv
         Etotal = Etotal * factor_inv
         kappa = kappa * factor_inv
         rcut = rcut  * factor
         rcut2 = rcut**2.
         v0 = v_new
         a(1) = side(1)
         b(2) = side(2)
         c(3) = side(3)
         r_unit = side
         ekm2 = ekm2 *  factor **2.
         km2 = km2 * factor_inv**2.
         dospix = dospix * factor_inv
         dospiy = dospiy * factor_inv
         dospiz = dospiz * factor_inv
         rcpcut = 1.05*Min(dospix*kmx,dospiy*kmy,dospiz*kmz)
         rcpcut2 = rcpcut**2
  Else
         side(1) = side(1) * factor_inv
         side(2) = side(2) * factor_inv
         side(3) = side(3) * factor_inv
         rc2 = rc2 * factor_inv**2.
  Endif

 End subroutine move_volume_cubic


 Subroutine move_volume_ortho
!
!   If the ratio between the edges of the box is maintained
!     we could do some scaling, but we need to define kappa_x, kappa_y and kappa_z
!     that possibility is not yet implemented
!
  implicit none
  Real (wp)  ::  xi, deltaEt, deltaE_sr
  Real (wp)  ::  Etotal_old, E_sr_old, E_coulomb_old, E_fourier_old, v_old
  Real (wp)  , dimension(ndim) :: a_old, b_old, c_old, side_old, r_unit_old
  Real (wp)  ::  deltaE_Fourier, rcpcut2_old, factor, rcpcut_old
  Real (wp)  ::  dospix_old, dospiy_old, dospiz_old
  Complex(wp) :: rhokk_old(kmt), ekm2_old(kmt), km2_old(0:kmt)
  Complex(wp) :: eix_old(1:natoms,-kmx:kmx),eiy_old(1:natoms,-kmy:kmy),eiz_old(1:natoms,0:kmz)

  E_sr_old = E_sr
  side_old=side
  a_old=a
  b_old=b
  c_old=c
  r_unit_old=r_unit
  v_old=v0

  If (elect) Then
    Etotal_old = Etotal  !Etot and Ecoul are not updated in move
    E_coulomb_old = E_coulomb
    E_fourier_old = E_fourier
    ekm2_old = ekm2
    km2_old = km2
    rcpcut2_old = rcpcut2  ! just in case there are large volume changes
    rcpcut_old = rcpcut  ! just in case there are large volume changes
    rhokk_old = rhokk   ! 
    eix_old = eix   ! for rigid molecules, we need also to actualize eix,eiy,eiz
    eiy_old = eiy   ! or for flexible molecules that are not scale with volume
    eiz_old = eiz   ! for atomic systems it is actually not necessary
    dospix_old = dospix
    dospiy_old = dospiy
    dospiz_old = dospiz
  EndIf

  Call Random_number(xi)
  side(1) = side(1) + (2*xi-1)*vdmax

  If (scaling == 'ortho') Then  ! the edges change independently

     Call Random_number(xi)
     side(2) = side(2) + (2*xi-1)*vdmax
     Call Random_number(xi)
     side(3) = side(3) + (2*xi-1)*vdmax

  Else   ! the edges change by the same factor

     factor = side(1) / side_old (1)
     side(2) = side(2) *factor
     side(3) = side(3) *factor

  Endif

  v0=side(1)*side(2)*side(3)

  a(1) = side(1)
  b(2) = side(2)
  c(3) = side(3)
  r_unit(:) = side(:)

  If (elect) call Init_fourier
  if (use_cell) then    
     Call Energ_cell
  else
     Call Energ
  endif

  deltaE_sr = E_sr - E_sr_old
  If (elect) Then
        deltaE_Fourier = E_fourier - E_fourier_old
  Else 
        deltaE_Fourier = 0.d0
  Endif
  deltaEt =  deltaE_sr + deltaE_Fourier 

     Call Random_number(xi)
     If (Exp(-(deltaEt+pres*(v0-v_old))/kT+dble(natoms)*Log(v0/v_old)) .Gt. xi) Then
         nvaccept = nvaccept + 1
     Else   ! reject
         Etotal = Etotal_old  ! etot and ecoulomb are not updated in Move
         E_sr = E_sr_old
         a = a_old
         b = b_old
         c = c_old
         side = side_old
         r_unit = r_unit_old
         v0 = v_old
         If (elect) Then
           E_coulomb = E_coulomb_old
           E_fourier = E_fourier_old
           ekm2 = ekm2_old
           km2 = km2_old
           rcpcut2 = rcpcut2_old
           rcpcut = rcpcut_old
           rhokk = rhokk_old
           eix = eix_old   ! for rigid molecules, we also need to update eix
           eiy = eiy_old
           eiz = eiz_old
           dospix = dospix_old
           dospiy = dospiy_old
           dospiz = dospiz_old
         Endif
  Endif

 End subroutine move_volume_ortho
 
End Module VolumeChange
