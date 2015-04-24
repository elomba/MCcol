
Module Moving 
  Use set_precision
  Use Thermo
  Use rundata
  Use properties
Contains 
  Subroutine moven
    !
    ! Brute force move routine for displacement moves.
    !
    Implicit None
    Integer :: ntest, i, iti, itj, nit
    Real (wp), External :: fpot
    Real (wp) :: rdd(ndim), rddn(ndim), rp(ndim), harvest(0:ndim), eng0, eng1
    Real (wp) :: rr, rd2, rdn2, deltaE, deltaFour, deltaEt, xi
    Real (wp), external :: dist2
    ntrial = ntrial+1
    ! Choose randomly particle to move
    Call Random_number(harvest(0:ndim))
    ntest = Int(natoms*harvest(0))+1
    !
    !
    !
    rp(1:ndim) = r(ntest,1:ndim)+rdmax(1:ndim)*(2*harvest(1:ndim)-1)/r_unit(1:ndim)
    !
    ! Recalculate trial movement according to periodic boundary conditions
    where (rp(:) < -0.5d0) rp(:)=rp(:)+1
    where (rp(:) > 0.5d0 ) rp(:)=rp(:)-1
    !
    ! Calculate test particle energy in the old and new configurations.
    !
    eng0 = 0.0d0
    eng1 = 0.0d0
    i=0
    do while (i<natoms)
       i=i+1
       if (i .eq. ntest) cycle
       rdd(:) = r(i,:)-r(ntest,:)
       rddn(:) = r(i,:)-rp(:)
       rd2 = dist2(rdd)
       iti = iatype(i)
       itj = iatype(ntest)
       nit = itp(iti,itj)
       If (rd2 < rc2(nit)) Then
          eng0 = eng0 + fpot(rd2,nit)
       End If
       rdn2 = dist2(rddn)
       If (rdn2 < rc2(nit)) Then
          eng1 = eng1 + fpot(rdn2,nit)
       End If
    End Do
    deltaE = eng1-eng0
    if (elect) then
       Call change_fourier(ntest,rp,deltaFour)
    else
       deltaFour = 0.0d0
    endif
    deltaEt = deltaE+deltaFour
    !
    ! Check acceptance
    !
    If (deltaEt < 0) Then
       r(ntest,:) = rp(:)
       E_sr = E_sr+deltaE
       if (elect) then
          E_Fourier = E_Fourier+deltaFour
          rhokk(1:kmt) = rhokk(1:kmt)+deltann(1:kmt)
          eix(ntest,:) = einx(:)
          eiy(ntest,:) = einy(:)
          eiz(ntest,:) = einz(:)
       endif
       naccept = naccept+1
    Else 
       Call Random_number(xi)
       If (Exp(-deltaEt/kT) .Gt. xi) Then
          r(ntest,:) = rp(:)
          E_sr = E_sr+deltaE
          if (elect) then
             E_Fourier = E_Fourier+deltaFour
             rhokk(1:kmt) = rhokk(1:kmt)+deltann(1:kmt)
             eix(ntest,:) = einx(:)
             eiy(ntest,:) = einy(:)
             eiz(ntest,:) = einz(:)
          endif
          naccept = naccept+1
       Endif
    End If
  End Subroutine moven

  Subroutine move_linkcell(f)
    !
    ! Move routine implementing link cell method
    !
    use linkcell
    use interp
    Implicit None
    Integer :: ntest, i, j, k, icell, iti, itj, cell, nit
    real (wp), external :: f
    Real (wp) :: rdd(ndim), rddn(ndim), rp(ndim), harvest(0:ndim), eng0, eng1
    Real (wp) :: rr, rd2, rdn2, deltaE, deltaFour, deltaEt, xi
    Real (wp), external :: dist2
    ntrial = ntrial+1
    ! Choose randomly particle to move
    Call Random_number(harvest(0:ndim))
    ntest = Int(natoms*harvest(0))+1
    !
    !
    !
    rp(1:ndim) = r(ntest,1:ndim)+rdmax(1:ndim)*(2*harvest(1:ndim)-1)/r_unit(1:ndim)
    !
    ! Recalculate trial movement according to periodic boundary conditions
    where (rp(:) < -0.5d0) rp(:)=rp(:)+1
    where (rp(:) > 0.5d0) rp(:)=rp(:)-1
    !
    ! Locate cell
    !
    i=int((rp(1)+0.5d0)/cellx)
    j=int((rp(2)+0.5d0)/celly)
    k=int((rp(3)+0.5d0)/cellz)
    cell=(i*maxj+j)*maxk+k

    !
    ! Calculate test particle energy in the old and new configurations
    !
    eng0 = 0.0d0
    eng1 = 0.0d0
    !
    ! Loop over neighbouring cells
    !
    Do icell = 1, nn
       i = head(neigh(cell,icell))
       do while (i .ne. 0)
          if (i.ne.ntest) then
             iti = iatype(i)
             itj = iatype(ntest)
             nit = itp(iti,itj)
             rdd(:) = r(i,:)-r(ntest,:)
             rddn(:) = r(i,:)-rp(:)
             rd2 = dist2(rdd)
             If (rd2 < rmin2(nit)) Return
             If (rd2 < rc2(nit)) Then
                eng0 = eng0 + fpot(rd2,nit)
             End If
             rdn2 = dist2(rddn)
             If (rdn2 < rc2(nit)) Then
                eng1 = eng1 + fpot(rdn2,nit)
             End If
          end if
          i=list(i)
       End Do
    end Do
    deltaE = eng1-eng0
    if (elect) then
       Call change_fourier(ntest,rp,deltaFour)
    else
       deltaFour = 0.0d0
    endif
    deltaEt = deltaE+deltaFour
    !
    ! Check acceptance
    !
    If (deltaEt < 0) Then
       r(ntest,:) = rp(:)
       E_sr = E_sr+deltaE
       if (elect) then
          E_Fourier = E_Fourier+deltaFour
          rhokk(1:kmt) = rhokk(1:kmt)+deltann(1:kmt)
          eix(ntest,:) = einx(:)
          eiy(ntest,:) = einy(:)
          eiz(ntest,:) = einz(:)
       endif
       naccept = naccept+1
    Else 
       Call Random_number(xi)
       If (Exp(-deltaEt/kT) .Gt. xi) Then
          r(ntest,:) = rp(:)
          E_sr = E_sr+deltaE
          if (elect) then
             E_Fourier = E_Fourier+deltaFour
             rhokk(1:kmt) = rhokk(1:kmt)+deltann(1:kmt)
             eix(ntest,:) = einx(:)
             eiy(ntest,:) = einy(:)
             eiz(ntest,:) = einz(:)
          endif
          naccept = naccept+1
       Endif
    End If
  End Subroutine move_linkcell

  Subroutine move_link_int
    !
    ! Move routine with link cell method and potential interpolation
    !
    use linkcell
    use interp
    use cells, only : update_cell_list, build_cells
    Implicit None
    Integer :: ntest, i, j, k, icell, iti, itj, cell, nit, ir
    Real (wp) :: rdd(ndim), rddn(ndim), rp(ndim), harvest(0:ndim), eng0, eng1
    Real (wp) :: y0, y1, y2, y3, a0, a1, a2, a3
    Real (wp) :: rr, rd2, rdn2, deltaE, deltaFour, deltaEt, xi, mu, mu2, xmu
    Real (wp), external :: dist2
    ntrial = ntrial+1
    ! Choose randomly particle to move
    Call Random_number(harvest(0:ndim))
    ntest = Int(natoms*harvest(0))+1
    !
    ! rescale displacement to box length units
    !
    rp(1:ndim) = r(ntest,1:ndim)+rdmax(1:ndim)*(2*harvest(1:ndim)-1)/r_unit(1:ndim)
    !
    ! Recalculate trial movement according to periodic boundary conditions
    where (rp(:) < -0.5d0) rp(:)=rp(:)+1
    where (rp(:) > 0.5d0) rp(:)=rp(:)-1
    !
    ! Locate cell
    !
    i=int((rp(1)+0.5d0)/cellx)
    j=int((rp(2)+0.5d0)/celly)
    k=int((rp(3)+0.5d0)/cellz)
    cell=(i*maxj+j)*maxk+k

    !
    ! Calculate test particle energy in the old and new configurations
    !
    eng0 = 0.0d0
    eng1 = 0.0d0
    !
    ! Loop over neighbouring cells
    !
    Do icell = 1, nn
       i = head(neigh(cell,icell))
!       write(18,*) i,j
       do while (i .ne. 0)
          if (i.ne.ntest) then
             iti = iatype(i)
             itj = iatype(ntest)
             nit = itp(iti,itj)
             rdd(:) = r(i,:)-r(ntest,:)
             rddn(:) = r(i,:)-rp(:)
!!$             Where (rdd(:) > side2(:)) rdd(:) = rdd(:) -side(:)
!!$             Where (rdd(:) < -side2(:)) rdd(:) = rdd(:) + side(:)
!!$             rd2 = Dot_product(rdd(:),rdd(:))
             rd2 = dist2(rdd)
             if (rd2 < rmin2(nit)) return
             If (rd2 < rc2(nit)) Then
                rr =Sqrt(rd2)
!
! Use smooth cubic interpolation with coefficients as suggested by
! Paul Breeuwsma (http://paulbourke.net/miscellaneous/interpolation/)
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
                eng0 = eng0+a0*mu*mu2+a1*mu2+a2*mu+a3
             End If
!!$             Where (rddn(:) > side2(:)) rddn(:) = rddn(:) -side(:)
!!$             Where (rddn(:) < -side2(:)) rddn(:) = rddn(:) + side(:)
!!$             rdn2 = Dot_product(rddn(:),rddn(:))
             rdn2 = dist2(rddn)
             If (rdn2 < rc2(nit)) Then
                rr =Sqrt(rdn2)
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
                eng1 = eng1+a0*mu*mu2+a1*mu2+a2*mu+a3
             End If
          end if
          i=list(i)
       End Do
    end Do
    deltaE = eng1-eng0
    if (elect) then
       Call change_fourier(ntest,rp,deltaFour)
    else
       deltaFour = 0.0d0
    endif
    deltaEt = deltaE+deltaFour
    !
    ! Check acceptance
    !
    If (deltaEt < 0) Then
       ! Check for particle leaving the cell and update cell list
       call update_cell_list(ntest,cell)
       r(ntest,:) = rp(:)
       E_sr = E_sr+deltaE
       if (elect) then
          E_Fourier = E_Fourier+deltaFour
          rhokk(1:kmt) = rhokk(1:kmt)+deltann(1:kmt)
          eix(ntest,:) = einx(:)
          eiy(ntest,:) = einy(:)
          eiz(ntest,:) = einz(:)
       endif
       naccept = naccept+1
    Else 
       Call Random_number(xi)
       If (Exp(-deltaEt/kT) .Gt. xi) Then
          ! Check for particle leaving cell and update cell list
          call update_cell_list(ntest,cell)
          r(ntest,:) = rp(:)
          E_sr = E_sr+deltaE
          if (elect) then
             E_Fourier = E_Fourier+deltaFour
             rhokk(1:kmt) = rhokk(1:kmt)+deltann(1:kmt)
             eix(ntest,:) = einx(:)
             eiy(ntest,:) = einy(:)
             eiz(ntest,:) = einz(:)
             naccept = naccept+1
          endif
       Endif
    End If
  !     call build_cells

  End Subroutine move_link_int

  Subroutine change_fourier(n,rtry,deltaf)
    !
    ! Compute change of Fourier component of potential energy when
    ! moving one particle 
    !
    Implicit None
    Integer, Intent(in) :: n
    Real(wp), Intent(in) :: rtry(ndim)
    Real(wp), Intent(out) :: deltaf
    Integer :: kx, ky, kz, k
    Real (wp) :: rold(ndim)
    Complex (wp) :: deltannp, deltannpm
    rold(1:ndim) = r(n,1:ndim)
    !
    ! Initialize matrices containing Exp(i*2pi*(kx/Lx+ky/Ly+kz/Lz)) for test particle
    ! Unreduce test position coordinates

    einx(1) = Exp(ii*dospix*rtry(1)*r_unit(1))
    einy(1) = Exp(ii*dospiy*rtry(2)*r_unit(2))
    einz(1) = Exp(ii*dospiz*rtry(3)*r_unit(3))
    einx(-1) =  Conjg(einx(1))
    einy(-1) =  Conjg(einy(1))
    !
    ! Calculate using iterative approach
    !
    Do kx=2,kmx
       einx(kx)=einx(kx-1)*einx(1)
       einx(-kx) =  Conjg(einx(kx))
       einy(kx)=einy(kx-1)*einy(1)
       einy(-kx) =  Conjg(einy(kx))
       einz(kx)=einz(kx-1)*einz(1)
    End Do
    Do ky=kmx+1,kmy
       einy(ky)=einy(ky-1)*einy(1)
       einy(-ky) =  Conjg(einy(ky))
       einz(ky)=einz(ky-1)*einz(1)
    End Do
    Do kz=kmy+1,kmz
       einz(kz)=einz(kz-1)*einz(1)
    End Do
    !
    ! End initializations
    !
    deltaf=0.0d0
    k=1
    Do kx = 1, kmx
       If (km2(k) .Le. rcpcut2) Then
          deltannp = q(n)*(einx(kx)-eix(n,kx))
          deltann(k) = deltannp
          deltannpm = Conjg(deltannp)
          deltaf = deltaf + ekm2(k)*(2*Real(rhokk(k)*deltannpm)&
               &+deltannp*deltannpm)
       End If
       k=k+1
    End Do
    Do ky = 1, kmy
       Do kx = -kmx, kmx
          If (km2(k) .Le. rcpcut2) Then
             deltannp = q(n)*(einx(kx)*einy(ky)-eix(n,kx)*eiy(n,ky))
             deltann(k) = deltannp
             deltannpm = Conjg(deltannp)
             deltaf = deltaf + ekm2(k)*(2*Real(rhokk(k)*deltannpm)+deltannp*deltannpm)
         End If
          k=k+1
       End Do
    End Do
    Do kz = 1, kmz
       Do ky = -kmy, kmy
          Do kx = -kmx, kmx
             If (km2(k) .Le. rcpcut2) Then
                deltannp = q(n)*(einx(kx)*einy(ky)*einz(kz)-eix(n,kx)*eiy(n,ky)*eiz(n,kz))
                deltann(k) = deltannp
                deltannpm = Conjg(deltannp)
                deltaf = deltaf + ekm2(k)*(2*Real(rhokk(k)*deltannpm)+deltannp*deltannpm)
             End If
             k=k+1
          End Do
       End Do
    End Do
    !
    ! A factor 2 appears due to the symmetry simplifications
    !
    deltaf = 2*deltaf*pi2*ctr/v0
  End Subroutine change_fourier

End Module Moving



Subroutine move_natoms(natoms)
  !
  ! Routine to move natoms sequentially 
  !
  Use set_precision
  use Moving, only : move_link_int, moven
  use linkcell, only : use_cell
  implicit none
  integer, intent(IN) :: natoms
  integer :: i
  do i=1,natoms
     if (use_cell) then
        !
        ! Use link cells and interpolations in the energy
        !
        call move_link_int
     else
        !
        ! Simple move with explicit evaluations
        ! 
        call moven
     end if
  end do
end subroutine move_natoms
