Module Energy
  Use set_precision
  Use configuration
  Use potential
  Use properties
  Use rundata
Contains 
  Subroutine energ
    !
    !  Calculate energy
    !
    Implicit None
    Real(wp) :: rdd(ndim), rd2, rr, eng,  eng_f
    Integer :: i,j, iti, itj, nit
    Real(wp), external :: dist2, fpot
    !
    ! Calculate pairwise energies (Morse and real part of Ewald)
    ! Full calculation
    !
    eng = 0
    Do i = 1, natoms
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
    E_sr = eng
    !Write(*,'(80("-")/" vdW Energy=",g15.7,1x,a2)')eng,units
    eng_f = 0
    if (elect) then
       !
       !  Calculate Fourier component of Ewald sum
       !
       call fourier(eng_f)
       Etotal = E_sr+eng_f+selfe
       E_coulomb = eng_f+selfe
       E_Fourier = eng_f
       !write(*, '(" K cutoff",f10.5," 1/A")')rcpcut
!       write(*, '(" Coulomb energy (real space) ",g15.7,1x,a2)')  units
       !write(*, '(" Coulomb energy (self energy) ",g15.7,1x,a2)') selfe, units
       !write(*, '(" Coulomb energy (Fourier term)",g15.7,1x,a2)') eng_f, units
       !write(*, '(" Coulomb energy ",g15.7,1x,a2)') selfe+eng_f,units
       !write(*, '(" E_sr =",g15.7,1x,a2)') E_sr, units
       !write(*, '(" Deviation from charge neutrality =",g15.7)') qtotal
    endif
    !Write(*, '(" Etotal =",g15.7,1x,a2/80("-"))') selfe+eng_f+eng,units

  End Subroutine energ

  Subroutine energ_cell
    use linkcell
    use interp
    !
    !  Calculate energy using linked cells. This version uses cubic
    !  interpolation for short range interactions, hence the short range part
    !  of the Coulomb potential cannot be separated.
    Implicit None
    Real(wp) :: rdd(ndim), rd2, rr, eng,  eng_f
    Real (wp) :: y0, y1, y2, y3, a0, a1, a2, a3, mu, mu2, xmu
    Integer :: i,j, iti, itj, ix, jy, kz, cell, icell, nit, ir, nop=0
    Real(wp), external :: dist2, fpot
    !
    ! Calculate pairwise energies 
    !
    eng = 0
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
!!$                Where (rdd(:) > side2(:)) rdd(:) = rdd(:) -side(:)
!!$                Where (rdd(:) < -side2(:)) rdd(:) = rdd(:) + side(:)
!!$                rd2 = Dot_product(rdd(:),rdd(:))
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
                nop=nop+1
             end if
             j = list(j)
          End Do
       End Do
    end Do
    eng = eng/2
    E_sr = eng
    !Write(*,'(80("-")/" vdW Energy=",g15.7,1x,a2)')eng,units
    eng_f = 0
    if (elect) then
       !
       !  Calculate Fourier component of Ewald sum
       !
       call fourier(eng_f)
       Etotal = E_sr+eng_f+selfe
       E_coulomb = eng_f+selfe
       E_Fourier = eng_f
       !write(*, '(" K cutoff",f10.5," 1/A")')rcpcut
       !write(*, '(" Coulomb energy (self energy) ",g15.7,1x,a2)') selfe, units
       !write(*, '(" Coulomb energy (Fourier term)",g15.7,1x,a2)') eng_f, units
       !write(*, '(" Coulomb energy ",g15.7,1x,a2)') selfe+eng_f,units
       !write(*, '(" E_sr =",g15.7,1x,a2)') E_sr, units
       !write(*, '(" Deviation from charge neutrality =",g15.7)') qtotal
    endif
    !Write(*, '(" Etotal =",g15.7,1x,a2/80("-"))') selfe+eng_f+eng,units

  End Subroutine energ_cell


  Subroutine fourier(eng_f)
    !
    !  Routine to determine Fourier component of Ewald summation
    !
    use configuration
    use potential
    Implicit None
    Real(wp), Intent(out) :: eng_f
    Real(wp) :: ax, bx, cx, eng_f1, eng_f2
    Complex (wp) :: rhok
    Integer :: i, k, kx, ky, kz
    ! Use only on orthogonal cells !!!!
    ! Initialize matrices containing Exp(i*2pi*(kx/Lx+ky/Ly+kz/Lz))
    ! Unreduce atomic coordinates

    ! eix, eiy, eiz only change with volume for rigid molecules  
    ! or flexible if we do not scale the molecule with the volume
    ! otherwise they remain constant with volume changes, in this case they only need to be computed once
    eix(1:natoms,1) = Exp(ii*dospix*r(1:natoms,1)*r_unit(1))   
    eiy(1:natoms,1) = Exp(ii*dospiy*r(1:natoms,2)*r_unit(2))
    eiz(1:natoms,1) = Exp(ii*dospiz*r(1:natoms,3)*r_unit(3))
    eix(1:natoms,-1) =  Conjg(eix(1:natoms,1))
    eiy(1:natoms,-1) =  Conjg(eiy(1:natoms,1))
    !
    ! Calculate using iterative approach
    !
    Do kx=2,kmx
       eix(1:natoms,kx)=eix(1:natoms,kx-1)*eix(1:natoms,1)
       eix(1:natoms,-kx) =  Conjg(eix(1:natoms,kx))
       eiy(1:natoms,kx)=eiy(1:natoms,kx-1)*eiy(1:natoms,1)
       eiy(1:natoms,-kx) =  Conjg(eiy(1:natoms,kx))
       eiz(1:natoms,kx)=eiz(1:natoms,kx-1)*eiz(1:natoms,1)
    End Do
    Do ky=kmx+1,kmy
       eiy(1:natoms,ky)=eiy(1:natoms,ky-1)*eiy(1:natoms,1)
       eiy(1:natoms,-ky) =  Conjg(eiy(1:natoms,ky))
       eiz(1:natoms,ky)=eiz(1:natoms,ky-1)*eiz(1:natoms,1)
    End Do
    Do kz=kmy+1,kmz
       eiz(1:natoms,kz)=eiz(1:natoms,kz-1)*eiz(1:natoms,1)
    End Do
    !
    ! End initializations
    !

    eng_f=0.0d0
    k=1
    Do kx = 1, kmx
       If (km2(k) .Le. rcpcut2) Then
          rhok = Sum(q(:)*eix(:,kx))
          rhokk(k) = rhok
          eng_f = eng_f + ekm2(k)*rhok*Conjg(rhok)
       End If
       k=k+1
    End Do
    Do ky = 1, kmy
       Do kx = -kmx, kmx
          If (km2(k) .Le. rcpcut2) Then
             rhok = Sum(q(:)*eix(:,kx)*eiy(:,ky))
             rhokk(k) = rhok
             eng_f = eng_f + ekm2(k)*rhok*Conjg(rhok)
          End If
          k=k+1
       End Do
    End Do
    Do kz = 1, kmz
       Do ky = -kmy, kmy
          Do kx = -kmx, kmx
             If (km2(k) .Le. rcpcut2) Then
                rhok = Sum(q(:)*eix(:,kx)*eiy(:,ky)*eiz(:,kz))
                rhokk(k) = rhok
                eng_f = eng_f + ekm2(k)*rhok*Conjg(rhok)
             End If
             k=k+1
          End Do
       End Do
    End Do
    !
    ! A factor 2 appears due to the symmetry simplifications
    !
    eng_f = 2*eng_f*pi2*ctr/v0
  End Subroutine fourier


End Module Energy
