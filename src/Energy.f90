!===============================================================================
! Module: Energy
!
! Purpose:
!   Computes total and partial potential energies of the atomistic system:
!   - Short-range pair interactions (Lennard-Jones or Morse)
!   - Real-space screened electrostatic Ewald sum
!   - Reciprocal-space (Fourier) Ewald sum
!   - Electrostatic self-energy correction
!
! Evaluation methods:
!   1. energ          : Direct all-pairs O(N^2) evaluation.
!   2. energ_cell     : Link-cell O(N) evaluation with analytical pair functions.
!   3. energ_cell_int : Link-cell O(N) evaluation with cubic spline interpolation.
!   4. fourier        : Reciprocal space Ewald sum using complex exponentials.
!===============================================================================
Module Energy
  Use set_precision
  Use configuration
  Use potential
  Use properties
  Use rundata
  Use interfaces, Only : dist2
Contains 

  !-----------------------------------------------------------------------------
  ! Subroutine: energ
  !
  ! Purpose:
  !   Computes total energy by direct all-pairs O(N^2) double summation.
  !   Used for baseline initialization and when link cells cannot be formed.
  !
  ! Arguments:
  !   f (external function) : Pair potential evaluation function.
  !-----------------------------------------------------------------------------
  Subroutine energ(f)
    Implicit None
    Real(dkind) :: eng, eng_f
    Real(wp) :: rdd(ndim), rd2, rr
    Integer :: i,j, iti, itj, nit
    Real(wp), external :: f

    ! Direct all-pairs loop over unique pairs (i < j)
    eng = 0
    Do i = 1, natoms
       Do j = i+1, natoms
          iti = iatype(i)
          itj = iatype(j)
          nit = itp(iti,itj)
          rdd(:) = R(:,i)-R(:,j)
          rd2 = dist2(rdd)
          If (rd2 < rc2(nit)) Then
             eng = eng + f(rd2,nit)
          End If
       End Do
    End Do
    E_sr = eng

    eng_f = 0.0d0
    if (elect) then
       ! Calculate Fourier reciprocal space component
       if (fou_type == 1) then
          call fourier(eng_f)
       else
          call four_pme(eng_f)
       endif
       Etotal = E_sr + eng_f + selfe
       E_coulomb = eng_f + selfe
       E_Fourier = eng_f
    endif
  End Subroutine energ

  !-----------------------------------------------------------------------------
  ! Subroutine: energ_cell
  !
  ! Purpose:
  !   Computes short-range pair energy using the 3D link-cell method (O(N)),
  !   calling the analytical pair potential function f.
  !
  ! Arguments:
  !   f (external function) : Pair potential evaluation function.
  !-----------------------------------------------------------------------------
  Subroutine energ_cell(f)
    use linkcell
    Implicit None
    real(dkind) :: eng, eng_f
    Real(wp), external :: f

    call Eshort_r(f, eng)
    E_sr = eng

    eng_f = 0.0d0
    if (elect) then
       if (fou_type == 1) then
          call fourier(eng_f)
       else
          call four_pme(eng_f)
       endif
       Etotal = E_sr + eng_f + selfe
       E_Fourier = eng_f
    endif
  End Subroutine energ_cell

  !-----------------------------------------------------------------------------
  ! Subroutine: energ_cell_int
  !
  ! Purpose:
  !   Computes short-range energy using link cells coupled with Paul Breeuwsma
  !   cubic spline interpolation from pre-tabulated potential tables.
  !-----------------------------------------------------------------------------
  Subroutine energ_cell_int
    use linkcell
    use interp
    !
    !  Calculate energy using linked cells. This version uses cubic
    !  interpolation for short range interactions, hence the short range part
    !  of the Coulomb potential cannot be separated.
    Implicit None
    real(dkind) ::  eng,  eng_f
    Real(wp) :: rdd(ndim), rd2, rr
    Real (wp) :: y(0:3), a(0:3),  mu, mu2, xmu
    Integer :: i,j, iti, itj, ix, jy, kz, cell, icell, nit, ir, nop=0
    Real(wp), external :: fpot
    !
    ! Calculate pairwise energies 
    !
    eng = 0
    Do i = 1, natoms
       !
       ! Locate cell
       !
       ix=int((R(1,i)+0.5d0)/cellx)
       jy=int((R(2,i)+0.5d0)/celly)
       kz=int((R(3,i)+0.5d0)/cellz)
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
                rdd(:) = R(:,i)-R(:,j)
                rd2 = dist2(rdd)
                If (rd2 < rc2(nit)) Then
                   rr =Sqrt(rd2)
!
! Use smooth cubic interpolation with coefficients as suggested by
! Paul Breeuwsma
!
                   xmu = rr*idr
                   ir = int(xmu)
                   mu = xmu-ir
                   mu2 = mu*mu
                   y(0:3) = utab(ir-1:ir+2,nit)
                   a= matmul(am,y)
                   eng = eng+(a(0)*mu+a(1))*mu2+a(2)*mu+a(3)
                End If
                nop=nop+1
             end if
             j = list(j)
          End Do
       End Do
    end Do
    eng = eng/2
    E_sr = eng
    !Write(*,'(80("-")/" **vdW Energy=",g15.7,1x,a2)')eng,units
    eng_f = 0
    if (elect) then
       !
       !  Calculate Fourier component of Ewald sum
       !
       if (fou_type == 1) then
          call fourier(eng_f)
       else
          call four_pme(eng_f)
       endif
       Etotal = E_sr+eng_f+selfe
       E_coulomb = eng_f+selfe
       E_Fourier = eng_f
       !write(*, '(" K cutoff",f10.5," 1/A")')rcpcut
       !write(*, '(" Coulomb energy (self energy) ",g15.7,1x,a2)') selfe, units
       !write(*, '(" **Coulomb energy (Fourier term)",g15.7,1x,a2)') eng_f, units
       !write(*, '(" Coulomb energy ",g15.7,1x,a2)') selfe+eng_f,units
       !write(*, '(" E_sr =",g15.7,1x,a2)') E_sr, units
       !write(*, '(" Deviation from charge neutrality =",g15.7)') qtotal
    endif
    !Write(*, '(" Etotal =",g15.7,1x,a2/80("-"))') selfe+eng_f+eng,units

  End Subroutine energ_cell_int


  !-----------------------------------------------------------------------------
  ! Subroutine: fourier
  !
  ! Purpose:
  !   Computes the reciprocal-space (Fourier) component of the Ewald summation:
  !     E_Four = (1 / 2*epsilon_0*V) * sum_{k /= 0} [ exp(-k^2 / 4*kappa^2) / k^2 ] * |rho(k)|^2
  !   where rho(k) = sum_i q_i * exp(i * k * r_i).
  !
  ! Optimization:
  !   Takes advantage of inversion symmetry rho(-k) = rho(k)*. The sum is
  !   restricted to half-space (kz > 0; or kz=0, ky > 0; or kz=ky=0, kx > 0)
  !   and multiplied by a factor of 2.
  !   Complex exponentials are updated recursively:
  !     exp(i * kx * x) = exp(i * (kx-1) * x) * exp(i * x).
  !
  ! Arguments:
  !   eng_f (out) : Reciprocal space electrostatic energy.
  !-----------------------------------------------------------------------------
  Subroutine fourier(eng_f)
    use configuration
    use potential
    Implicit None
    Real(dkind), Intent(out) :: eng_f
    Real(wp) :: ax, bx, cx
    Complex (wp) :: rhok
    Integer :: i, k, kx, ky, kz

    ! Initialize base complex exponentials for k = 1 along each axis
    eix(1:natoms,1) = Exp(ii*dospix*R(1,1:natoms)*r_unit(1))   
    eiy(1:natoms,1) = Exp(ii*dospiy*R(2,1:natoms)*r_unit(2))
    eiz(1:natoms,1) = Exp(ii*dospiz*R(3,1:natoms)*r_unit(3))
    eix(1:natoms,-1) = Conjg(eix(1:natoms,1))
    eiy(1:natoms,-1) = Conjg(eiy(1:natoms,1))

    ! Recursively generate higher harmonic phase factors: exp(i * k * r)
    Do kx=2,kmx
       eix(1:natoms,kx) = eix(1:natoms,kx-1)*eix(1:natoms,1)
       eix(1:natoms,-kx) = Conjg(eix(1:natoms,kx))
       eiy(1:natoms,kx) = eiy(1:natoms,kx-1)*eiy(1:natoms,1)
       eiy(1:natoms,-kx) = Conjg(eiy(1:natoms,kx))
       eiz(1:natoms,kx) = eiz(1:natoms,kx-1)*eiz(1:natoms,1)
    End Do
    Do ky=kmx+1,kmy
       eiy(1:natoms,ky) = eiy(1:natoms,ky-1)*eiy(1:natoms,1)
       eiy(1:natoms,-ky) = Conjg(eiy(1:natoms,ky))
       eiz(1:natoms,ky) = eiz(1:natoms,ky-1)*eiz(1:natoms,1)
    End Do
    Do kz=kmy+1,kmz
       eiz(1:natoms,kz) = eiz(1:natoms,kz-1)*eiz(1:natoms,1)
    End Do

    ! Sum reciprocal energy over unique symmetry-reduced k-vectors
    eng_f = 0.0d0
    k = 1
    ! Sector 1: kz = 0, ky = 0, kx > 0
    Do kx = 1, kmx
       If (km2(k) .Le. rcpcut2) Then
          rhok = Sum(q(:)*eix(:,kx))
          rhokk(k) = rhok
          eng_f = eng_f + ekm2(k)*rhok*Conjg(rhok)
       End If
       k=k+1
    End Do
    ! Sector 2: kz = 0, ky > 0, all kx
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
    ! Sector 3: kz > 0, all ky, all kx
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

    ! Multiply by 2 for inversion symmetry [-k, k] and volume prefactor
    eng_f = 2.0d0 * eng_f * pi2 * ctr / v0
  End Subroutine fourier

  !-----------------------------------------------------------------------------
  ! Subroutine: four_pme
  !
  ! Purpose:
  !   Placeholder interface for Particle Mesh Ewald (PME) reciprocal sum.
  !-----------------------------------------------------------------------------
  Subroutine four_pme(eng_f)
    use configuration
    use potential
    Implicit None
    Real(dkind), Intent(out) :: eng_f
    eng_f = 0.0d0
    eng_f = 2.0d0 * eng_f * pi2 * ctr / v0
  End Subroutine four_pme

  !-----------------------------------------------------------------------------
  ! Subroutine: Eshort_r
  !
  ! Purpose:
  !   Computes total short-range energy using the 3D link-cell linked lists,
  !   looping over atoms and their 27 neighboring cells.
  !
  ! Arguments:
  !   f   (external function) : Pair interaction energy function.
  !   eng (inout)             : Output accumulated energy (divided by 2 for double count).
  !-----------------------------------------------------------------------------
  Subroutine Eshort_r(f,eng)
    use linkcell
    implicit none
    Integer :: i,j, iti, itj, ix, jy, kz, cell, icell, nit
    real(dkind), intent(INOUT) ::  eng
    Real(wp) :: rdd(ndim), rd2, rr
    Real(wp), external :: f

    eng = 0.0d0
    Do i = 1, natoms
       ix = int((R(1,i)+0.5d0)/cellx)
       jy = int((R(2,i)+0.5d0)/celly)
       kz = int((R(3,i)+0.5d0)/cellz)
       cell = (ix*maxj+jy)*maxk+kz
       ! Loop over 27 neighbouring cells
       do icell = 1, nn
          j = head(neigh(cell,icell))
          do while (j .ne. 0)
             if (j .ne. i) then
                iti = iatype(i)
                itj = iatype(j)
                nit = itp(iti,itj)
                rdd(:) = R(:,i)-R(:,j)
                rd2 = dist2(rdd)
                If (rd2 < rc2(nit)) Then
                   eng = eng + f(rd2,nit)
                End If
             end if
             j = list(j)
          End Do
       End Do
    end Do
    ! Divide by 2 to prevent double counting
    eng = eng / 2.0d0
  end Subroutine Eshort_r

End Module Energy

