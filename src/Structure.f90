!===============================================================================
! Subroutine: structure
!
! Purpose:
!   Calculates and normalizes the partial Radial Distribution Functions (RDF),
!   g_ij(r), for multi-component systems under periodic boundary conditions.
!
! Mathematical Definition:
!   The partial pair correlation function between species j and l is:
!     g_jl(r) = (V / N_j * N_l) * < sum_{i in j} sum_{m in l, m /= i} delta(r - r_im) >
!   Normalized by the spherical shell volume at distance r:
!     Delta_V(r) = (4/3) * pi * [ (r + dr/2)^3 - (r - dr/2)^3 ]
!   and the number density rho = N / V and mole fraction x_l = N_l / N.
!
! Outputs:
!   gmix(bin, j, l) : Normalized pair correlation function stored in properties module.
!===============================================================================
Subroutine structure
  Use set_precision
  Use configuration, Only : r, iatype, ndim, nsp, ntype, natoms, v0
  Use properties, Only : nmaxgr, gmix, histomix, deltagr
  Use potential, Only : rcut2, pi
  Use interfaces, only : dist2
  Implicit None
  Real(wp) :: rd(ndim), rr2, ri, xfj, densty, deltaV
  Integer :: i, j, l, ind
  Integer, save :: naver=0

  naver = naver + 1

  ! Bin pair distances into histogram
  Do i = 1, natoms - 1
     Do j = i + 1, natoms
        rd(:) = R(:,i) - R(:,j)
        rr2 = dist2(rd)
        If (rr2 < rcut2) Then
           ind = Nint(Sqrt(rr2)/deltagr)
           if (ind >= 1 .and. ind <= nmaxgr) then
              histomix(ind,iatype(i),iatype(j)) = histomix(ind,iatype(i),iatype(j)) + 1
           endif
        Endif
     End Do
  End Do

  ! Normalize histograms to obtain g_ij(r)
  densty = natoms / v0
  Do i = 1, nmaxgr
     ri = i * deltagr
     If (ndim == 3) Then
        deltaV = 4.0d0*pi*((ri+deltagr/2.0d0)**3 - (ri-deltagr/2.0d0)**3) / 3.0d0
     Else
        deltaV = pi*((ri+deltagr/2.0d0)**2 - (ri-deltagr/2.0d0)**2)
     End If
     Do j = 1, nsp
        Do l = j, nsp
           xfj = Real(ntype(l), wp) / Real(natoms, wp)
           gmix(i,j,l) = (j/l+1) * histomix(i,j,l) / &
                (deltaV * ntype(j) * Real(naver, wp) * densty * xfj)
        End Do
     End Do
  End Do
End Subroutine structure
