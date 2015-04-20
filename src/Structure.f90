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
  naver = naver+1
  Do i=1, natoms-1
     Do j=i+1, natoms
        rd(:) = r(i,:)-r(j,:)
        rr2 = dist2(rd)
        If (rr2 < rcut2) Then
           ind = Nint(Sqrt(rr2)/deltagr)
           histomix(ind,iatype(i),iatype(j)) =  histomix(ind,iatype(i),iatype(j)) + 1
        Endif
     End Do
  End Do
  densty = natoms/v0
  Do i=1, nmaxgr
     ri = i*deltagr
     If (ndim == 3) Then
        deltaV = 4*pi*((ri+deltagr/2)**3-(ri-deltagr/2)**3)/3.0
     Else
        deltaV = pi*((ri+deltagr/2)**2-(ri-deltagr/2)**2)
     End If
     Do j=1,nsp
        Do l=j,nsp
           xfj = Real(ntype(l))/Real(natoms)
           gmix(i,j,l) = (j/l+1)*histomix(i,j,l)/(deltaV*ntype(j)*Naver&
                &*densty*xfj)
        End Do
     End Do
  End Do
End Subroutine structure
