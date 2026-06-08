Module readconf
    use set_precision
    use configuration, only : a, b, c, side, side2, v0, r_unit, &
        & q, r, natoms, ntype, nsp, ndim, atoms, iatype
    use rundata, only : iocfg, initcf
contains
  subroutine dlplmp_readconf
    implicit none
    Integer :: keytrj, imcon, iatm, i, j, dumm
    real(wp) :: dumx, dumy, dumz, rdum(3), xl, xh, yl, yh, zl, zh
    Character :: atms*8
    if ( initcf == "dlp") then
       Open (iocfg,file='data/CONFIG',status='old')
       Read(iocfg,'(1x)')
       Read(iocfg,*) keytrj, imcon
       Read(iocfg,*) a
       Read(iocfg,*) b
       Read(iocfg,*) c
       If (a(2)**2+a(3)**2+b(1)**2+b(3)**3+c(1)**2+c(2)**(2) > 1.d-7) Then
          Print *, " *** Error: Non orthorombic cells not implemented "
          Stop
       End If
       ntype(:) = 0
       Do i=1,natoms
          Read(iocfg,"(a8,i10)")atms, iatm
          Do j=1, nsp
             If (adjustl(trim((adjustl(atms)))) == trim((adjustl(atoms(j))))) Then
                iatype(iatm) = j
                ntype(j) = ntype(j)+1
             End If
          End Do
          If (keytrj.Eq.0) Then
             Read(iocfg,*) R(1:ndim,iatm)
          Else If (keytrj.Eq.1) Then
             Read(iocfg,*) R(1:ndim,iatm)
             Read(iocfg,*) dumx,dumy,dumz
          Else If (keytrj.Eq.2) Then
             Read(iocfg,*) R(1:ndim,iatm)
             Read(iocfg,*) dumx,dumy,dumz
             Read(iocfg,*) dumx,dumy,dumz
          Endif
       End Do
    else if (initcf == "lmp") then
       a(:)=0
       b(:)=0
       c(:)=0
       Open (iocfg,file='data/data.atoms',status='old')
       read (iocfg,'(1x,//////)')
       !
       !Only orthogonal cells
       !
       read (iocfg,*) xl, xh
       read (iocfg,*) yl, yh
       read (iocfg,*) zl, zh
       a(1) = xh-xl
       b(2) = yh-yl
       c(3) = zh-zl

       read(iocfg,'(1x,///////)')
       do i = 1, natoms
          read(iocfg,*) j, dumm, iatm,  q(j),  R(:,j) 
          iatype(j) = iatm
          ntype(iatm) = ntype(iatm)+1
       enddo
    endif
    close(iocfg)
  end subroutine dlplmp_readconf

end module readconf
