Module readconf
    use set_precision
    use configuration, only : a, b, c, side, side2, v0, r_unit, &
        & r, natoms, ntype, nsp, ndim, atoms, iatype
    use rundata, only : iocfg
contains
    subroutine dlp_readconf
        implicit none
        Integer :: keytrj, imcon, iatm, i, j
        real(wp) :: dumx, dumy, dumz
        Character :: atms*8
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
                Read(iocfg,*) R(iatm,1:ndim)
            Else If (keytrj.Eq.1) Then
                Read(iocfg,*) R(iatm,1:ndim)
                Read(iocfg,*) dumx,dumy,dumz
            Else If (keytrj.Eq.2) Then
                Read(iocfg,*) R(iatm,1:3)
                Read(iocfg,*) dumx,dumy,dumz
                Read(iocfg,*) dumx,dumy,dumz
            Endif
        End Do

        close(iocfg)
    end subroutine dlp_readconf

end module readconf
