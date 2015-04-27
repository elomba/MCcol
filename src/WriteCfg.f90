module WriteCfg
    !
    ! Dump final configuration and trajectory files in different formats
    !
    use configuration, only : a, b, c, side, r, natoms, ntype,&
        &atoms, iatype, nsp, r_unit, ndim
    use rundata, only: iocfg, iotrj
contains
    Subroutine writecfg_dlp
    !
    ! Dump last configuration in DLPOLY format (usable as starting configuration)
    !
        Implicit None
        Integer :: keytrj=0, imcon=1, iatm, i, j
        Open (iocfg,file='results/CONFIG.last')
        Write(iocfg,'(1x)')
        Write(iocfg,*) keytrj, imcon
        Write(iocfg,*) a
        Write(iocfg,*) b
        Write(iocfg,*) c
        do i=1, natoms
            write(iocfg,"(a8,i10)")atoms(iatype(i)),i
            write(iocfg,"(3g20.10)")r(i,1:ndim)*r_unit(:)
        enddo
        close(iocfg)
    End Subroutine writecfg_dlp

    Subroutine writecfg_xyz(istep)
    !
    ! Dump last conf in .xyz format
    !
        Implicit None
        integer,intent(IN) :: istep
        Integer :: i
        Open (iocfg,file='results/last.xyz')
        Write(iocfg,*) natoms
        Write(iocfg,'("# step=",i10," cell dims",3f15.7,&
           &" tstep ",f10.7,3(" no. part. tipo ",i1,i6:))'        )istep,&
            & side(1:ndim),0.0,(i,ntype(i),i=1,nsp)
        Do i=1, natoms
            Write(iocfg,'(i4,3f15.7)')iatype(i),R(i,:)*r_unit(:)
        enddo
        close(iocfg)
    End Subroutine writecfg_xyz

    subroutine dump_trj(istep)
    !
    ! dump trajectory in .xyz format (can be read by VMD)
    !
        implicit none
        integer,intent(IN) :: istep
        integer :: i
        Write(iotrj,*) natoms
        Write(iotrj,'("# step=",i10," cell dims",3f15.7,&
           &" tstep ",f10.7,3(" no. part. tipo ",i1,i6:))'        )istep,&
            & side(1:ndim),0.0,(i,ntype(i),i=1,nsp)
        Do i=1, natoms
            Write(iotrj,'(i4,3f15.7)')iatype(i),R(i,:)*r_unit(:)
        enddo
    end subroutine dump_trj

end module WriteCfg
