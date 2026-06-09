module WriteCfg
  !
  ! Dump final configuration and trajectory files in different formats
  !
  use configuration, only : a, b, c, side, r, natoms, ntype,&
       &atoms, iatype, nsp, r_unit, ndim
  use rundata, only: iocfg, iotrj
contains
  Subroutine writecfg_dlp
    ! DLPOLY format
    Implicit None
    Integer :: keytrj=0, imcon=1, iatm, i, j
    Open (iocfg,file='CONFIG.last')
    Write(iocfg,'(1x)')
    Write(iocfg,*) keytrj, imcon
    Write(iocfg,*) a
    Write(iocfg,*) b
    Write(iocfg,*) c
    do i=1, natoms
       write(iocfg,"(a8,i10)")atoms(iatype(i)),i
       write(iocfg,"(3g20.10)")R(1:ndim,i)*r_unit(1:ndim)
    enddo
    close(iocfg)
  End Subroutine writecfg_dlp

  Subroutine writecfg_lmp(istep)
    ! LAMMPS format dump custom id type mol x y z 
    Implicit None
    Integer :: keytrj=0, imcon=1, iatm, i, j, istep
    Open (iocfg,file='last.lammpstrj')
    call dump_trj(istep,iocfg)
    close(iocfg)
  End Subroutine writecfg_lmp


  Subroutine dump_trj(istep,iocfg)
    implicit none
    integer,intent(IN) :: istep, iocfg
    integer :: i
    Write(iocfg,'("ITEM: TIMESTEP")')
    Write(iocfg,*) istep
    Write(iocfg,'("ITEM: NUMBER OF ATOMS")')
    Write(iocfg,*) natoms
    Write(iocfg,'("ITEM: BOX BOUNDS pp pp pp")')
    Write(iocfg,*) -a(1)/2,a(1)/2
    Write(iocfg,*) -b(2)/2,b(2)/2
    Write(iocfg,*) -c(3)/2,c(3)/2
    Write(iocfg,'("ITEM: ATOMS id type mol x y z")')
    do i=1, natoms
       write(iocfg,"(i10,i2,i2,3g20.10)")i,iatype(i),iatype(i),R(1:ndim,i)*r_unit(1:ndim)
    enddo
  end subroutine dump_trj

end module WriteCfg
