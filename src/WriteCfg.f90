!===============================================================================
! Module: WriteCfg
!
! Purpose:
!   Writes final configuration snapshots and continuous trajectory frames in
!   DL_POLY 2 or LAMMPS custom formats.
!
! Supported formats:
!   - writecfg_dlp : Writes CONFIG.last in standard DL_POLY 2 format.
!   - writecfg_lmp : Writes last.lammpstrj in LAMMPS custom dump format.
!   - dump_trj     : Writes an individual trajectory frame to an open unit.
!===============================================================================
module WriteCfg
  use configuration, only : a, b, c, side, r, natoms, ntype,&
       &atoms, iatype, nsp, r_unit, ndim
  use rundata, only: iocfg, iotrj
contains

  !-----------------------------------------------------------------------------
  ! Subroutine: writecfg_dlp
  !
  ! Purpose:
  !   Writes the final configuration to 'CONFIG.last' in standard DL_POLY 2 format:
  !     Line 1: Header / Title record
  !     Line 2: keytrj (0), imcon (1)
  !     Lines 3-5: Box lattice vectors a, b, c in Angstroms
  !     Per-atom records: Species label and atom index, followed by unscaled (x, y, z).
  !-----------------------------------------------------------------------------
  Subroutine writecfg_dlp
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

  !-----------------------------------------------------------------------------
  ! Subroutine: writecfg_lmp
  !
  ! Purpose:
  !   Writes final configuration to 'last.lammpstrj' in LAMMPS custom dump format.
  !-----------------------------------------------------------------------------
  Subroutine writecfg_lmp(istep)
    Implicit None
    Integer :: keytrj=0, imcon=1, iatm, i, j, istep
    Open (iocfg,file='last.lammpstrj')
    call dump_trj(istep,iocfg)
    close(iocfg)
  End Subroutine writecfg_lmp

  !-----------------------------------------------------------------------------
  ! Subroutine: dump_trj
  !
  ! Purpose:
  !   Writes a trajectory snapshot to the specified logical unit (iocfg) in
  !   standard LAMMPS dump format ("ITEM: TIMESTEP", "ITEM: NUMBER OF ATOMS",
  !   "ITEM: BOX BOUNDS", "ITEM: ATOMS id type mol x y z").
  !
  ! Arguments:
  !   istep (in) : Current simulation timestep / sweep counter.
  !   iocfg (in) : Logical file unit number to write to.
  !-----------------------------------------------------------------------------
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
