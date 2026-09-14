!===============================================================================
! Module: readconf
!
! Purpose:
!   Provides configuration reading routines for initial particle coordinates,
!   box geometries, and atom type assignments.
!
! Supported Formats:
!   1. DL_POLY 2 format (initcf = "dlp", file 'CONFIG'):
!      - Reads record title, keytrj, imcon, lattice vectors a, b, c.
!      - Checks orthorhombic condition (diagonal lattice matrix).
!      - Reads atom labels and coordinates.
!   2. LAMMPS data format (initcf = "lmp", file 'data.atoms'):
!      - Parses header lines for 'atoms', 'atom types', and box bounds:
!        [xlo, xhi], [ylo, yhi], [zlo, zhi].
!      - Verifies that atom types and count match system.dat definitions.
!      - Supports LAMMPS full (id mol type q x y z), molecular (id mol type x y z),
!        and atomic (id type x y z) styles.
!      - Re-centers coordinates relative to box midpoint.
!===============================================================================
Module readconf
    use set_precision
    use configuration, only : a, b, c, side, side2, v0, r_unit, &
        & q, r, natoms, ntype, nsp, ndim, atoms, iatype
    use rundata, only : iocfg, initcf, c_red, c_reset 
contains

  !-----------------------------------------------------------------------------
  ! Subroutine: dlplmp_readconf
  !
  ! Purpose:
  !   Dispatcher subroutine that reads initial configuration files based on the
  !   `initcf` variable ("dlp" for DL_POLY CONFIG or "lmp" for LAMMPS data.atoms).
  !-----------------------------------------------------------------------------
  subroutine dlplmp_readconf
    implicit none
    Integer :: keytrj, imcon, iatm, i, j, dumm
    real(wp) :: dumx, dumy, dumz, rdum(3), xl, xh, yl, yh, zl, zh, rcf(3)
    Character :: atms*8
    Character(len=256) :: line
    Integer :: ios, n_types_in_file, natms, imol
    Real(wp) :: q_atom
    if ( initcf == "dlp") then
       Open (iocfg,file='CONFIG',status='old')
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
      Open (iocfg,file='data.atoms',status='old')
      !
      Do
         Read(iocfg,'(A)', iostat=ios) line
         If (ios /= 0) Exit 
         If (len_trim(line) == 0) Cycle
         If (index(line, 'atom types') > 0) Then
            ! Read the integer at the beginning of the line
            Read(line, *) n_types_in_file
      
            ! Compare it with the namelist variable
            If (n_types_in_file /= nsp) Then
               Print *, c_red
               Print *, "=========================================================================="
               Print *, " *** FATAL SECURITY ERROR: Mismatch in Particle Types! ***"
               Print *, "=========================================================================="
               Print *, " -> nsp defined in input system.dat : ", nsp
               Print *, " -> 'atom types' defined in data.atoms    : ", n_types_in_file
               Print *, ""
               Print *, " Please ensure your input namelist matches your configuration file."
               Print *, " Simulation aborted to prevent reckless behavior."
               Print *, "=========================================================================="
               Print *, c_reset
               Stop
            End If
         End If 

         If (index(line, 'atoms') > 0 .and. index(line, 'atom types') == 0) Then
            read(line, *, iostat=ios) natms
            if (ios == 0) then
               if (natms .ne. natoms) then
                  write(*,"(A,' ** Error natoms in data.atoms .ne. to system.dat !!',A)")c_red, c_reset
                  stop
               endif
            endif
         End If

         If (index(line, 'xlo xhi') > 0) Then; read(line, *) xl, xh; End If
         If (index(line, 'ylo yhi') > 0) Then; read(line, *) yl, yh; End If
         If (index(line, 'zlo zhi') > 0) Then; read(line, *) zl, zh;  End If
         ntype(:) = 0
         If (index(adjustl(line), 'Atoms') == 1) Then
            !Only orthogonal cells
            !
            a(1) = xh-xl
            b(2) = yh-yl
            c(3) = zh-zl
            Do i = 1, natoms
               Read(iocfg, '(A)') line
               Do while (len_trim(line) == 0)
                  Read(iocfg, '(A)') line
               End do
               ! Try reading as full style: id mol type q x y z (7 values)
               read(line, *, iostat=ios) j, imol, iatm, q_atom, Rcf(1), Rcf(2), Rcf(3)
               if (ios /= 0) then
                  ! Try reading as molecular style: id mol type x y z (6 values)
                  read(line, *, iostat=ios) j, imol, iatm, Rcf(1), Rcf(2), Rcf(3)
                  if (ios /= 0) then
                     ! Fallback to atomic style: id type x y z (5 values)
                     read(line, *) j, iatm, Rcf(1), Rcf(2), Rcf(3)
                  endif
               else
                  q(j) = q_atom
               endif
               iatype(j) = iatm 
               ntype(iatm) = ntype(iatm)+1
               ! Center coordinates
               R(1,j) = Rcf(1)-(xh+xl)/2
               R(2,j) = Rcf(2)-(yh+yl)/2
               R(3,j) = Rcf(3)-(zh+zl)/2
            End Do
 
         End If
 
      end do
   endif
   close(iocfg)
  end subroutine dlplmp_readconf

end module readconf
