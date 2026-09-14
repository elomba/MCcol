!===============================================================================
! Program: gpMC (General Purpose Monte Carlo)
! Project: MCcol
!
! Authors:
!   Enrique Lomba (enrique.lomba@csic.es)
!   Eva G. Noya   (eva.noya@iqf.csic.es)
!   Instituto de Química Física Rocasolano / IQFR-CSIC, Madrid, Spain
!
! Description:
!   General Purpose atomistic Monte Carlo code for condensed-matter systems.
!   Simulates bulk multicomponent and ionic soft spherical systems with
!   arbitrary stoichiometry.
!
! Key Capabilities:
!   - Ensembles:
!       NVT : Canonical ensemble (fixed N, V, T)
!       NpT : Isobaric-Isothermal ensemble (fixed N, P, T) with isotropic
!             ('isotr') or anisotropic ('ortho') orthorhombic cell fluctuations.
!   - Interactions:
!       Morse Potential (keyp = 1)
!       Lennard-Jones 12-6 Potential (keyp = 2)
!       Potential Truncation & Shifting
!   - Long-Range Electrostatics:
!       Full Ewald summation for periodic boundaries:
!         * Real-space complementary error function screening (erfc(kappa*r)/r)
!         * Reciprocal space Fourier sum over half-space wavevectors (symmetry optimized)
!         * Electrostatic self-energy correction
!   - Algorithmic Speedups:
!       * 3D Link-Cell domain decomposition for O(N) neighbor searches
!       * Paul Breeuwsma smooth cubic spline potential interpolation
!       * Incremental structure factor update for single-particle reciprocal moves O(K)
!   - Interoperability:
!       * Input: DL_POLY 2 CONFIG format (initcf="dlp"), LAMMPS data.atoms format (initcf="lmp")
!       * Output: LAMMPS custom trajectory dump (gpMC.lammpstrj, last.lammpstrj),
!                 DL_POLY CONFIG.last, binary checkpoint restart.dmp
!       * Analysis: Instantaneous thermo (thermoins.dat), block averages (thermoaver.dat),
!                   multicomponent partial radial distribution functions g_ij(r) (gmix.dat).
!===============================================================================
Program gpMC
    use set_precision

    ! System potential parameters and electrostatics control
    Use potential, Only : keyp, kint, elect

    ! System particle count and topology
    Use configuration, Only : natoms

    ! Run control data, step limits, ensemble flags, and I/O units
    Use rundata, Only : kT, restart, nequil, nstep, nb, ensemble, npgr,&
         & s_cput, ntraj, istep, istep_ini, iotrj, ilong 

    ! Accumulated and instantaneous thermodynamic properties
    Use properties

    ! 3D Link-Cell configuration flag
    Use linkcell, Only : use_cell

    ! Explicit interface blocks
    Use interfaces, Only : move_natoms, fpot_elecLJ, fpot_elecMorse, fpot_Morse, fpot_LJ

    ! Initialization procedures
    Use Init, Only : Init_conf, Init_pot, Init_rundata, Init_interp

    ! Thermodynamic averaging routines
    Use Thermo, Only : Averages

    ! Output and logging routines
    Use Output, Only : Printout, init_printout, initout, run_info,&
        & printgr, end_printout, print_ener
    Use WriteCfg, only : dump_trj

    ! Energy calculation subroutines (direct sum and link-cell)
    Use Energy, Only : Energ, Energ_cell, Eshort_r

    ! Link-cell spatial grid builder
    Use Cells, Only : build_cells, init_cell

    ! Isobaric-isothermal volume change driver
    Use VolumeChange, Only : move_volume

    ! High-resolution timing utility
    Use Util, Only : cputime

    Implicit None
    Integer :: j, ntest
    Real (dkind) :: esrold, fourold

    ! Intercept POSIX signals (SIGTERM, SIGINT) for orderly shutdown & checkpointing
    call catch()
    ! Get initial CPU time
    s_cput = cputime()
    !
    ! Initialize particle configuration
    !
    Call Init_conf
    !
    ! Read in parameters for run
    !
    Call Init_rundata
    !
    ! Program information
    !
    call initout
    if (.not. restart) then
        !
        ! Initialize potential parameters
        !
        Call Init_pot
        !
        ! Initialize interpolation tables for pair potentials
        !
        if (elect) then
           if (kint == 1) then
              Call Init_interp(fpot_elecMorse)
           else
              Call Init_interp(fpot_elecLJ)
           endif
        else
           if (kint == 1) then
              Call Init_interp(fpot_Morse)
           else
              Call Init_interp(fpot_LJ)
           endif
        endif
    endif
    !
    ! Initialize and build link cells (if possible), controlled by use_cell)
    !
    if (use_cell) Call Init_cell
    !
    ! Print run specific info
    !
    call run_info
    if (use_cell) then
       call build_cells
       if (elect) then
           if (kint == 1) then
              Call Energ_cell(fpot_elecMorse)
              Call Eshort_r(fpot_Morse,Evdw)
           else
              Call Energ_cell(fpot_elecLJ)
              Call Eshort_r(fpot_LJ,Evdw)
           endif
        else
           if (kint == 1) then
              Call Energ_cell(fpot_Morse)
           else
              Call Energ_cell(fpot_LJ)
           endif
        endif
    else
       if (elect) then
           if (kint == 1) then
              Call Energ(fpot_elecMorse)
           else
              Call Energ(fpot_elecLJ)
           endif
        else
           if (kint == 1) then
              Call Energ(fpot_Morse)
           else
              Call Energ(fpot_LJ)
           endif
        endif
     endif
    call print_ener
    !
    ! Initialize printout
    !
    call init_printout
    !
    ! Store initial short range energy and Fourier component values
    !
    Esrold = E_sr
    Fourold = E_Fourier
    !
    ! Run nstep configuration generations (natoms*nstep atom
    ! displacements at present )
    !
    Do istep=1+istep_ini, nstep+istep_ini
       !   if (use_cell) call build_cells
       Call move_natoms(natoms)

        !
        ! Insert here particle insertions/deletions, volume changes, etc ..
        !
        if(Ensemble == 'npt') Call move_volume
        !
        ! Perform averages when equilibration has been reached.
        !
        If (istep >= nequil) Then
            If (Mod(istep-istep_ini,npgr).Eq.0) Then
                Call printgr
            Endif
            !
            ! Dump trajectory file if needed
            if (ntraj .ne. 0) then
                if (mod(istep-istep_ini,ntraj) .Eq. 0) then
                    call dump_trj(istep,iotrj)
                endif
            endif
        End If
        If (Mod(istep-istep_ini,nb).Eq.0) Then
           if (kint==1) then
              call Eshort_r(fpot_Morse,Evdw)
           else
              call Eshort_r(fpot_LJ,Evdw)
           endif
            If (istep > nequil) Then
                Call Averages
                Call structure
            Else
                !
                ! Print out instantaneous values.
               !
                Call Printout(.false.)
               !
            End If
        End If
     End Do
    !
    ! Calculate energy from last configuration (consistency check).
    !
    if (use_cell) then
       if (elect) then
           if (kint == 1) then
              Call Energ_cell(fpot_elecMorse)
           else
              Call Energ_cell(fpot_elecLJ)
           endif
        else
           if (kint == 1) then
              Call Energ_cell(fpot_Morse)
           else
              Call Energ_cell(fpot_LJ)
           endif
        endif
    else
       if (elect) then
           if (kint == 1) then
              Call Energ(fpot_elecMorse)
           else
              Call Energ(fpot_elecLJ)
           endif
        else
           if (kint == 1) then
              Call Energ(fpot_Morse)
           else
              Call Energ(fpot_LJ)
           endif
        endif
    Endif
    call print_ener
    Call end_printout
    ! Dump restart file
    call cierra(1)
End Program gpMC






