!
! General purpose Monte Carlo code for atomistic simulations. At
! present simulates bulk systems composed of soft spherical particles
! (single or multicomponent) with or without charges. Charges are
! treated using Ewald sums. The program implements a link cell
! algorithm if the sample size allows it.
! Implemented ensembles: NVT
!
! Input data files (see the files for parameter specifications) :
!           system.dat : contains description of the system to be simulated
!           runMC.dat :  contains specific parameters that control
!                        the run
!           CONFIG:      initial configuration in DLPOLY 2 format (to
!                         be generalized)
! Output files:
!           thermoaver.dat : thermodynamic averages
!           thermoins.dat  : instantaneous thermodynamic quantities
!           gmix.dat       : pair correlation functions
!
! Program units:
!         Energy: "eV" electronVolts
!                 "K" Lennard-Jones (Energy in Kelvin)
!         Distance: Angstrom (internal units in box length)
!
!  General modules are contained in Definitions.f90
!  Internal control variables:
!           use_cell (logical): if .true. link cells are used
!
!           elect    (logical): if .true. Ewald electrostatics is computed
!
!   An interpolation algorithm is used to evaluate interactions.
!
! E. Lomba, April, 2015
!
!
Program gpMC
    !
    ! Module to define generic precision
    !
    use set_precision
    !--------- general Definitions --------------------------
    !
    ! Interaction potential parameters.
    !       keyp(nit) : integer array                      = 1 Morse
    !                   (one value for every interaction)  = 2 Lennard-Jones
    !
    Use potential, Only : keyp
    ! Variables defining the system configuration
    Use configuration, Only : natoms
    ! Control parameters for the run
    Use rundata, Only : nequil, nstep, nb, ensemble, npgr, s_cput, ntraj
    ! Initialization routies (including input of data)
    ! System properties
    Use properties
    ! Parameters that define link cells
    Use linkcell, Only : use_cell
    ! Interfaces to routines
    Use interfaces, Only : move_natoms, histograms
    !
    !----------- Specific routines --------------------------
    !
    Use Init, Only : Init_conf, Init_pot, Init_rundata, Init_interp
    ! Thermodynamics averages
    Use Thermo, Only : Averages
    ! Output routines
    Use Output, Only : Printout, init_printout, initout, run_info,&
        & printgr, end_printout
    Use WriteCfg, only : dump_trj
    ! Routines to calculate the energy (with and without link cell)
    Use Energy, Only : Energ, Energ_cell
    ! Link cell routines
    Use Cells, Only : build_cells, init_cell
    ! Utility routines
    Use Util, Only : cputime
    Implicit None
    Integer :: istep, j, ntest
    Real (wp) :: esrold, fourold
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
    !
    ! Initialize potential parameters
    !
    Call Init_pot
    !
    ! Initialize interpolation tables for pair potentials
    !
    Call Init_interp
    !
    ! Initialize and build link cells (if possible), controlled by use_cell)
    !
    Call Init_cell
    !
    ! Print run specific info
    !
    call run_info
    if (use_cell) then
        call build_cells
        Call Energ_cell
    else
        Call Energ
    endif
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
    Do istep=1, nstep
        !   if (use_cell) call build_cells
        Call move_natoms(natoms)

        !
        ! Insert here particle insertions/deletions, volume changes, etc ..
        !

        !
        ! Perform averages when equilibration has been reached.
        !
        If (istep >= nequil) Then
            Call histograms(ensemble)
            If (Mod(istep,npgr).Eq.0) Then
                Call printgr
            Endif
            !
            ! Dump trajectory file if needed
            if (ntraj .ne. 0) then
                if (mod(istep,ntraj)) then
                    call dump_trj(istep)
                endif
            endif
        End If
        If (Mod(istep,nb).Eq.0) Then
            If (istep >= nequil) Then
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
        Call energ_cell
    else
        Call energ
    Endif
    Call end_printout

End Program gpMC






