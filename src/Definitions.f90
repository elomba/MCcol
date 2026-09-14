!===============================================================================
! File: Definitions.f90
!
! Purpose:
!   Contains primary global data definitions, simulation parameters, and
!   explicit interface blocks for the gpMC Monte Carlo engine.
!
! Modules defined:
!   - configuration : Particle coordinates, box geometry, atom types, species.
!   - potential     : Potential parameters (LJ / Morse), cutoffs, Ewald terms,
!                     physical constants, and unit conversion factors.
!   - rundata       : Run control parameters, I/O unit numbers, acceptance counts,
!                     ensemble flags, and ANSI terminal color formatting.
!   - properties    : Instantaneous and accumulated thermodynamic properties
!                     (energies, volume, box sides, radial distribution data).
!   - linkcell      : Data structures for the 3D link-cell / cell-list domain
!                     decomposition method (O(N) neighbor searching).
!   - interp        : Look-up tables and Paul Breeuwsma cubic spline matrix for
!                     rapid pair potential interpolation.
!   - interfaces    : Explicit Fortran 90 interfaces for external procedures.
!===============================================================================

!-------------------------------------------------------------------------------
! Module: configuration
! Purpose: Stores structural configuration, box dimensions, and atom attributes.
!-------------------------------------------------------------------------------
module configuration
    use set_precision
    implicit none

    integer, parameter :: ndim = 3         ! Spatial dimensions (always 3 for bulk)
    integer            :: length           ! Number of integer words required for RNG seed
    integer            :: nsp              ! Number of distinct chemical species
    integer            :: natoms           ! Total number of atoms in the system

    integer, allocatable :: ntype(:)       ! Number of atoms per species [1:nsp]
    integer, allocatable :: iatype(:)      ! Species ID for each atom [1:natoms]
    integer, allocatable :: seed(:)        ! Random number generator seed array [1:length]

    real(wp), dimension(:), allocatable :: q   ! Partial charge per atom [1:natoms] (electrons)
    real(wp), dimension(:), allocatable :: qsp ! Reference charge per species [1:nsp] (electrons)
    real(wp), allocatable :: r(:,:)        ! Scaled atom coordinates [1:ndim, 1:natoms] in [-0.5, 0.5)

    ! Simulation box vectors and dimensions:
    real(wp), dimension(ndim) :: a         ! First box vector (a_x, 0, 0) in Angstroms
    real(wp), dimension(ndim) :: b         ! Second box vector (0, b_y, 0) in Angstroms
    real(wp), dimension(ndim) :: c         ! Third box vector (0, 0, c_z) in Angstroms
    real(wp), dimension(ndim) :: side      ! Box side lengths (Lx, Ly, Lz) in Angstroms
    real(wp), dimension(ndim) :: side2     ! Half-box side lengths (Lx/2, Ly/2, Lz/2)
    real(wp), dimension(ndim) :: r_unit    ! Box edge scaling factors for unscaling coordinates
    real(wp)                  :: v0        ! Total box volume V = Lx * Ly * Lz in Angstroms^3

    character, dimension(:), allocatable :: atoms*8  ! Alphanumeric labels for species [1:nsp]
end module configuration


!-------------------------------------------------------------------------------
! Module: potential
! Purpose: Force field parameters, electrostatics constants, and Ewald sums.
!-------------------------------------------------------------------------------
module potential
    use set_precision
    implicit none

    character :: units*8                   ! Energy units label: 'eV', 'K', or 'kcal/mol'
    character, dimension(:), allocatable :: pot*4  ! Potential model label per interaction ('lj  ', 'mors')
    integer, dimension(:), allocatable   :: keyp   ! Numeric key: 1 = Morse, 2 = Lennard-Jones
    integer :: kint                        ! Primary interaction type key (1 = Morse, 2 = LJ)
    integer :: nitmax                      ! Total number of unique pairs = nsp*(nsp+1)/2
    integer, dimension(:,:), allocatable :: itp    ! Symmetric pair interaction mapping matrix [nsp, nsp] -> nit

    ! Ewald summation parameters and cutoffs:
    integer :: kmx, kmy, kmz               ! Maximum k-vector index along x, y, z
    integer :: kmt                         ! Total number of reciprocal space k-vectors
    real(wp) :: dospix, dospiy, dospiz     ! 2*pi / L_x, 2*pi / L_y, 2*pi / L_z
    real(wp) :: rcpcut, rcpcut2            ! Reciprocal-space spherical cutoff and squared cutoff (1/Angstrom)
    real(wp) :: rcut, rcut2                ! Real-space spherical potential cutoff and squared cutoff (Angstrom)
    real(wp) :: kappa                      ! Ewald Gaussian screening parameter (1/Angstrom)
    real(wp) :: qtotal                     ! Total net charge in system (should be ~0 for neutrality)
    real(dkind) :: selfe = 0.0d0           ! Ewald electrostatic self-energy term
    real(dkind) :: pi2                     ! 2 * pi constant

    ! Pairwise force field parameters (matrices and condensed vectors):
    real(wp), dimension(:,:), allocatable :: aa  ! Well depth: epsilon (LJ) or D_e (Morse)
    real(wp), dimension(:,:), allocatable :: bb  ! Sigma (LJ) or stiffness alpha (Morse)
    real(wp), dimension(:,:), allocatable :: cc  ! Equilibrium distance r_0 (Morse)
    real(wp), dimension(:,:), allocatable :: rc  ! Cutoff radius r_c per pair in Angstroms
    real(wp), dimension(:),   allocatable :: rc2 ! Squared cutoff r_c^2 per pair index [1:nitmax]
    real(wp), dimension(:),   allocatable :: al  ! Scaled well depth divided by kT
    real(wp), dimension(:),   allocatable :: bl  ! Interaction parameter (stiffness / sigma)
    real(wp), dimension(:),   allocatable :: cl  ! Morse r_0 parameter
    real(wp), dimension(:),   allocatable :: bl2 ! Squared sigma (bl^2) for LJ
    real(wp), dimension(:),   allocatable :: qprod ! Charge product q_i * q_j per pair
    real(wp), dimension(:),   allocatable :: ucut  ! Value of potential at cutoff for shifting

    ! Reciprocal space wave-vectors and phase factors:
    real(wp), dimension(:),   allocatable :: km2   ! Squared magnitude |k|^2 for reciprocal vectors
    real(wp), dimension(:),   allocatable :: ekm2  ! Ewald reciprocal weight factor: exp(-k^2/(4*kappa^2)) / k^2
    real(wp), dimension(:),   allocatable :: kr    ! Current k-vector components
    complex(wp), dimension(:), allocatable :: rhokk   ! Total structure factor rho(k)
    complex(wp), dimension(:), allocatable :: deltann ! Structure factor change delta_rho(k)
    complex(wp), dimension(:), allocatable :: einx, einy, einz ! Single-particle trial phase factors
    complex(wp), dimension(:,:), allocatable :: eix, eiy, eiz  ! Multi-particle phase factors [natoms, -k:k]

    ! Physical Constants:
    ! ctreV   : e^2 / (4 * pi * epsilon_0) in eV * Angstrom = 14.39964361
    ! pi      : Circle constant pi = 3.141592653589793
    ! kbeV    : Boltzmann constant in eV/K = 8.6173324e-5
    ! ev2K    : Conversion factor from eV to Kelvin = 11604.58857702
    ! ev2Kcal : Conversion factor from eV to kcal/mol = 23.06054195
    ! kbKcal  : Boltzmann constant in kcal/(mol*K) = 0.001987191686
    real(wp) :: ctr                        ! Coulomb prefactor scaled by 1/kT
    real(dkind), parameter :: ctreV   = 14.39964361d0
    real(dkind), parameter :: pi      = 3.141592653589793d0
    real(dkind), parameter :: kbeV    = 8.6173324d-5
    real(dkind), parameter :: ev2K    = 11604.58857702d0
    real(dkind), parameter :: ev2Kcal = 23.06054195d0
    real(dkind), parameter :: kbKcal  = 0.001987191686d0

    ! Pressure conversion factors from bar to internal units:
    ! bar2eV   : bar -> eV / Angstrom^3
    ! bar2K    : bar -> Kelvin / Angstrom^3
    ! bar2Kcal : bar -> kcal / (mol * Angstrom^3)
    real(wp), parameter :: bar2eV   = 0.624150932d-6
    real(wp), parameter :: bar2K    = 0.72429715652d-2
    real(wp), parameter :: bar2Kcal = 1.4393258781725146d-5

    complex(wp), parameter :: ii = (0.0d0, 1.0d0) ! Imaginary unit sqrt(-1)
    logical :: elect   ! True if Ewald electrostatics is activated
    logical :: pshift  ! True if short-range potential is shifted to zero at r_cut
    integer :: fou_type ! Reciprocal method: 1 = Standard Ewald sum, 2 = PME
end module potential


!-------------------------------------------------------------------------------
! Module: rundata
! Purpose: Controls run parameters, MC step counters, and file I/O unit numbers.
!-------------------------------------------------------------------------------
module rundata
    use configuration
    implicit none

    ! Logical I/O unit numbers
    integer, parameter :: iosys = 7  ! Input system specification (system.dat)
    integer, parameter :: iocfg = 8  ! Configuration I/O (CONFIG, data.atoms)
    integer, parameter :: iotrj = 9  ! Trajectory output (gpMC.lammpstrj)
    integer, parameter :: iorun = 10 ! Run parameters input (runMC.dat)
    integer, parameter :: ioth  = 11 ! Thermodynamic block averages (thermoaver.dat)
    integer, parameter :: iothi = 12 ! Instantaneous thermodynamics (thermoins.dat)
    integer, parameter :: igr   = 13 ! Radial distribution function output (gmix.dat)
    integer, parameter :: ilong = 20 ! Max length of file path strings

    logical :: restart               ! True if resuming simulation from restart.dmp
    integer :: nstep                 ! Total number of Monte Carlo production sweeps
    integer :: nequil                ! Number of equilibration sweeps before averaging
    integer :: nb                    ! Sampling / block average interval (sweeps)
    integer :: npgr                  ! RDF g(r) accumulation/output interval (sweeps)
    integer :: ntraj                 ! Trajectory dump interval (sweeps, 0 = disable)
    integer :: istep                 ! Current MC sweep counter
    integer :: istep_ini = 0         ! Initial sweep index (offset if restarting)

    ! Acceptance and trial move statistics:
    integer :: ntrial = 0            ! Total particle displacement trial moves attempted
    integer :: naccept = 0           ! Displacements accepted via Metropolis criterion
    integer :: nvaccept = 0          ! Volume trial moves accepted (NpT ensemble)
    integer :: naver = 0             ! Number of thermodynamic average samples collected

    real(wp) :: temp                 ! Simulation temperature in Kelvin
    real(wp) :: kT                   ! Thermal energy k_B * T in current energy units
    real(wp) :: pres                 ! Imposed pressure in bar (NpT ensemble)
    real(wp) :: rdmax(1:ndim)        ! Maximum trial displacement along x, y, z (Angstrom)
    real(wp) :: vdmax                ! Maximum trial volume change (Angstrom or volume)
    real(wp) :: s_cput, e_cput       ! Simulation start and end CPU times (seconds)

    character :: ensemble*3          ! 'nvt' = Canonical, 'npt' = Isobaric-Isothermal
    character :: initcf*3            ! Initial config format: 'dlp' (DL_POLY) or 'lmp' (LAMMPS)
    character :: scaling*5           ! Volume scaling mode: 'isotr' (isotropic) or 'ortho' (anisotropic)
    character :: stat*10 = "sequential" ! File access mode ("sequential" or "append")

    ! ANSI Terminal escape codes for colorized logging:
    character(len=8), parameter :: c_blue   = char(27)//'[1;34m'
    character(len=8), parameter :: c_cyan   = char(27)//'[1;36m'
    character(len=8), parameter :: c_green  = char(27)//'[1;32m'
    character(len=8), parameter :: c_yellow = char(27)//'[1;33m'
    character(len=8), parameter :: c_red    = char(27)//'[1;31m'
    character(len=8), parameter :: c_reset  = char(27)//'[0m'
end module rundata


!-------------------------------------------------------------------------------
! Module: properties
! Purpose: Accumulates instantaneous and running average thermodynamic properties.
!-------------------------------------------------------------------------------
module properties
    use set_precision
    use configuration, only : ndim
    implicit none

    ! Running thermodynamic averages (divided by naver at output):
    real(dkind) :: Etav    ! Average total energy
    real(dkind) :: E_sav   ! Average short-range energy
    real(dkind) :: E_lrav  ! Average long-range energy
    real(dkind) :: E_vdwav ! Average van der Waals / dispersion energy

    ! Instantaneous energy components:
    real(dkind) :: Etotal     ! Total system potential energy = E_sr + E_Fourier + selfe
    real(dkind) :: E_sr       ! Short-range pairwise energy (real space Coul + vdW)
    real(dkind) :: E_coulomb  ! Total electrostatic energy = E_Fourier + selfe + E_sr(Coul)
    real(dkind) :: E_fourier  ! Reciprocal space (Fourier) Ewald electrostatic energy
    real(dkind) :: Evdw       ! Pure van der Waals / dispersion energy
    real(dkind) :: virial     ! Internal virial accumulator

    ! Geometry averages for NpT simulations:
    real(dkind) :: Vol_av              ! Accumulated box volume average
    real(dkind) :: side_av(1:ndim)     ! Accumulated box dimensions average (Lx, Ly, Lz)

    ! Pair correlation / Radial Distribution Function g_ij(r):
    real(wp) :: deltagr                ! RDF bin width (Angstrom)
    integer  :: nmaxgr                 ! Total number of distance bins = rcut / deltagr
    integer,  dimension(:,:,:), allocatable :: histomix ! Distance pair histogram [bins, nsp, nsp]
    real(wp), dimension(:,:,:), allocatable :: gmix     ! Normalized RDF g_ij(r) [bins, nsp, nsp]
end module properties


!-------------------------------------------------------------------------------
! Module: linkcell
! Purpose: Manages 3D cell lists for accelerated O(N) neighbor searching.
!-------------------------------------------------------------------------------
module linkcell
    use set_precision
    implicit none

    real(wp) :: cellx, celly, cellz    ! Normalized cell dimensions along x, y, z
    integer, dimension(:,:), allocatable :: neigh ! 27 neighbor cell IDs per cell [0:ncell-1, 27]
    integer, dimension(:),   allocatable :: head  ! First atom index in each cell [0:ncell-1]
    integer, dimension(:),   allocatable :: list  ! Next atom index in linked list [1:natoms]
    integer :: ncell                   ! Total number of spatial cells = maxi * maxj * maxk
    integer :: nn                      ! Number of neighboring cells checked (3^ndim = 27)
    integer :: maxi, maxj, maxk        ! Number of cells along x, y, z axes
    logical :: use_cell = .true.       ! Active flag: true = use link cells, false = direct sum
end module linkcell


!-------------------------------------------------------------------------------
! Module: interp
! Purpose: Cubic spline potential interpolation tables and Breeuwsma matrix.
!-------------------------------------------------------------------------------
module interp
    use set_precision
    implicit none

    real(wp), dimension(:,:), allocatable :: utab ! Tabulated pair potential values [ncut, nitmax]
    real(wp), parameter :: dr = 0.001d0           ! Distance grid resolution (0.001 Angstrom)
    real(wp), parameter :: idr = 1000.0d0         ! Inverse grid resolution (1 / dr)
    real(wp), dimension(:), allocatable :: rmin2  ! Squared minimum allowable distance per pair

    ! Paul Breeuwsma smooth cubic spline interpolation matrix:
    ! Coefficients ensure continuous function and derivative values across grid boundaries.
    real(wp), dimension(0:3,0:3) :: am = reshape((/ &
        -0.5d0,  1.0d0, -0.5d0, 0.0d0, &
         1.5d0, -2.5d0,  0.0d0, 1.0d0, &
        -1.5d0,  2.0d0,  0.5d0, 0.0d0, &
         0.5d0, -0.5d0,  0.0d0, 0.0d0 /), (/ 4, 4 /))
    integer :: ncut                               ! Number of grid points up to cutoff
end module interp


!-------------------------------------------------------------------------------
! Module: interfaces
! Purpose: Explicit interface signatures for core calculation routines.
!-------------------------------------------------------------------------------
module interfaces
    implicit none

    interface
        subroutine move_natoms(natom)
            use set_precision
            integer, intent(in) :: natom
        end subroutine move_natoms

        subroutine histograms(ensemble)
            character, intent(in) :: ensemble*3
        end subroutine histograms

        subroutine printout(inst)
            logical, intent(in) :: inst
        end subroutine printout

        function fpot_Morse(rx, nit)
            use set_precision
            real(wp) :: fpot_Morse, rx
            integer  :: nit
        end function fpot_Morse

        function fpot_elecMorse(rx, nit)
            use set_precision
            real(wp) :: fpot_elecMorse, rx
            integer  :: nit
        end function fpot_elecMorse

        function fpot_LJ(rx, nit)
            use set_precision
            real(wp) :: fpot_LJ, rx
            integer  :: nit
        end function fpot_LJ

        function fpot_elecLJ(rx, nit)
            use set_precision
            real(wp) :: fpot_elecLJ, rx
            integer  :: nit
        end function fpot_elecLJ

        function dist2(r)
            use set_precision
            use configuration, only : ndim
            real(wp) :: dist2
            real(wp), dimension(ndim) :: r
        end function dist2
    end interface

end module interfaces
