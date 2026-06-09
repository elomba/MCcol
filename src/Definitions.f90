module configuration
    !
    ! Define system configuration
    !
    use set_precision
    implicit none
    Integer, Parameter :: ndim=3
    Integer :: length
    Integer :: nsp, natoms
    Integer, Allocatable ::ntype(:), iatype(:), seed(:)
    Real(wp), Dimension(:), Allocatable :: q, qsp
    Real (wp), Allocatable :: r(:,:)
    Real (wp), Dimension(ndim) :: a, b, c, side, side2, r_unit
    Real (wp) :: v0
    Character, Dimension(:), Allocatable :: atoms*8
end module configuration


module potential
    !
    ! Define potential parameters, units, and physical constants.
    !
    Use set_precision
    Implicit None
    Character :: units*8
    Character, dimension(:), allocatable :: pot*4
    Integer, dimension(:), allocatable :: keyp
    Integer :: kmt, kmx, kmy, kmz, nitmax, kint
    Integer, Dimension(:,:), Allocatable :: itp
    Real (wp) :: dospix, dospiy, dospiz, rcpcut, rcpcut2
    Real(wp), Dimension(:,:), Allocatable :: aa, bb, cc, rc
    Real(wp), Dimension(:), Allocatable :: km2, ekm2, kr, rc2, al, bl, cl, bl2, qprod, ucut
    Complex(wp), Dimension(:), Allocatable :: rhokk, deltann, einx, einy, einz
    Complex(wp), Dimension(:,:), Allocatable :: eix, eiy, eiz
    !
    !  ctreV = e^2/(4pi epsilon_0) in eVxAmstrong, kb = Boltzmann's
    !   constant in ev/K, ev2k, conversion factor from eV to K
    !
    Real(wp) :: ctr
    Real(dkind), Parameter :: ctreV= 14.39964361d0, pi&
        &=3.141592653589793d0, kbeV=8.6173324d-5, ev2K=11604.58857702&
        &, ev2Kcal=23.06054195, kbKcal=0.001987191686
          ! Eva's change: parameters to convert pressure
    !      bar2eV: change units from bar to eV/Angstrom^3
    !      bar2K: change from bar to K/Angstrom^3
    !      bar2Kcal: change from bar to Kcal/mol/Å³
    Real (wp), Parameter :: bar2eV = 0.624150932d-6, bar2K =&
         & 0.72429715652d-2, bar2Kcal=1.4393258781725146d-5   
    Real (wp) :: rcut, rcut2, kappa,  qtotal
    Real (dkind) :: selfe=0.0d0, pi2
    Complex (wp), Parameter :: ii=(0.0d0,1.0d0)
    logical :: elect, pshift
    Integer :: fou_type
end module potential

Module rundata
    !
    ! Specific data for the run
    !
    Use configuration
    integer, parameter :: iosys=7 ! Unit for input system config
    integer, parameter :: iocfg=8 ! Unit for i/o configuration
    integer, parameter :: iotrj=9 ! Unit for outout trajectory file
    integer, parameter :: iorun=10 ! Unitr for run parameters
    integer, parameter :: ioth=11 ! Thermal averages output
    integer, parameter :: iothi=12 ! Instantaneous therm quantities
    integer, parameter :: igr=13   ! output g(r)
    integer, parameter :: ilong=20  ! length of data and results filenames
    logical :: restart
    Integer :: nb, nstep, istep, istep_ini=0, ntrial, naccept, nvaccept, naver, nequil, npgr, ntraj
    Real(wp) :: temp, kT, s_cput, e_cput
    Real(wp) :: rdmax(1:ndim)
    ! NpT ensemble (Eva's change)
    Real(wp) :: pres, vdmax
    character :: ensemble*3, initcf*3, scaling*5, stat*10&
         &="sequential"  
    Character(len=8), Parameter :: c_blue   = char(27)//'[1;34m'
    Character(len=8), Parameter :: c_cyan   = char(27)//'[1;36m'
    Character(len=8), Parameter :: c_green  = char(27)//'[1;32m'
    Character(len=8), Parameter :: c_yellow = char(27)//'[1;33m'
    Character(len=8), Parameter :: c_red    = char(27)//'[1;31m'
    Character(len=8), Parameter :: c_reset  = char(27)//'[0m'
End Module rundata

Module properties
    !
    !  System properties as evaluated during the run
    !
    Use set_precision
    use configuration, only : ndim
    Real(dkind) :: Etav, E_sav, E_lrav, E_vdwav
    Real(dkind) :: Etotal, E_sr, E_coulomb, E_fourier, Evdw
    Real(dkind) :: virial
    Real(dkind) :: Vol_av, side_av(1:ndim)
    Real(wp) :: deltagr
     ! NpT ensemble (Eva's change)
    Real(wp), Dimension(:,:,:), Allocatable :: gmix
    Integer, Dimension(:,:,:), Allocatable :: histomix
    Integer :: nmaxgr
End Module properties

Module linkcell
    !
    ! Shared components in the link cell method
    !
    Use set_precision
    Real (wp) :: cellx, celly, cellz
    Integer, Dimension(:,:), Allocatable:: neigh
    Integer, Dimension(:), Allocatable :: head, list
    Integer :: ncell, nn, maxi,maxj,maxk
    logical :: use_cell=.true.
End Module linkcell
  
module interp
    !
    ! Parameters and matrix to store interpolated pair interactions
    !
    Use set_precision
    real (wp), dimension(:,:), allocatable :: utab
    real (wp), parameter :: dr=0.001d0, idr=1000.0d0
    real (wp), dimension(:), allocatable :: rmin2
    ! Interpolation matrix
    real (wp), dimension(0:3,0:3) :: am=(/-0.5d0,1.0d0,-0.5d0,0.0d0&
         &,1.5d0,-2.5d0,0.0d0,1.0d0,-1.5d0,2.0d0,0.5d0,0.0d0,0.5d0,&
         &-0.5d0,0.0d0,0.0d0/)
    integer :: ncut
end module interp

module interfaces
    !
    ! Interfaces to routine that use input parameters.
    !
    interface
        Subroutine move_natoms(natom)
            use set_precision
            Integer, Intent(IN) :: natom
        End Subroutine move_natoms
        subroutine histograms(ensemble)
            character, intent(IN) :: ensemble*3
        end subroutine histograms
        subroutine printout(inst)
            logical, intent(IN) :: inst
        end subroutine printout
        function fpot_Morse(rx,nit)
            use set_precision
            real(wp) :: fpot_Morse, rx
            integer :: nit
        end function fpot_Morse
        function fpot_elecMorse(rx,nit)
            use set_precision
            real(wp) :: fpot_elecMorse, rx
            integer :: nit
        end function fpot_elecMorse
        function fpot_LJ(rx,nit)
            use set_precision
            real(wp) :: fpot_LJ, rx
            integer :: nit
        end function fpot_LJ
        function fpot_elecLJ(rx,nit)
            use set_precision
            real(wp) :: fpot_elecLJ, rx
            integer :: nit
        end function fpot_elecLJ
        Function dist2(r)
            Use set_precision
            Use configuration, Only : ndim
            Real(wp) :: dist2
            Real(wp), dimension(ndim) :: r
        End Function dist2
    end interface
  
end module interfaces
