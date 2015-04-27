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
    Character :: units*2
    Character, dimension(:), allocatable :: pot*4
    Integer, dimension(:), allocatable :: keyp
    Integer :: kmt, kmx, kmy, kmz, nitmax
    Integer, Dimension(:,:), Allocatable :: itp
    Real (wp) :: dospix, dospiy, dospiz, rcpcut, rcpcut2
    Real(wp), Dimension(:,:), Allocatable :: aa, bb, cc, rc
    Real(wp), Dimension(:), Allocatable :: km2, ekm2, kr, rc2, al, bl, cl, bl2, qprod
    Complex(wp), Dimension(:), Allocatable :: rhokk, deltann, einx, einy, einz
    Complex(wp), Dimension(:,:), Allocatable :: eix, eiy, eiz
    !
    !  ctreV = e^2/(4pi epsilon_0) in eVxAmstrong, kb = Boltzmann's
    !   constant in ev/K, ev2k, conversion factor from eV to K
    !
    Real(wp) :: ctr
    Real(wp), Parameter :: ctreV= 14.39964361d0, pi&
        &=3.141592653589793d0, kbeV=8.6173324d-5, ev2K=11604.58857702
          ! Eva's change: parameters to convert pressure
    !      bar2eV: change units from bar to eV/Angstrom^3
    !      bar2K: change from bar to K/Angstrom^3
    Real (wp), Parameter :: bar2eV = 0.624150932d-6, bar2K = 0.72429715652d-2
    Real (wp) :: rcut, rcut2, kappa, selfe=0, qtotal, pi2
    Complex (wp), Parameter :: ii=(0.0d0,1.0d0)
    logical :: elect
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
    integer, parameter :: ioh=14   ! histograms
    logical :: restart
    Integer :: nb, nstep, istep, istep_ini=0, ntrial, naccept, nvaccept, naver, nequil, npgr, ntraj
    Real(wp) :: temp, kT, s_cput, e_cput
    Real(wp) :: rdmax(1:ndim)
    ! NpT ensemble (Eva's change)
    Real(wp) :: pres, vdmax
    character :: ensemble*3, initcf*3, scaling*5, stat*10="sequential"
End Module rundata

Module properties
    !
    !  System properties as evaluated during the run
    !
    Use set_precision
    use configuration, only : ndim
    Integer, Parameter :: nde=50000, nrmax=5000
    Real(wp) :: Etotal, E_sr, E_coulomb, E_fourier, Eref, deltaEng, deltar, deltagr
    Real(wp) :: Etav, E_sav, E_lrav
     ! NpT ensemble (Eva's change)
    Real(wp) :: Vol_av, side_av(1:ndim)
    Real(wp) :: virial
    Real(wp), Dimension(-nde:nde) :: Ehisto
    Real(wp), Dimension(0:nrmax) :: rhisto
    Real(wp), Dimension(:,:,:), Allocatable :: gmix
    Integer, Dimension(:,:,:), Allocatable :: histomix
    Integer :: nmaxgr
End Module properties

Module linkcell
    !
    ! Shared components in the link cell method
    !
    Use set_precision
    Real (wp) :: cellx, celly, cellz, cellxo, cellyo, cellzo
    Integer, Dimension(:,:), Allocatable:: neigh, neigho
    Integer, Dimension(:), Allocatable :: head, list, heado, listo
    Integer :: ncell, ncellmax=0, nn, maxi,maxj,maxk, ncello, nno, maxio,maxjo,maxko
    logical :: use_cell=.true.
End Module linkcell
  
module interp
    !
    ! Parameters and matrix to store interpolated pair interactions
    !
    Use set_precision
    real (wp), dimension(:,:), allocatable :: utab
    real (wp), parameter :: dr=0.001d0
    real (wp), dimension(:), allocatable :: rmin2
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
        function fpot(rx,nit)
            use set_precision
            real(wp) :: fpot, rx
            integer :: nit
        end function fpot
        Function dist2(r)
            Use set_precision
            Use configuration, Only : ndim
            Real(wp) :: dist2
            Real(wp), dimension(ndim) :: r
        End Function dist2
    end interface
  
end module interfaces
