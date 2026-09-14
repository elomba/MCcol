!===============================================================================
! Module: Init
!
! Purpose:
!   Orchestrates simulation initialization, parameter parsing from input files,
!   force field setup, Ewald reciprocal grid construction, and pair potential
!   interpolation tabulations.
!
! Routines:
!   - Init_conf    : Reads system.dat, loads configuration, sets box metrics.
!   - Init_rundata : Reads runMC.dat, parses ensemble, steps, RNG seeds.
!   - read_potpars : Parses pair potential parameters (Morse or LJ).
!   - Init_pot     : Converts units, defines cutoffs, shifts potentials.
!   - Init_selfe   : Evaluates Ewald self-energy, allocates k-space arrays.
!   - Init_fourier : Constructs 3D reciprocal wavevectors and Ewald weight factors.
!   - Init_interp  : Precalculates spline interpolation table utab on radial grid.
!===============================================================================
Module Init
    Use set_precision
    Use configuration
    Use potential
    Use properties
    Use rundata
contains

    !---------------------------------------------------------------------------
    ! Subroutine: Init_conf
    !
    ! Purpose:
    !   Reads system.dat to obtain species count, atom counts, charges, units,
    !   and forcefield parameters. Loads coordinates via dlplmp_readconf,
    !   computes box volume, and normalizes particle coordinates to [-0.5, 0.5).
    !---------------------------------------------------------------------------
    Subroutine Init_conf
        use readconf, only : dlplmp_readconf
        Implicit None
        Integer :: keytrj, imcon, iatm, i, j, nit
        Real(wp) :: dumx, dumy, dumz, qsp2
        Open (iosys,file='system.dat')
        read(iosys,*) restart
        if (restart) then
            ! Load binary dump file for restart
            call load
            return
        endif
        read(iosys,*) initcf
        Read(iosys,*) nsp, natoms
        
        nitmax = (nsp*nsp+nsp)/2
        Allocate(ntype(nsp),atoms(nsp),qsp(nsp),q(natoms),&
            & qprod(nitmax))
        q(:) = 0.0_wp
        do i=1, nsp
           Read(iosys,*) j, atoms(j),qsp(j)
        end do
        Read(iosys,*) units

        Allocate(R(ndim,1:natoms),iatype(natoms))
        nit = 1
        do i=1, nsp
            do j=i, nsp
                qprod(nit) = qsp(i)*qsp(j)
                nit = nit+1
            end do
        end do
        !
        !  Input potential parameters
        !
        call read_potpars
        !
        ! Read in Ewald convergence parameter and no. of k vectors in each
        ! direction. Warning, place the z-axis along the longest unit cell direction.
        !
        qsp2 = dot_product(qsp(1:nsp),qsp(1:nsp))
        if (qsp2 > 1.d-6) then
           read(iosys,*) elect, pshift, fou_type
           if (elect) Read(iosys,*) kappa, kmx, kmy, kmz
        else
            elect = .false.
        endif
        !
        ! When reading or generating the initial configuration
        ! the numbers of particles of each type (ntype) are defined
        !
        if (initcf == "dlp" .or. initcf == "lmp") then
            ! Read in DLPOLY CONFIG or LAMMPS data.atoms file
            call dlplmp_readconf
        else
            print *, "*** Input error ",initcf," not supported as input configuration"
            stop
        endif
        !
        ! This is valid only for orthorombic cells
        !
        side(1) = a(1)
        side(2) = b(2)
        side(3) = c(3)
        side2(:) = side(:)/2
        !
        ! Internal length units are defined in terms of the sides
        ! of the simulation box
        !
        r_unit(1) = Sqrt(Dot_product(a,a))
        r_unit(2) = Sqrt(Dot_product(b,b))
        r_unit(3) = Sqrt(Dot_product(c,c))
        side(:) = r_unit(:)
        ! Simulation box volume
        v0 = a(1)*b(2)*c(3)+a(2)*b(3)*c(1)+a(3)*b(1)*c(2)-a(3)*b(2)*c(1)&
            &-a(2)*b(1)*c(3)-a(1)*b(3)*c(2)
        !
        ! Rescale atomic coordinates to box length units [-0.5, 0.5)
        !
        Forall (i=1:natoms) R(1:ndim,i) = R(1:ndim,i)/r_unit(1:ndim)
    End Subroutine Init_conf


    !---------------------------------------------------------------------------
    ! Subroutine: Init_rundata
    !
    ! Purpose:
    !   Parses runMC.dat for ensemble ('nvt' or 'npt'), production/equilibration
    !   sweep counts, averaging intervals, maximum trial displacements,
    !   temperature, pressure, RDF binning, and seeds the random number generator.
    !---------------------------------------------------------------------------
    Subroutine Init_rundata
        Implicit None
        Open(iorun,file='runMC.dat')
        ! Read in name of results directory
        read(iorun,*) ensemble
        If (ensemble .Ne. "nvt" .and. ensemble.Ne. "npt") Then
            Print *, " *** Input error:", ensemble," not implemented .."
            stop
        End If
        ! No. of steps, no. of squilibration steps, average every nb, print
        ! gr ever ngr, no. of step to dump trajectory (0 no dumps)
        !
        Read(iorun,*) nstep, nequil, nb, npgr, ntraj
        ! Max. displacement and Temperature
        If (ensemble == 'nvt') Then
            Read(iorun,*) rdmax(1:ndim)
            Read(iorun,*) temp
            Read(iorun,*) deltagr
        else if(ensemble == 'npt') Then
            ! Max displacement and volumen change, scaling, temperature and pressure
            Read(iorun,*) rdmax(1:ndim), vdmax, scaling
            If (Trim(Adjustl(scaling)) .Ne. "ortho".and. Trim(Adjustl(scaling)).Ne."isotr") Then
               Print *, " *** Input error:", scaling," not implemented .."
               stop
            Endif
            Read(iorun,*) temp,pres
            Read(iorun,*) deltagr
            
        Endif
        !
        ! Read and set seed for random number generator
        if (.not. restart) then
            call random_seed
            call random_seed(size=length)
!            print *, "** Seed length in words:",length
            allocate(seed(1:length))
            Read(iorun,*) seed(1:length)
            call random_seed(put=seed(1:length))
            !
            ! Set counters and accumulators to 0
            !
            naccept = 0
            ntrial = 0
            naver = 0
            Etotal = 0
            Etav = 0
            E_sav = 0
            E_vdwav = 0
        endif
        !
        !
        !
        if (ntraj .ne. 0) then
            open(iotrj,file="gpMC.lammpstrj")
        endif
    End Subroutine Init_rundata


    !---------------------------------------------------------------------------
    ! Subroutine: read_potpars
    !
    ! Purpose:
    !   Reads interaction potential types and parameter matrices from system.dat.
    !   Supports Morse ('mors', keyp=1) and Lennard-Jones ('lj', keyp=2).
    !   Populates symmetric indexing matrix itp(i, j) -> nit.
    !---------------------------------------------------------------------------
    subroutine read_potpars
        implicit none
        integer :: nit, i, j, k, l

        Allocate(aa(nsp,nsp),cc(nsp,nsp),bb(nsp,nsp),rc(nsp&
            &,nsp),itp(nsp,nsp),rc2(nitmax),pot(nitmax)&
            &,keyp(nitmax),al(nitmax),bl(nitmax),&
            & cl(nitmax),bl2(nitmax),ucut(nitmax))
        ucut(:)=0.0d0
        nit = 1

        ! Populate symmetric mapping matrix itp: (species_i, species_j) -> pair_index
        do i = 1, nsp
           do j = i, nsp
              itp(i,j) = nit
              itp(j,i) = nit
              nit = nit+1
           enddo
        enddo

        Do k= 1,nsp
            Do l= k, nsp
               Read(iosys,*) i, j, pot(itp(i,j))
               nit = itp(i,j)
               If (Trim((Adjustl(pot(nit)))) == Trim(Adjustl("mors")) ) Then
                  keyp(nit) = 1
                  Read(iosys,*) aa(i,j),bb(i,j),cc(i,j),rc(i,j)
                  cl(nit) = cc(i,j)
               Else If (Trim((Adjustl(pot(nit)))) == Trim(Adjustl("lj")) ) Then
                  keyp(nit) = 2
                  Read(iosys,*) aa(i,j),bb(i,j),rc(i,j)
                  bl2(nit) = bb(i,j)**2
               Else
                  Print *, "*** Input error:", pot(nit),"  not implemented yet.."
                  stop
               End If
               al(nit) = aa(i,j)
               bl(nit) = bb(i,j)
               rc(j,i) = rc(i,j)
               rc2(nit) = rc(i,j)**2
            End Do
         End Do
         if (nitmax.gt.1) then
            if (keyp(1) .ne. keyp(2)) then
               print *, " Error : only one type of interaction implemented, Morse or LJ"
               stop
            endif
         endif
         kint = keyp(1)
    end subroutine read_potpars


    !---------------------------------------------------------------------------
    ! Subroutine: Init_pot
    !
    ! Purpose:
    !   Initializes physical constants and unit conversion factors for the chosen
    !   energy unit system ('eV', 'K', or 'kcal/mol').
    !   Determines overall cutoff rcut, computes potential shifts ucut if pshift=.true.,
    !   and initializes Ewald and RDF data structures.
    !---------------------------------------------------------------------------
    Subroutine Init_pot
      Implicit None
      real(wp), external :: fpot_LJ, fpot_Morse
      logical :: old_elect
      Integer :: i,j, nit

      ! Adjust thermal energy kT and conversion factors based on energy units
      If (Trim(Adjustl(units)) == Trim(Adjustl("eV"))) Then
         kT = kbev*temp
         ctr = ctreV/kT
         pres = pres*bar2eV/kT
      Else If (Trim(Adjustl(units)) == Trim(Adjustl("K"))) Then
         kT = temp
         ctr = ctreV*ev2k/kT
         pres = pres*bar2k/kT
      Else if (Trim(Adjustl(units)) == Trim(Adjustl("kcal/mol"))) Then
         kT = kbKcal*temp
         ctr = ctreV*ev2Kcal/kT
         pres = pres*bar2Kcal/kT
      Else
         Print *, " *** Input error:",units," not implemented as energy unit"
         Stop
      End If

      ! Scale well-depth parameter al with thermal energy 1/kT
      if (kint == 2) then
         al(1:nitmax) = 4.0d0*al(1:nitmax)/kT
      else
         al(1:nitmax) = al(1:nitmax)/kT
      endif

      ! Set global cutoff to maximum of pairwise cutoffs
      rcut = Maxval(rc(:,:))
      rcut2 = rcut**2

      ! Compute potential shift at cutoff if requested
      if (pshift) then
         old_elect = elect
         elect = .false.
         nit = 1
         do i=1,nsp
            do j=i,nsp
               if (kint == 1) then
                  ucut(nit) = fpot_Morse(rc(i,j)**2,nit)
               else
                  ucut(nit) = fpot_LJ(rc(i,j)**2,nit)
               endif
               nit=nit+1
            enddo
         enddo
         elect = old_elect
      endif

      ! Initialize Ewald electrostatics
      if (elect) then
         call init_selfe
         call init_fourier
      Endif

      ! Allocate RDF bins if fresh simulation run
      if (.not. restart) then
         nmaxgr = Nint(rcut/deltagr)
         Allocate(histomix(nmaxgr,nsp,nsp),gmix(nmaxgr,nsp,nsp))
         histomix(:,:,:) = 0
      endif
    End Subroutine Init_pot


    !---------------------------------------------------------------------------
    ! Subroutine: Init_selfe
    !
    ! Purpose:
    !   Calculates Ewald electrostatic self-energy:
    !     E_self = - (kappa / sqrt(pi)) * [e^2 / (4*pi*epsilon_0*kT)] * sum_i q_i^2
    !   Allocates reciprocal space arrays (eix, eiy, eiz, km2, ekm2, rhokk).
    !---------------------------------------------------------------------------
    subroutine Init_selfe
        implicit none
        Integer :: i

        selfe = 0.0d0
        qtotal = 0.0d0
        ! Initialize individual atomic charges from species charges
        Do i = 1, natoms
           if (abs(q(i)) < 1.0d-12) q(i) = qsp(iatype(i))
        End Do
        selfe = Dot_product(q(1:natoms),q(1:natoms))
        qtotal = Sum(q(1:natoms))
        selfe = - ctr*(kappa/Sqrt(pi))*selfe

        kmt = (2*kmx+1)*(2*kmy+1)*(kmz+1)
        Allocate(eix(1:natoms,-kmx:kmx),eiy(1:natoms,-kmy:kmy),eiz(1:natoms,0:kmz))
        Allocate(einx(-kmx:kmx),einy(-kmy:kmy),einz(0:kmz))
        Allocate(kr(ndim),km2(0:kmt),ekm2(kmt),rhokk(kmt),deltann(kmt))
        Write(*, '("** Charged system: Init Fourier terms with Ewald parameters:",f10.5,3i4)') kappa, kmx, kmy, kmz
    end subroutine Init_selfe


    !---------------------------------------------------------------------------
    ! Subroutine: Init_fourier
    !
    ! Purpose:
    !   Pre-calculates 3D reciprocal wavevectors k = 2*pi*(kx/Lx, ky/Ly, kz/Lz),
    !   their squared magnitudes km2, and Ewald reciprocal weights:
    !     ekm2 = exp( -k^2 / (4*kappa^2) ) / k^2
    !   Exploits inversion symmetry to store only non-redundant k-vectors.
    !---------------------------------------------------------------------------
    subroutine Init_fourier
        implicit none
        Integer :: i, ind, kx, ky, kz

        pi2 = 2.0d0 * pi
        rhokk(:) = 0.0d0
        dospix = 2.0d0*pi/side(1)
        dospiy = 2.0d0*pi/side(2)
        dospiz = 2.0d0*pi/side(3)

        ! Determine reciprocal cutoff
        rcpcut = 1.05d0*Min(dospix*kmx, dospiy*kmy, dospiz*kmz)
        rcpcut2 = rcpcut**2
        ind = 1

        ! Symmetry partition 1: kz=0, ky=0, kx=1..kmx
        Do kx = 1, kmx
            kr(1) = dospix*kx
            km2(ind) = kr(1)*kr(1)
            ekm2(ind) = Exp(-km2(ind)/(4.0d0*kappa**2))/km2(ind)
            ind = ind+1
        End Do
        ! Symmetry partition 2: kz=0, ky=1..kmy, kx=-kmx..kmx
        Do ky = 1, kmy
            kr(2) = dospiy*ky
            Do kx = -kmx, kmx
                kr(1) = dospix*kx
                km2(ind) = Dot_product(kr(1:2),kr(1:2))
                ekm2(ind) = Exp(-km2(ind)/(4.0d0*kappa**2))/km2(ind)
                ind = ind+1
            End Do
        End Do
        ! Symmetry partition 3: kz=1..kmz, ky=-kmy..kmy, kx=-kmx..kmx
        Do kz = 1, kmz
            kr(3) = dospiz*kz
            Do ky = -kmy, kmy
                kr(2) = dospiy*ky
                Do kx = -kmx, kmx
                    kr(1) = dospix*kx
                    km2(ind) = Dot_product(kr(:),kr(:))
                    ekm2(ind) = Exp(-km2(ind)/(4.0d0*kappa**2))/km2(ind)
                    ind = ind+1
                End Do
            End Do
        End Do

        ! Zero-frequency initial phase components
        eix(1:natoms,0) = (1.0d0, 0.0d0)
        eiy(1:natoms,0) = (1.0d0, 0.0d0)
        eiz(1:natoms,0) = (1.0d0, 0.0d0)
        einx(0) = 1.0d0
        einy(0) = 1.0d0
        einz(0) = 1.0d0
    end subroutine Init_fourier


    !---------------------------------------------------------------------------
    ! Subroutine: Init_interp
    !
    ! Purpose:
    !   Tabulates short-range pair potentials on a fine 1D radial grid with
    !   spacing dr = 0.001 Angstrom up to the cutoff radius.
    !   Tabulated arrays are used by Paul Breeuwsma cubic spline interpolation.
    !
    ! Arguments:
    !   f (external function) : Pair potential evaluation function.
    !---------------------------------------------------------------------------
    subroutine Init_interp(f)
        use interp, only : utab, dr, rmin2, ncut
        implicit none
        integer :: iti, itj, i, j, k, nit
        real (wp) :: upot, rr, upmax=80.0d0
        real(wp), external :: f

        ncut = nint((rcut+0.5d0)/dr)
        allocate(utab(ncut,nitmax),rmin2(nitmax))

        nit = 1
        do i = 1, nsp
            do j = i, nsp
                upot = 0.0d0
                k = ncut
                do while (k >= 1 .and. upot/kT < upmax)
                   rr = (k*dr)**2
                   upot = f(rr,nit)
                   utab(k,nit) = upot
                   k = k - 1
                End Do
                rmin2(nit) = rr
                utab(1:k,nit) = utab(k+1,nit)
                nit = nit+1
            end do
         end do
    end subroutine Init_interp

End Module Init
