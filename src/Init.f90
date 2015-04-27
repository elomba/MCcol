Module Init
    !
    ! This module contains initialization and input routines
    !
    Use set_precision
    Use configuration
    Use potential
    Use properties
    Use rundata
contains
    Subroutine Init_conf
        !
        ! Read in system configuration and initialize storage.
        ! Initialize particle configuration
        !
        use readconf, only : dlp_readconf
        use linkcell, only : use_cell, list, listo
        !
        !   Unit cell is tetragonal !!
        !
        Implicit None
        Integer :: keytrj, imcon, iatm, i, j, nit
        Real(wp) :: dumx, dumy, dumz, qsp2
        Open (iosys,file='data/system.dat')
        read(iosys,*) restart
        if (restart) then
            ! Load dump file
            call load
            return
        endif
        read(iosys,*) initcf
        Read(iosys,*) nsp, natoms
        nitmax = (nsp*nsp+nsp)/2
        Allocate(ntype(nsp),atoms(nsp),qsp(nsp),q(natoms),&
            & qprod(nitmax))
        do i=1, nsp
            Read(iosys,*) atoms(i),qsp(i)
        enddo
        Read(iosys,*) units

        Allocate(r(1:natoms,ndim),iatype(natoms))
        if (use_cell) Allocate(list(natoms),listo(natoms))
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
            Read(iosys,*) kappa, kmx, kmy, kmz
            elect = .true.
        else
            elect = .false.
        endif
        !
        ! When reading or generating the initial configuration
        ! the numbers of particles of each type (ntype) are defined
        !
        if (initcf == "dlp") then
            ! Read in DLPOLY CONFIG file with particle positions
            call dlp_readconf
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
        !
        ! Rescale atomic coordinates to box length units
        !
        Forall (i=1:natoms) R(i,1:ndim) = R(i,1:ndim)/r_unit(1:ndim)
    End Subroutine Init_conf


    Subroutine Init_rundata
        !
        ! Read data specific for the run
        !
        Implicit None
        Open(iorun,file='data/runMC.dat')
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
            Read(iorun,*) deltaEng, deltagr
        else if(ensemble == 'npt') Then
            ! Max displacement and volumen change, scaling, temperature and pressure
            Read(iorun,*) rdmax(1:ndim), vdmax, scaling
            If (Trim(Adjustl(scaling)) .Ne. "ortho".and. Trim(Adjustl(scaling)).Ne."isotr") Then
                Print *, " *** Input error:", scaling," not implemented .."
                stop
            Endif

            Read(iorun,*) temp,pres
            Read(iorun,*) deltar, deltagr
        Endif
        !
        ! Adjust T and conversion factors depending on units used.
        !
        If (units == "eV") Then
            ctr = ctreV
            kT = kbev*temp
            pres = pres*bar2eV

        Else If (Trim(Adjustl(units)) == Trim(Adjustl("K"))) Then
            !
            ! energy units are K
            !
            pres = pres*bar2k

            ctr = ctreV*ev2k
            kT = temp
        Else
            Print *, " *** Input error:",units," not implemented as energy un&
          &it"            
            Stop
        End If
        !
        !
        if (.not. restart) then
            call random_seed
            call random_seed(size=length)
            allocate(seed(1:length))
            !
            ! Commented out, useful to run parallel runs with different seeds
            !!$            print *, "** Seed length in words:",length
            !!$            Read(iorun,*) seed(1:length)
            !!$            call random_seed(put=seed(1:length))
            !
            ! Set counters and accumulators to 0
            !
            naccept = 0
            ntrial = 0
            naver = 0
            Etotal = 0
            Etav = 0
            E_sav = 0
            Ehisto(:) = 0.0d0
            rhisto(:) = 0.0d0
        endif
        !
        !
        !
        if (ntraj .ne. 0) then
            open(iotrj,file="results/traj.xyz",status=stat)
        endif
    End Subroutine Init_rundata


    subroutine read_potpars
        implicit none
        integer :: nit, i, j
        !
        !  Read parameters for interaction potential
        !
        Allocate(aa(nsp,nsp),cc(nsp,nsp),bb(nsp,nsp),rc(nsp&
            &,nsp),itp(nsp,nsp),rc2(nitmax),pot(nitmax)&
            &,keyp(nitmax),al(nitmax),bl(nitmax),&
            & cl(nitmax),bl2(nitmax))
        nit = 1
        !
        ! Morse (keyp=1) and LJ (keyp=2) implemented
        !
        Do i= 1,nsp
            Do j= i, nsp
                read(iosys,*) pot(nit)
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
                !
                ! itp controls transforms (i,j) notation for the
                ! interaction into a vector type
                !
                itp(i,j) = nit
                itp(j,i) = nit
                nit = nit+1
            End Do
        End Do
    end subroutine read_potpars

    Subroutine Init_pot
        !
        !  Initialize potential parameters
        !
        Implicit None
        Integer :: i,j, nit

        ! Set global cutoff to the maximum of the site-site cutoffs
        rcut = Maxval(rc(:,:))
        rcut2 = rcut**2
        !
        ! Initialized Ewald method
        !
        if (elect) then
            call init_selfe
            call init_fourier
        Endif
        if (.not. restart) then
            !
            ! With the cutoff defined, allocate arrays for g(r)
            !
            nmaxgr = Nint(rcut/deltagr)
            Allocate(histomix(nmaxgr,nsp,nsp),gmix(nmaxgr,nsp,nsp))
            histomix(:,:,:) = 0
        endif
    End Subroutine Init_pot

    subroutine Init_selfe
        implicit none
        Integer :: i
        !
        ! Self energy of Ewald contribution and allocate arrays for Ewald sums
        !
        selfe = 0
        qtotal = 0
        !
        ! Initialize q(1:natoms) with the atomic charges
        !
        Do i = 1, natoms
            q(i) = qsp(iatype(i))
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

    subroutine Init_fourier
        implicit none
        Integer :: i, ind, kx, ky, kz
        !
        ! Initialize Fourier components of Ewald contributions
        !
        pi2 = 2*pi
        rhokk(:) = 0.0d0
        dospix = 2*pi/side(1)
        dospiy = 2*pi/side(2)
        dospiz = 2*pi/side(3)
        !
        ! Determine the cutoff in Fourier space
        !
        rcpcut = 1.05*Min(dospix*kmx,dospiy*kmy,dospiz*kmz)
        rcpcut2 = rcpcut**2
        ind = 1
        !
        !  Initialize k-vectors for Fourier component of Ewald summation.
        !  Use k-space symmetry, i.e. kz = 0..kmz, if kz=0, ky=0..kmy, kx=-kmx..kmz
        !  and ky = -kmy..kmy otherwise, if kz=ky=0, kx=1..kmx and kx=-kmx..kmx
        !  otherwise
        !
        Do kx = 1, kmx
            kr(1) = dospix*kx
            km2(ind) =  kr(1)*kr(1)
            ekm2(ind) = Exp(-km2(ind)/(4*kappa**2))/km2(ind)
            ind = ind+1
        End Do
        Do ky = 1, kmy
            kr(2) = dospiy*ky
            Do kx = -kmx, kmx
                kr(1) = dospix*kx
                km2(ind) =  Dot_product(kr(1:2),kr(1:2))
                ekm2(ind) = Exp(-km2(ind)/(4*kappa**2))/km2(ind)
                ind = ind+1
            End Do
        End Do
        Do kz = 1, kmz
            kr(3) = dospiz*kz
            Do ky = -kmy, kmy
                kr(2) = dospiy*ky
                Do kx = -kmx, kmx
                    kr(1) = dospix*kx
                    km2(ind) =  Dot_product(kr(:),kr(:))
                    ekm2(ind) = Exp(-km2(ind)/(4*kappa**2))/km2(ind)
                    ind = ind+1
                End Do
            End Do
        End Do
        eix(1:natoms,0) = (1.0d0, 0.0d0)
        eiy(1:natoms,0) = (1.0d0, 0.0d0)
        eiz(1:natoms,0) = (1.0d0, 0.0d0)
        einx(0) = 1
        einy(0) = 1
        einy(0) = 1
    end subroutine Init_fourier



    subroutine Init_interp
        !
        !  Initialize interpolation tables for short range pair interactions
        !
        use interp, only : utab, dr, rmin2, ncut
        implicit none
        integer :: iti, itj, i, j, k, nit
        real (wp) :: upot, rr, upmax=80.0d0
        real(wp), external :: fpot
        ncut = nint ((rcut+0.5)/dr)
        allocate(utab(ncut,nitmax),rmin2(nitmax))
        !  allocate(itp(nsp,nsp))
        nit = 1
        do i = 1, nsp
            do j = i, nsp
                upot = 0.0
                k = ncut
                do while (k >= 1 .and. upot/kT < upmax)
                    rr = (k*dr)**2
                    upot = fpot(rr,nit)
                    utab(k,nit) = upot
                    k=k-1
                End Do
                rmin2(nit) = rr
                !        Print *, i,j,Sqrt(rmin2(nit)),upot,upot/kT
                utab(1:k,nit) = utab(k+1,nit)
                nit = nit+1
            end do
        end do
    end subroutine Init_interp

End Module Init
