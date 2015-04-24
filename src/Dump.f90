Subroutine Cierra(clean)
    Use rundata
    Use configuration
    Use properties
    Use potential
    Use linkcell
    Use interp
    Use Util, only : cputime
    Implicit None
    integer, intent(in) :: clean
    Character ::  Dia*10, Hora*10, Bhora*10, BDia*10, dfname*30
    !
    ! Dummy subrutine. Here the user should take appropriate actions for an orderly
    ! program exit in case of system crash with SIGTERM signals
    !
    Call Date_and_time(Dia,Hora)
    BHora = Hora(1:2) // ':' //  Hora(3:4) // ':' // Hora(5:6)
    BDia  = Dia(7:8) // '-' // Dia(5:6) // '-' // Dia(1:4)
    dfname = "results/dump"//dia(1:8)//hora//".dmp"
    call random_seed(size=length)
    call random_seed(get=seed(1:length))
    open (1000,file=dfname,form="unformatted")
    write(1000) seed,dia,hora, istep
    write(1000) ensemble, scaling, units
    write(1000) natoms, nitmax, nmaxgr
    write(1000) nsp, nb, nstep, ntrial, naccept, nvaccept, naver, nequil, npgr, ntraj
    write(1000) r, q, qsp, iatype, ntype, atoms, itp
    write(1000) aa, bb, cc, rc
    write(1000) dospix, dospiy, dospiz, rcpcut, rcpcut2
    write(1000) kmx, kmy, kmz, kmt
    write(1000) rc2, al, bl, cl, bl2
    write(1000) elect, ctr
    if (elect) then
        write(1000) km2, ekm2, kr
        write(1000) rhokk, deltann, einx, einy, einz, eix, eiy, eiz
    endif
    write(1000) rcut, rcut2, kappa, selfe, qtotal, pi2
    write(1000) Etotal, E_sr, E_coulomb, E_fourier, Eref, deltaEng, deltar, deltagr
    write(1000) Etav, E_sav, E_lrav
    write(1000) temp, kT, pres, deltaEng, deltagr, rdmax, vdmax
    write(1000) Vol_av, side_av(1:ndim), r_unit, a, b, c, side, side2, v0
    write(1000) histomix, gmix
    write(1000) ncut
    write(1000) utab, rmin2
    if (ensemble == "nvt") then
        write(1000) Ehisto
    else if (ensemble == "npt") then
        write(1000) rhisto
    endif
    close(1000)
    if (clean == 0) then
        Write(*,1000)BHora, Bdia, (cputime()-s_cput)/60
    else
        write(*,1100)BHora, Bdia, (cputime()-s_cput)/60
1000    Format("   ##################    WARNING      ##############################"/ &
            & "   #                                                               #"/ &
            & "   #       Abnormal exit  of Program gpMC                          #"/ &
            & "   #",19x,"  Time : ",a10,4x," Date : ",a10,"   #"/,  &
            & "   #",19x,"  CPU usage =",f15.3," mins.          #"/ &
            & "   #                                                               #"/ &
            & "   #################################################################")
1100    Format("   #################################################################"/ &
            & "   #                                                               #"/ &
            & "   #    dumping restart file Program gpMC                          #"/ &
            & "   #",19x,"  Time : ",a10,4x," Date : ",a10,"   #"/,  &
            & "   #",19x,"  CPU usage =",f15.3," mins.          #"/ &
            & "   #                                                               #"/ &
            & "   #################################################################")
    endif
    Stop
End Subroutine Cierra

Subroutine Load
    Use rundata
    Use configuration
    Use properties
    Use linkcell
    Use potential
    Use interp
    Use Util, only : cputime
    Implicit None
    Character ::  Dia*10, Hora*10, Bhora*10, BDia*10, dfname*30
    open (1000,file="results/restart.dmp",form="unformatted")
    call random_seed
    call random_seed(size=length)
    allocate(seed(1:length))
    Read(1000) seed(1:length), dia, hora, istep_ini
    call random_seed(put=seed(1:length))
    Read(1000) ensemble, scaling, units
    Read(1000) natoms, nitmax, nmaxgr
    Read(1000) nsp, nb, nstep, ntrial, naccept, nvaccept, naver, nequil, npgr, ntraj
    Allocate(r(1:natoms,ndim),iatype(natoms))
    Allocate(ntype(nsp),atoms(nsp),qsp(nsp),q(natoms),qprod(nitmax))
    Allocate(aa(nsp,nsp),cc(nsp,nsp),bb(nsp,nsp),rc(nsp&
        &,nsp),itp(nsp,nsp),rc2(nitmax),pot(nitmax)&
        &,keyp(nitmax),al(nitmax),bl(nitmax),&
        & cl(nitmax),bl2(nitmax))
    Read(1000) r, q, qsp, iatype, ntype, atoms, itp
    Read(1000) aa, bb, cc, rc
    Read(1000) dospix, dospiy, dospiz, rcpcut, rcpcut2
    read(1000) kmx, kmy, kmz, kmt
    read(1000) rc2, al, bl, cl, bl2
    read(1000) elect, ctr
    if (elect) then
        Allocate(eix(1:natoms,-kmx:kmx),eiy(1:natoms,-kmy:kmy),eiz(1:natoms,0:kmz))
        Allocate(einx(-kmx:kmx),einy(-kmy:kmy),einz(0:kmz))
        Allocate(kr(ndim),km2(0:kmt),ekm2(kmt),rhokk(kmt),deltann(kmt))
        Read(1000) km2, ekm2, kr, rc2, al, bl, cl, bl2
        Read(1000) rhokk, deltann, einx, einy, einz, eix, eiy, eiz
    endif
    Read(1000) rcut, rcut2, kappa, selfe, qtotal, pi2
    Read(1000) Etotal, E_sr, E_coulomb, E_fourier, Eref, deltaEng, deltar, deltagr
    Read(1000) Etav, E_sav, E_lrav
    Read(1000) temp, kT, pres, deltaEng, deltagr, rdmax, vdmax
    Read(1000) Vol_av, side_av(1:ndim), r_unit, a, b, c, side, side2, v0
    Allocate(histomix(nmaxgr,nsp,nsp),gmix(nmaxgr,nsp,nsp))
    Read(1000) histomix, gmix
    read(1000) ncut
    allocate(utab(ncut,nitmax),rmin2(nitmax))
    Read(1000) utab, rmin2
    if (ensemble == "nvt") then
        Read(1000) Ehisto
    else if (ensemble == "npt") then
        Read(1000) rhisto
    endif
    close(1000)
    stat = "append"
    BHora = Hora(1:2) // ':' //  Hora(3:4) // ':' // Hora(5:6)
    BDia  = Dia(7:8) // '-' // Dia(5:6) // '-' // Dia(1:4)
    write(*,1100)BHora, Bdia
1100 Format("   #################################################################"/ &
          & "   #                                                               #"/ &
          & "   #    Restarting from program dump made at                       #"/ &
          & "   #",19x,"  Time : ",a10,4x," Date : ",a10,"   #"/,  &
          & "   #                                                               #"/ &
          & "   #################################################################")
end subroutine load
