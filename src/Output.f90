module output
    use rundata, only : ioth, iothi, igr
    public :: printout
contains
    subroutine initout
        implicit none
        character :: date*8, time*10, name*20
        call date_and_time(date,time)
        call hostnm(name)
        Write(*,'(/20x," ***** Program gpMC ****** "//" Running &
         &on ",a20/" Date: ",a2,"/",a2,"/",a4," at ",a2,":",a2,":"&
         &,a2/)'        )name,date(7:8),date(5:7),date(1:4),time(1:2)&
            &,time(3:4),time(5:6)
    end subroutine initout

    Subroutine end_printout
        Use util, Only : cputime
        Use rundata, Only : s_cput, e_cput
        implicit none
        character :: date*8, time*10, name*20
        call date_and_time(date,time)
        e_cput = cputime()
        Write(*,'(/20x," ***** End Program gpMC ****** "//1x,&
         &" Date ",a2,"/",a2,"/",a4," at ",a2,":",a2,":"&
         &,a2/"  Elapsed CPU time",f15.2," s")'        )date(7:8),date(5:7),date(1:4),time(1:2)&
            &,time(3:4),time(5:6), e_cput-s_cput
    End Subroutine end_printout


    subroutine run_info
        use potential, only : units, elect
        use rundata, only : ensemble
        implicit none
        write(*,'(/" Simulating ",a3," ensemble. Energy units (",a8,")")')ensemble,units
        if (elect) then
            write(*,'(" Ewald electrostatics "/)')
        else
            write(*,'(" No electrostatics "/)')
        endif

    end subroutine run_info
    subroutine init_printout
        use potential, only : elect
        use rundata, only : ensemble, stat
        implicit none
        Open(ioth,file='thermoaver.dat',access=stat)
        Open(iothi,file='thermoins.dat',access=stat)
        if (ensemble == 'nvt') Then
            Write(ioth,1000)
1000        format("# No. moves  % accept.       <E_tot>         <E_sr&
                 &>        <E_vdw>      <E_coul>   "/"#",80("-"))
            if (elect) then
                Write(iothi,1010)
                write(*,1020)
1010            format(/" No. moves  % accept.        E_tot           E_sr           E_vdw         E_coul "/1x,100("-"))
1020            format("# No. moves  % accept.        E_tot           &
                     &E_sr           E_vdw         E_Four       E_self&
                     &        E_coul "/"#",120("-"))
            else
                Write(iothi,1030)
1030            format(/" No. moves  % accept.        E_tot           E_sr         E_vdw       "/1x,100("-"))
                Write(*,1030)
            endif
        elseif (ensemble == 'npt' ) Then
            Write(ioth,1040)
1040        format("# No. moves  % accept.   % accept. vol.   <E_tot>      <E_sr>       <E_vdw>        <Vol>      <Lx>     <Ly>    <Lz>"&
                /"#",120("-"))
            if (elect) then
                Write(iothi,1050)
                write(*,1060)
1050            format("# No. moves  % accept.  % acc. vol.     E_tot           E_sr        E_vdw        E_coul     Vol"&
                    /"#",120("-"))
1060            format("# No. moves  % accept.   % acc. vol.     E_tot           E_sr           E_vdw        E_Four       E_self     Vol"&
                    &/"#",120("-"))
            else
                Write(iothi,1070)
                Write(*,1080)
1070            format("# No. moves  % accept.   % acc. vol.    E_tot          E_sr       E_vdw       Vol     Lx      Ly    Lz"/"#",120("-"))
1080            format("# No. moves  % accept.    % acc. vol.    E_tot           E_sr         E_vdw       Vol"/"#",&
                    &100("-"))
            endif
        endif
    end subroutine init_printout

    Subroutine Printout(inst)
        use set_precision
        use properties, only : e_fourier, E_sr,  Etotal, Evdw,&
             & E_vdwav, E_sav, E_coulomb, Etav, side_av, vol_av 
        use potential, only : selfe, elect
        use configuration, only : natoms, v0, side
        use rundata, only : ntrial, naccept, naver, nvaccept, ensemble, kT
        Implicit None
        logical, intent(IN) :: inst
        Real(wp), Save ::  Etavq=0
        Real(wp)      :: prob_vol
        Integer, save :: npeq = 0
        Integer       :: ncycles
        Etotal = E_sr+E_Fourier+selfe
        ncycles = ntrial/natoms
        prob_vol = 100*Dble(nvaccept)/dble(ncycles)

        If (inst) Then
            If (ensemble == 'nvt' ) then
                Write(ioth,'(i9,3x,f8.4,3x,4g15.7)')ntrial/natoms, 100*Dble(naccept)/Dble(ntrial),&
                &Etav*kT/naver, E_sav*kT/naver, E_vdwav*kT/naver, (Etav - E_vdwav)*kT/naver 
            Elseif (ensemble == 'npt' ) then
                Write(ioth,'(i9,3x,2f9.4,3x,7g15.7)')ntrial/natoms, 100*Dble(naccept)/Dble(ntrial),&
                &prob_vol,Etav*kT/naver, E_sav*kT/naver, E_vdwav*kT/naver, vol_av/naver, side_av(:)/naver
            Endif

        Else
            !
            ! Set energy reference to average value during equilibration.
            !
            npeq = npeq+1
            Etavq = Etavq+Etotal
        Endif
        If (ensemble == 'nvt' ) then

            if (elect) then
                Write(iothi,'(i9,3x,f8.4,3x,5g15.7)') ntrial/natoms, 100*Dble(naccept)/dble(ntrial),&
                &Etotal*kT, E_sr*kT, Evdw*kT, (E_fourier+selfe+E_sr-Evdw)*kt
                Write(*,'(i9,3x,f8.4,3x,7g15.7)') ntrial/natoms, 100*Dble(naccept)/Dble(ntrial),&
                &Etotal*kT, E_sr*kT, Evdw*kT, E_Fourier*kT, selfe*kT,&
                & (E_fourier+selfe+E_sr-Evdw)*kt
            else
                Write(iothi,'(i9,3x,f8.4,3x,5g15.7)') ntrial/natoms, 100*Dble(naccept)/dble(ntrial),Etotal*kT, E_sr*kT, (E_Fourier+selfe)*kT 
                Write(*,'(i9,3x,f8.4,3x,7g15.7)') ntrial/natoms, 100*Dble(naccept)/dble(ntrial),&
                &Etotal*kT, E_sr*kT, Evdw*kT
            endif
        ElseIf (ensemble == 'npt' ) then
            if (elect) then
                Write(iothi,'(i9,3x,2f9.4,3x,10g15.7)') ncycles, 100*Dble(naccept)/dble(ntrial), &
                &prob_vol, Etotal*kT, E_sr*kT, Evdw*kT, (E_fourier+selfe+E_sr-Evdw)*kt, v0, side(:)
                Write(*,'(i9,3x,2f9.4,3x,11g15.7)') ncycles, 100*Dble(naccept)/Dble(ntrial), prob_vol, &
                &Etotal*kT, E_sr*kT, Evdw*kT, E_Fourier*kT, selfe*kT,&
                &  v0
            else
                Write(iothi,'(i9,3x,2f9.4,3x,10g15.7)') ncycles, 100*Dble(naccept)/dble(ntrial), &
                &prob_vol, Etotal*kT, E_sr*kT, v0, side(:)
                Write(*,'(i9,3x,2f9.4,3x,10g15.7)') ncycles, 100*Dble(naccept)/dble(ntrial), &
                &prob_vol, Etotal*kT, E_sr*kT, Evdw*kT, v0
            endif
        Endif

    End Subroutine Printout

    Subroutine printgr
        Use set_precision
        Use configuration, only : nsp
        Use properties, Only : nmaxgr, deltagr, gmix
        Implicit None
        real(wp) :: ri
        Integer :: i, j, l
        Open(igr,file="gmix.dat")
        Write(igr,'("#      r(Å)   ")',advance='no')
        Do i=1, nsp
            Do j=i,nsp
                Write(igr,'(9x,"g(",2i1,")",1x)',advance="no")i,j
            End Do
        End Do
        Write(igr,'(/"#",100("-"))')
        Do i=1, nmaxgr
            ri = i*deltagr
            Write(igr,'(18f15.7)')i*deltagr,((gmix(i,j,l),l=j,nsp),j=1,nsp)
        End Do
        close(igr)
    End Subroutine printgr


    Subroutine print_ener
        !    new subroutine: added by Eva
        Use set_precision
        Use potential, Only : selfe, elect, units, rcpcut,qtotal
        Use properties, Only: E_sr,E_Fourier, Evdw
        Use rundata, Only: kT
        Write(*,'(80("-")/" vdW Energy=",g15.7,1x,a8)')Evdw*kT,units
        If (elect) Then
            write(*, '(" K cutoff",f10.5," 1/A")')rcpcut
            write(*, '(" Coulomb energy (self energy) ",g15.7,1x,a8)') selfe*kT, units
            write(*, '(" Coulomb energy (Fourier term)",g15.7,1x,a8)') E_Fourier*kT, units
            write(*, '(" Coulomb energy ",g15.7,1x,a8)') (selfe&
                 &+E_Fourier+E_sr-Evdw)*kT,units
            write(*, '(" E_sr (SR Coul+vdW) =",g15.7,1x,a8)') E_sr*kT, units
            write(*, '(" Deviation from charge neutrality =",g15.7)') qtotal
        EndIf
        Write(*, '(" Etotal =",g15.7,1x,a8/80("-"))') (selfe+E_Fourier+E_sr)*kT,units

    End Subroutine print_ener

end module output

