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
        write(*,'(/" Simulating ",a3," ensemble. Energy units (",a2,")")')ensemble,units
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
        Open(ioth,file='results/thermoaver.dat',access=stat)
        Open(iothi,file='results/thermoins.dat',access=stat)
        if (ensemble == 'nvt') Then
            Write(ioth,'(/" No. moves  % accept.       <E_tot>         <E_sr>"/1x,80("-"))')
            if (elect) then
                Write(iothi,'(/" No. moves  % accept.        E_tot           E_sr         E_coul "/1x,80("-"))')
                write(*,'(/" No. moves  % accept.        E_tot           E_sr            E_Four       <E_self> "/1x,80("-"))')
            else
                Write(iothi,'(/" No. moves  % accept.        E_tot          E_sr   "/1x,80("-"))')
                Write(*,'(/" No. moves  % accept.        E_tot           E_sr    "/1x,80("-"))')
            endif
        elseif (ensemble == 'npt' ) Then
            Write(ioth,'(/" No. moves  % accept.   % accept. vol.   <E_tot>      <E_sr>    <Vol>      <Lx>     <Ly>    <Lz>"/1x,120("-"))')
            if (elect) then
                Write(iothi,'(/" No. moves  % accept.  % acc. vol.     E_tot           E_sr         E_coul     Vol"/1x,120("-"))')
                write(*,'(/" No. moves  % accept.   % acc. vol.     E_tot           E_sr            E_Four       <E_self>     Vol"/1x,120("-"))')
            else
                Write(iothi,'(/" No. moves  % accept.   % acc. vol.    E_tot          E_sr       Vol     Lx      Ly    Lz"/1x,120("-"))')
                Write(*,'(/" No. moves  % accept.    % acc. vol.    E_tot           E_sr      Vol"/1x,100("-"))')
            endif
        endif
    end subroutine init_printout

    Subroutine Printout(inst)
        use set_precision
        use properties, only : e_fourier, E_sr, Eref, Etotal, E_sav, E_coulomb, Etav, side_av, vol_av
        use potential, only : selfe, elect
        use configuration, only : natoms, v0, side
        use rundata, only : ntrial, naccept, naver, nvaccept, ensemble
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
                Write(ioth,'(i9,3x,f8.4,3x,2g15.7)')ntrial/natoms, 100*Dble(naccept)/Dble(ntrial),Etav/naver, E_sav/naver
            Elseif (ensemble == 'npt' ) then
                Write(ioth,'(i9,3x,2f9.4,3x,6g15.7)')ntrial/natoms, 100*Dble(naccept)/Dble(ntrial),prob_vol,Etav/naver, E_sav/naver, vol_av/naver, side_av(:)/naver
            Endif

        Else
            !
            ! Set energy reference to average value during equilibration.
            !
            npeq = npeq+1
            Etavq = Etavq+Etotal
            Eref = Etavq/npeq
        Endif
        If (ensemble == 'nvt' ) then

            if (elect) then
                Write(iothi,'(i9,3x,f8.4,3x,4g15.7)') ntrial/natoms, 100*Dble(naccept)/dble(ntrial),Etotal, E_sr, E_coulomb
                Write(*,'(i9,3x,f8.4,3x,6g15.7)') ntrial/natoms, 100*Dble(naccept)/Dble(ntrial),Etotal, E_sr, E_Fourier, selfe
            else
                Write(iothi,'(i9,3x,f8.4,3x,4g15.7)') ntrial/natoms, 100*Dble(naccept)/dble(ntrial),Etotal, E_sr
                Write(*,'(i9,3x,f8.4,3x,4g15.7)') ntrial/natoms, 100*Dble(naccept)/dble(ntrial),Etotal, E_sr
            endif
        ElseIf (ensemble == 'npt' ) then
            if (elect) then
                Write(iothi,'(i9,3x,2f9.4,3x,9g15.7)') ncycles, 100*Dble(naccept)/dble(ntrial), prob_vol, Etotal, E_sr, E_coulomb, v0, side(:)
                Write(*,'(i9,3x,2f9.4,3x,9g15.7)') ncycles, 100*Dble(naccept)/Dble(ntrial), prob_vol, Etotal, E_sr, E_Fourier, selfe, v0
            else
                Write(iothi,'(i9,3x,2f9.4,3x,9g15.7)') ncycles, 100*Dble(naccept)/dble(ntrial), prob_vol, Etotal, E_sr, v0, side(:)
                Write(*,'(i9,3x,2f9.4,3x,9g15.7)') ncycles, 100*Dble(naccept)/dble(ntrial), prob_vol, Etotal, E_sr, v0
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
        Open(igr,file="results/gmix.dat")
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
        Use properties, Only: E_sr,E_Fourier
        Write(*,'(80("-")/" vdW Energy=",g15.7,1x,a2)')E_sr,units
        If (elect) Then
            write(*, '(" K cutoff",f10.5," 1/A")')rcpcut
            write(*, '(" Coulomb energy (self energy) ",g15.7,1x,a2)') selfe, units
            write(*, '(" Coulomb energy (Fourier term)",g15.7,1x,a2)') E_Fourier, units
            write(*, '(" Coulomb energy ",g15.7,1x,a2)') selfe+eng_f,units
            write(*, '(" E_sr =",g15.7,1x,a2)') E_sr, units
            write(*, '(" Deviation from charge neutrality =",g15.7)') qtotal
        EndIf
        Write(*, '(" Etotal =",g15.7,1x,a2/80("-"))') selfe+E_Fourier+E_sr,units

    End Subroutine print_ener

end module output

