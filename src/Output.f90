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
         &,a2/)')name,date(7:8),date(5:7),date(1:4),time(1:2)&
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
         &,a2/"  Elapsed CPU time",f15.2," s")')date(7:8),date(5:7),date(1:4),time(1:2)&
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
    implicit none
    Open(ioth,file='thermoaver.dat')
    Open(iothi,file='thermoins.dat')
    Write(ioth,'(/" No. moves  % accept.       <E_tot>         <E_sr>"/1x,80("-"))')
    if (elect) then
       Write(iothi,'(/" No. moves  % accept.        E_tot           E_sr         E_coul "/1x,80("-"))')
       write(*,'(/" No. moves  % accept.        E_tot           E_sr            E_Four       <E_self> "/1x,80("-"))')
    else
       Write(iothi,'(/" No. moves  % accept.        E_tot          E_sr   "/1x,80("-"))')
       Write(*,'(/" No. moves  % accept.        E_tot           E_sr    "/1x,80("-"))')
    endif
  end subroutine init_printout

  Subroutine Printout(inst)
    use set_precision
    use properties, only : e_fourier, E_sr, Eref, Etotal, E_sav, E_coulomb, Etav
    use potential, only : selfe, elect
    use configuration, only : natoms
    use rundata, only : ntrial, naccept, naver
    Implicit None
    logical, intent(IN) :: inst
    Real(wp), Save ::  Etavq=0 
    Integer, save :: npeq = 0
    Etotal = E_sr+E_Fourier+selfe
    If (inst) Then
       Write(ioth,'(i9,3x,f8.4,3x,g15.7)')ntrial/natoms, 100*Dble(naccept)/Dble(ntrial),Etav/naver, E_sav/naver
    Else
       !
       ! Set energy reference to average value during equilibration.
       !
       npeq = npeq+1
       Etavq = Etavq+Etotal
       Eref = Etavq/npeq
    Endif
    if (elect) then
       Write(iothi,'(i9,3x,f8.4,3x,4g15.7)') ntrial/natoms, 100*Dble(naccept)/dble(ntrial),Etotal, E_sr, E_coulomb
       Write(*,'(i9,3x,f8.4,3x,6g15.7)') ntrial/natoms, 100*Dble(naccept)/Dble(ntrial),Etotal, E_sr, E_Fourier, selfe
    else
       Write(89,'(i9,3x,f8.4,3x,4g15.7)') ntrial/natoms, 100*Dble(naccept)/dble(ntrial),Etotal, E_sr
       Write(*,'(i9,3x,f8.4,3x,4g15.7)') ntrial/natoms, 100*Dble(naccept)/dble(ntrial),Etotal, E_sr
    endif
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
end module output

