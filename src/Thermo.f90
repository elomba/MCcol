Module Thermo
  use interfaces
  Use set_precision
  Use configuration
  Use potential
  Use properties
  Use rundata
Contains 

  Subroutine histo_energ
    Implicit None
    Integer :: nei
    Etotal = E_sr+E_Fourier+selfe
    nei = Nint((Etotal-Eref)/deltaE)
    If (Abs(nei) > nde) Then
       Print *, " Warning : histogram out of bounds. Increase (decrease) nde up (down) to ", nei 
    Else
       Ehisto(nei) = Ehisto(nei)+1
    End If
   End Subroutine histo_energ

  Subroutine averages
    Use rundata, only : ensemble
    use output, only : Printout
    Implicit None
    Real(wp) :: Esum
    Integer :: nei, nhstmin, nhstmax, i
    Logical, Save :: first=.True.
    If (first) Then
       Write(*,"(' *** End of equilibration ***')")
       If (elect) Then
          Write(*,'(/" No. moves  % accept.        E_tot           E_sr            E_Four       <E_self> "/1x,80("-"))')
       Else
          Write(*,'(/" No. moves  % accept.        E_tot           E_sr    "/1x,80("-"))')
       Endif
       first = .False.
    End If
    !
    ! Printout energy histogram for NVT (so far ..)
    !
    If (ensemble == "nvt") Then
       i=-nde
       Do While (Abs(Ehisto(i)) < 1.d-8 .And. i <= nde-1 )
          i=i+1
       End Do
       nhstmin = i
       i=nde
       Do While (Abs(Ehisto(i)) < 1.d-8 .And. i >= -nde+1)
          i=i-1
       End Do
       nhstmax = i
       Esum = Sum(Ehisto(-nde:nde))
       Open (14,file='results/Ehisto.dat')
       Do i=Max(nhstmin,-nde),Min(nhstmax,nde)
          Write(14,'(2f15.7)')(i+0.5)*deltaE+Eref,Ehisto(i)/Esum
       End Do
       Close(14)
    End If
    naver = naver+1
    Etotal = E_sr+E_Fourier+selfe
    Etav = Etotal+Etav
    E_sav = E_sr+E_sav
    call Printout(.true.)
  End Subroutine averages


End Module Thermo

subroutine histograms(ensemble)
  use Thermo, only : histo_energ
  implicit none
  character, intent(IN) :: ensemble*3
  if (ensemble == "nvt") then
     call histo_energ
  end if
end subroutine histograms
