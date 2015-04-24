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
        nei = Nint((Etotal-Eref)/deltaEng)
        If (Abs(nei) > nde) Then
            Print *, " Warning : E histogram out of bounds. Increase  nde up to ", Abs(nei)
        Else
            Ehisto(nei) = Ehisto(nei)+1
        End If
    End Subroutine histo_energ

    Subroutine histo_rho
        implicit none
        integer :: nrho
        real(wp) :: rhoins
        rhoins = natoms/v0
        nrho = Nint(rhoins/deltar)
        If (Abs(nrho) > nrmax) Then
            Print *, " Warning : rho histogram out of bounds. Increase  nrmax up  to ", nrho
        Else
            rhisto(nrho) = rhisto(nrho)+1

        End If
    end subroutine histo_rho

    Subroutine averages
        Use rundata, only : ensemble
        use output, only : Printout
        Implicit None
        Real(wp) :: Esum, rsum
        Integer :: nei, nhstmin, nhstmax, i
        Logical, Save :: first=.True.
        If (first) Then
            Write(*,"(' *** End of equilibration ***')")
            If (ensemble == 'nvt') Then

                If (elect) Then
                    Write(*,1000)
1000                format(/" No. moves  % accept.        E_tot           E_sr            E_Four       <E_self> "&
                        &/1x,80("-"))
                Else
                    Write(*,1010)
1010                format(/" No. moves  % accept.        E_tot           E_sr    "/1x,80("-"))
                Endif
            Elseif (ensemble == 'npt') Then
                If (elect) Then
                    Write(*,1020)
1020                format(/" No. moves  % accept.   % accept.vol     E_tot           E_sr ",&
                        &"           E_Four       <E_self>   Vol"/1x,100("-"))
                Else
                    Write(*,1030)
1030                format(/" No. moves  % accept.   % accept.vol     E_tot           E_sr    Vol"/1x,120("-"))
                Endif
            endif
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
            Open (ioh,file='results/Ehisto.dat')
            Do i=Max(nhstmin,-nde),Min(nhstmax,nde)
                Write(ioh,'(2f15.7)')(i+0.5)*deltaEng+Eref,Ehisto(i)/Esum
            End Do
            Close(ioh)
        ElseIf (ensemble == "npt") Then
            Vol_av = Vol_av+v0
            side_av(:) = side_av(:) + side(:)
            rsum = Sum(rhisto(0:nrmax))
            Open (ioh,file='results/rho_histo.dat')
            do i=0,nrmax
                if (rhisto(i) > 0) then
                    write(ioh,'(2f15.7)')(i+0.5)*deltar, rhisto(i)/rsum
                endif
            enddo
        End If
        close(ioh)
        naver = naver+1
        Etotal = E_sr+E_Fourier+selfe
        Etav = Etotal+Etav
        E_sav = E_sr+E_sav
        call Printout(.true.)
    End Subroutine averages


End Module Thermo

subroutine histograms(ensemble)
    use Thermo, only : histo_energ, histo_rho
    implicit none
    character, intent(IN) :: ensemble*3
    if (ensemble == "nvt") then
        call histo_energ
    else if (ensemble == "npt") then
        call histo_rho
    end if
end subroutine histograms
