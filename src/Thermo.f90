Module Thermo
    use interfaces
    Use set_precision
    Use configuration
    Use potential
    Use properties
    Use rundata
Contains 

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
1000                format(/" No. moves  % accept.        E_tot           E_sr            E_vdw          E_Four       <E_self>           <E_coul>"&
                        &/1x,120("-"))
                Else
                    Write(*,1010)
1010                format(/" No. moves  % accept.        E_tot           E_sr         E_vdw       "/1x,80("-"))
                Endif
            Elseif (ensemble == 'npt') Then
                If (elect) Then
                    Write(*,1020)
1020                format(/" No. moves  % accept.   % accept.vol     E_tot           E_sr       E_coul",&
                        &"           E_vdw         E_Four       <E_self>   Vol"/1x,100("-"))
                Else
                    Write(*,1030)
1030                format(/" No. moves  % accept.   % accept.vol     E_tot           E_sr        E_vdw       Vol"/1x,120("-"))
                Endif
            endif
            first = .False.
        End If
        naver = naver+1
        Etotal = E_sr+E_Fourier+selfe
        Etav = Etotal+Etav
        E_sav = E_sr+E_sav
        E_vdwav = E_vdwav + Evdw
        call Printout(.true.)
    End Subroutine averages


End Module Thermo


