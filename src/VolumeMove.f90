Module VolumeChange

    Use Set_precision
    Use Interfaces
    Use Linkcell, Only : use_cell
    Use Thermo
    Use Energy, Only : Energ, Energ_cell
    Use Init, Only : Init_fourier

Contains 

    Subroutine move_volume
        implicit none

        if (scaling == 'cubic') then

            print *, " ** Input error: cubic scaling not implemented yet ... "
            stop

        elseif (scaling == 'ortho'.or.scaling == 'isotr') then

            call move_volume_ortho

        Endif

    End subroutine move_volume




    Subroutine move_volume_ortho
        !
        !   If the ratio between the edges of the box is maintained
        !     we could do some scaling, but we need to define kappa_x, kappa_y and kappa_z
        !     that possibility is not yet implemented
        !
        use cells
        use linkcell
        implicit none
        Real (wp)  ::  xi, deltaEt, deltaE_sr
        Real (wp)  ::  Etotal_old, E_sr_old, E_coulomb_old, E_fourier_old, v_old
        Real (wp)  , dimension(ndim) :: a_old, b_old, c_old, side_old, r_unit_old
        Real (wp)  ::  deltaE_Fourier, rcpcut2_old, factor, rcpcut_old
        Real (wp)  ::  dospix_old, dospiy_old, dospiz_old
        Complex(wp) :: rhokk_old(kmt), ekm2_old(kmt), km2_old(0:kmt)
        Complex(wp) :: eix_old(1:natoms,-kmx:kmx),eiy_old(1:natoms,-kmy:kmy),eiz_old(1:natoms,0:kmz)
        if (use_cell) then
            call store_cell
        endif
        E_sr_old = E_sr
        side_old=side
        a_old=a
        b_old=b
        c_old=c
        r_unit_old=r_unit
        v_old=v0

        If (elect) Then
            Etotal_old = Etotal  !Etot and Ecoul are not updated in move
            E_coulomb_old = E_coulomb
            E_fourier_old = E_fourier
            ekm2_old = ekm2
            km2_old = km2
            rcpcut2_old = rcpcut2  ! just in case there are large volume changes
            rcpcut_old = rcpcut  ! just in case there are large volume changes
            rhokk_old = rhokk   !
            eix_old = eix   ! for rigid molecules, we need also to actualize eix,eiy,eiz
            eiy_old = eiy   ! or for flexible molecules that are not scale with volume
            eiz_old = eiz   ! for atomic systems it is actually not necessary
            dospix_old = dospix
            dospiy_old = dospiy
            dospiz_old = dospiz
        EndIf

        Call Random_number(xi)
        side(1) = side(1) + (2*xi-1)*vdmax

        If (scaling == 'ortho') Then  ! the edges change independently

            Call Random_number(xi)
            side(2) = side(2) + (2*xi-1)*vdmax
            Call Random_number(xi)
            side(3) = side(3) + (2*xi-1)*vdmax

        Else   ! the edges change by the same factor

            factor = side(1) / side_old (1)
            side(2) = side(2) *factor
            side(3) = side(3) *factor

        Endif

        v0=side(1)*side(2)*side(3)

        a(1) = side(1)
        b(2) = side(2)
        c(3) = side(3)
        r_unit(:) = side(:)

        If (elect) call Init_fourier
        if (use_cell) then
            call init_cell
            Call Energ_cell
        else
            Call Energ
        endif

        deltaE_sr = E_sr - E_sr_old
        If (elect) Then
            deltaE_Fourier = E_fourier - E_fourier_old
        Else
            deltaE_Fourier = 0.d0
        Endif
        deltaEt =  deltaE_sr + deltaE_Fourier

        Call Random_number(xi)
        If (Exp(-(deltaEt+pres*(v0-v_old))/kT+dble(natoms)*Log(v0/v_old)) .Gt. xi) Then
            nvaccept = nvaccept + 1
        Else   ! reject
            if (use_cell) call restore_cell
            Etotal = Etotal_old  ! etot and ecoulomb are not updated in Move
            E_sr = E_sr_old
            a = a_old
            b = b_old
            c = c_old
            side = side_old
            r_unit = r_unit_old
            v0 = v_old
            If (elect) Then
                E_coulomb = E_coulomb_old
                E_fourier = E_fourier_old
                ekm2 = ekm2_old
                km2 = km2_old
                rcpcut2 = rcpcut2_old
                rcpcut = rcpcut_old
                rhokk = rhokk_old
                eix = eix_old   ! for rigid molecules, we also need to update eix
                eiy = eiy_old
                eiz = eiz_old
                dospix = dospix_old
                dospiy = dospiy_old
                dospiz = dospiz_old
            Endif
        Endif

    End subroutine move_volume_ortho
 
End Module VolumeChange
