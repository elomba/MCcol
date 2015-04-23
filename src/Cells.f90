module cells
    use configuration, only : ndim, natoms, a, b, c, r
    use potential, only : rcut
contains
    Subroutine Init_cell
        Use linkcell
        use rundata, only : rdmax
        Implicit None
        !
        !   Only for orthogonal cells
        !
        Integer :: i,j,k,ii,jj,kk,l,in
        !
        ! Cell size extends up to rcut
        !
        maxi= Int(a(1)/(rcut+rdmax(1)))
        maxj= int(b(2)/(rcut+rdmax(2)))
        maxk= Int(c(3)/(rcut+rdmax(2)))
        ncell = maxi*maxj*maxk
        nn = 3**ndim
        Allocate(neigh(0:ncell-1,nn),head(0:ncell-1),list(natoms))
        cellx = 1.0d0/maxi
        celly = 1.0d0/maxj
        cellz = 1.0d0/maxk
        if (min(maxi,maxj,maxk) < 3) then
            print *, " *** Error: box too small for link cell method "
            use_cell = .false.
            return
        end if
        l=0
        Do i=0,maxi-1
            Do j=0,maxj-1
                Do k=0,maxk-1
                    in = 1
                    Do ii=-1,1
                        Do jj=-1,1
                            Do kk=-1,1
                                neigh(l,in) = fijk(i+ii,j+jj,k+kk)
                                in = in+1
                            End Do
                        End Do
                    End Do
                    l=l+1
                End Do
            End Do
        End Do
    end Subroutine Init_cell

    Integer Function fijk(ix,jx,kx)
        use linkcell, only : maxi, maxj, maxk
        Implicit None
        Integer :: ix, jx, kx, i, j, k
        i = ix
        j = jx
        k = kx
        !
        !  Use periodic boundary conditions
        !
        if (i < 0) i=i+maxi
        if (j < 0) j=j+maxj
        if (k < 0) k=k+maxk
        if (i >= maxi) i=i-maxi
        if (j >= maxj) j=j-maxj
        if (k >= maxk) k=k-maxk
        fijk=(i*maxj+j)*maxk+k
    end function fijk

    subroutine build_cells
        use linkcell
        implicit none
        integer :: n, i, j, k, icell
        head(:) = 0
        list(:) = 0
        do n = 1, natoms
            i=int((r(n,1)+0.5d0)/cellx)
            j=int((r(n,2)+0.5d0)/celly)
            k=int((r(n,3)+0.5d0)/cellz)
            icell=(i*maxj+j)*maxk+k
            list(n) = head(icell)
            head(icell) = n
        end do

    end subroutine build_cells

    subroutine update_cell_list(ntest,icell)
        !
        ! Check whether particle ntest exits cell ocell, update cell list
        ! if needed
        !
        use linkcell, only : maxk, maxj, cellx, celly, cellz, head, list
        implicit none
        integer, intent(IN) :: ntest, icell
        integer :: i,j,k, io, ocell
        ! r contains old particle coordinates
        i=int((r(ntest,1)+0.5d0)/cellx)
        j=int((r(ntest,2)+0.5d0)/celly)
        k=int((r(ntest,3)+0.5d0)/cellz)
        ocell=(i*maxj+j)*maxk+k
        ! Check if particle has left its cell
        if (icell .ne. ocell) then
            !
            ! Remove ntest from the list of ocell
            !
            io = ntest
            i = head(ocell)
            if (i == ntest) then
                head(ocell) = list(ntest)
            else
                do while ( i .ne. ntest )
                    io = i
                    i = list(i)
                enddo
                list(io) = list(ntest)
            endif
            !
            ! Add ntest to the top of the list of ncell.
            !
            io = head(icell)
            head(icell) = ntest
            list(ntest) = io
        !           print *, ' Updating list...'
        endif

    end subroutine update_cell_list
end module cells
