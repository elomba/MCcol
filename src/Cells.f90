!===============================================================================
! Module: cells
!
! Purpose:
!   Implements the 3D Link-Cell (Cell-List) spatial decomposition method.
!   Reduces short-range pair interaction computations from O(N^2) to O(N).
!
! Algorithm Overview:
!   1. The simulation box is partitioned into a 3D grid of subcells of edge
!      length >= (rcut + rdmax). This guarantees that any particle within rcut
!      of particle i in cell C must reside either in cell C or in one of its
!      26 immediate adjacent neighboring cells (total 27 cells).
!   2. Atoms in each cell are stored as a singly linked list via the arrays:
!      - head(icell) : Points to the first atom in cell icell (or 0 if empty).
!      - list(iatom) : Points to the next atom in the same cell (or 0 at end).
!   3. Periodic boundary conditions are handled by wrapping cell indices via fijk.
!   4. During MC moves, only atoms that cross a cell boundary trigger an O(1)
!      re-linking of the head and list arrays via update_cell_list.
!===============================================================================
module cells
    use configuration, only : ndim, natoms, a, b, c, r
    use potential, only : rcut
contains

    !---------------------------------------------------------------------------
    ! Subroutine: Init_cell
    !
    ! Purpose:
    !   Determines cell grid dimensions (maxi, maxj, maxk), allocates cell arrays,
    !   and precalculates the 27 neighbor cell indices for each cell.
    !---------------------------------------------------------------------------
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
        maxk= Int(c(3)/(rcut+rdmax(3)))
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

    !---------------------------------------------------------------------------
    ! Function: fijk
    !
    ! Purpose:
    !   Maps 3D grid cell indices (ix, jx, kx) to a flat 1D cell index in [0, ncell-1]
    !   with periodic wrapping along all three spatial dimensions.
    !---------------------------------------------------------------------------
    Integer Function fijk(ix,jx,kx)
        use linkcell, only : maxi, maxj, maxk
        Implicit None
        Integer :: ix, jx, kx, i, j, k
        i = ix
        j = jx
        k = kx
        ! Apply periodic boundary wrapping
        if (i < 0) i=i+maxi
        if (j < 0) j=j+maxj
        if (k < 0) k=k+maxk
        if (i >= maxi) i=i-maxi
        if (j >= maxj) j=j-maxj
        if (k >= maxk) k=k-maxk
        fijk=(i*maxj+j)*maxk+k
    end function fijk

    !---------------------------------------------------------------------------
    ! Subroutine: build_cells
    !
    ! Purpose:
    !   Rebuilds the entire link-cell structure from scratch for all atoms.
    !   Resets head and list arrays, assigns each atom to its respective cell.
    !---------------------------------------------------------------------------
    subroutine build_cells
        use linkcell
        implicit none
        integer :: n, i, j, k, icell
        head(:) = 0
        list(:) = 0
        do n = 1, natoms
            i=int((R(1,n)+0.5d0)/cellx)
            j=int((R(2,n)+0.5d0)/celly)
            k=int((R(3,n)+0.5d0)/cellz)
            icell=(i*maxj+j)*maxk+k
            list(n) = head(icell)
            head(icell) = n
        end do
    end subroutine build_cells

    !---------------------------------------------------------------------------
    ! Subroutine: update_cell_list
    !
    ! Purpose:
    !   Updates cell linked lists when a single particle (ntest) moves to a new cell.
    !   Removes ntest from its old cell list and prepends it to the new cell list in O(1) time.
    !
    ! Arguments:
    !   ntest (in) : Atom index that was displaced.
    !   icell (in) : Destination cell index.
    !---------------------------------------------------------------------------
    subroutine update_cell_list(ntest,icell)
        use linkcell, only : maxk, maxj, cellx, celly, cellz, head, list
        implicit none
        integer, intent(IN) :: ntest, icell
        integer :: i,j,k, io, ocell
        ! r contains old particle coordinates
        i=int((R(1,ntest)+0.5d0)/cellx)
        j=int((R(2,ntest)+0.5d0)/celly)
        k=int((R(3,ntest)+0.5d0)/cellz)
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
