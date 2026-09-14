!===============================================================================
! Function: dist2
!
! Purpose:
!   Calculates the squared Euclidean distance between two particles applying
!   the minimum image convention (MIC) under periodic boundary conditions (PBC)
!   for an orthorhombic simulation cell.
!
! Arguments:
!   r (Input)  : Real displacement vector [ndim] in reduced coordinates,
!                i.e. normalized by the box dimensions (r_scaled in [-1, 1]).
!   dist2      : Squared physical distance |r_ij|^2 in Angstroms^2.
!
! Method:
!   1. Applies minimum image convention by wrapping components outside [-0.5, 0.5):
!        if r_alpha >  0.5 -> r_alpha = r_alpha - 1
!        if r_alpha < -0.5 -> r_alpha = r_alpha + 1
!   2. Converts reduced displacement back to physical dimensions:
!        r_alpha = r_alpha * r_unit_alpha  (Angstroms)
!   3. Computes dot product: dist2 = sum(r_alpha^2).
!===============================================================================
Function dist2(r)
  use set_precision
  Use configuration, Only : r_unit, ndim
  Implicit None
  Real(wp) :: r(ndim)
  Real(wp) :: dist2

  ! Apply minimum image convention in reduced coordinates [-0.5, 0.5)
  Where (r(:) >  0.5d0) r(:) = r(:) - 1.0d0
  Where (r(:) < -0.5d0) r(:) = r(:) + 1.0d0

  ! Convert reduced coordinates to physical lengths in Angstroms
  r(:) = r_unit(:) * r(:)

  ! Return squared distance (Angstrom^2)
  dist2 = Dot_product(r(:), r(:))
End Function dist2
