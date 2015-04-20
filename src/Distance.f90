Function dist2(r)
  !
  ! Euclidean distance in an orthogonal cell. Input in box side
  ! units, output in program input units (A).
  !
  use set_precision
  Use configuration, Only : r_unit, ndim
  Implicit None
  Real(wp) :: r(ndim)
  Real(wp) :: dist2
  Where (r(:) > 0.5d0 ) r(:) = r(:) - 1
  Where (r(:) < -0.5d0 ) r(:) = r(:) + 1
  r(:) = r_unit(:)*r(:)
  dist2 = Dot_product(r(:),r(:))
End Function dist2
