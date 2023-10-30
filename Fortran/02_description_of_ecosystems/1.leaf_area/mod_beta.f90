module MathFunctions
  implicit none

  private

  public :: beta

  contains

  ! Beta function (using the gamma function)
  real(4) function beta(p, q) result(beta_value)
    implicit none
    real(4), intent(in) :: p, q
    beta_value = GAMMA(p) * GAMMA(q) / GAMMA(p + q)
  end function beta

end module MathFunctions
