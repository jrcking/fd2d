module common_parameter
  use kind_parameters
  implicit none 

  !! Discretisation related parameters
  integer(ikind), parameter :: npar=999999  !! Only used in source/gen/datclass.F90 (up to npar in a slice)
  integer(ikind) ,parameter :: dims = 2

  !! Numbers --------------------------------------------------------------------------------------
  real(rkind), parameter :: pi=3.141592653589793238462643383279502884197d0
  real(rkind), parameter :: oosqrt2 = 1.0d0/dsqrt(2.0d0)
  real(rkind), parameter :: zero = 0.0d0
  real(rkind), parameter :: one = 1.0d0
  real(rkind), parameter :: half = 0.5d0
  real(rkind), parameter :: verysmall = 1.0d-30
  real(rkind), parameter :: two = 2.0d0  
  

   
end module common_parameter
