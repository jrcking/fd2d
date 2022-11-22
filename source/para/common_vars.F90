module common_vars

  use iso_c_binding
  use common_parameter
  use kind_parameters      

  implicit none

  !! Control variables
  real(rkind) :: Lx,Ly  !! Domain size
  integer(ikind) :: nx,ny !! grid size
  real(rkind) :: Pr !! Prandtl number
  real(rkind) :: Ra !! Rayleigh number 
  real(rkind) :: time,time_end,dt_out !! Start time, end time, output period
  
  !! Non-Newtonian parameters
  real(rkind) :: power_law_index
  real(rkind) :: carreau_lambda,carreau_visc_inf,carreau_n

  !! Evolved fluid properties
  real(rkind), dimension(:,:), allocatable, target :: T,w
  
  !! Evaluated fluid properties
  real(rkind), dimension(:,:), allocatable, target :: u,v,psi
  
  !! Additional array for output
  real(rkind), dimension(:,:), allocatable :: alpha_out  
  
  !! Right hand sides
  real(rkind),dimension(:,:), allocatable :: rhs_u,rhs_v,rhs_T,rhs_w
  
  !! Secondary fluid properties
  real(rkind), dimension(:,:), allocatable, target :: visc
  
  !! Performance of SOR method:
  integer(ikind) :: itercount
  real(rkind) :: residual  
  
  !! Discretisation properties
  real(rkind), dimension(:,:,:), allocatable, target :: rp
  
  !! grid spacing, 1/grid spacing
  real(rkind) :: dx,dy,oodx,oody
  

  !! Parameters related to time and some forces etc
  real(rkind) :: dt
  real(rkind) :: umax,smax                    !! maximum velocity and node spacing   
  integer(ikind) :: itime

  !! Finite Difference weightings 
  integer(ikind) :: ij_count_fd ! Size of FD stencil  
  integer(ikind),dimension(:,:),allocatable :: fd_link_x,fd_link_y
  real(rkind),dimension(:,:),allocatable :: fd_grad_x,fd_grad2_x,fd_hyp_x
  real(rkind),dimension(:,:),allocatable :: fd_grad_y,fd_grad2_y,fd_hyp_y
   
  !! Profiling and openMP parallelisation
  real(rkind) ts_start,ts_end,t_run,t_per_dt,t_last_X
  integer(ikind) :: n_threads  
  real(rkind) :: cputimecheck
 
          
end module common_vars
