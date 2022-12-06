module setup
  !! This module contains routines to read in node distributions, modify them with shifting 
  !! pre-process to obtain boundary normal vectors, and write field data out to file.
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  use output
  
  implicit none
  
  real(rkind),dimension(:,:),allocatable :: rndshift_x,rndshift_y 
  integer(ikind) :: n_wvnmbrs
  real(rkind) :: wn_max,wn_min  

contains
!! ------------------------------------------------------------------------------------------------
  subroutine initial_setup  
     !! Initialises some key simulation parameters

     !! Time count begins at zero
     itime=0
  
     !! Initial data for profiling
     !$omp parallel
     n_threads = omp_get_num_threads()
     !$omp end parallel
     t_run = zero;t_last_X=zero

     !! Open a files for outputting
     open(unit=21,file='./data_out/time_out')
     open(unit=191,file='data_out/statistics/cputime.out')
     open(unit=192,file='data_out/statistics/dt.out')
     open(unit=195,file='data_out/statistics/ke.out')  
     open(unit=196,file='data_out/statistics/meantemp.out')
     open(unit=197,file='data_out/statistics/nusselt.out')
  
     !! Profiling:
     cputimecheck = zero
  
  end subroutine initial_setup
!! ------------------------------------------------------------------------------------------------
  subroutine setup_domain
     !! Reads in boundary patches, builds domain, calls decomposition and 
     !! boundary setup routines
     integer(ikind) i,j,ii,jj,ny_tmp,k,dummy_int,dummy_int2,nblock,nblock0,nm
     real(rkind) :: ns,dummy,prox,rad,radmin,xmin_local,xmax_local,x,y,xn,yn
     real(rkind),dimension(dims) :: rij
     real(rkind),dimension(:,:),allocatable :: tmp_vec
     integer(ikind) :: nl_ini,nl_end,nl_iniC,nl_endC,nl_ini2,nl_end2

     !! Load data
     open(unit=11,file='control.in')
     
     !! Read header
     read(11,*) 
     read(11,*)
     
     !! Read domain size
     read(11,*)
     read(11,*) Lx,Ly
     read(11,*)
     
     !! Read grid size
     read(11,*) 
     read(11,*) nx,ny
     read(11,*)
     
     !! Set the grid spacings
     dx = Lx/dble(nx-1)
     dy = Ly/dble(ny-1)  
     
     !! Read Prandtl and Rayleigh
     read(11,*) 
     read(11,*) Pr
     read(11,*)
     read(11,*) 
     read(11,*) Ra
     read(11,*)

     !! Start and end times
     read(11,*) 
     read(11,*) time,time_end
     read(11,*)     

     !! How often to output data?
     read(11,*) 
     read(11,*) dt_out 
     read(11,*)     
     
     !! Controls for power law
     read(11,*) 
     read(11,*) power_law_index
     read(11,*)
     
     !! Controls for Carreau fluid
     read(11,*)
     read(11,*) carreau_lambda,carreau_n,carreau_visc_inf
     read(11,*)

     !! Construct mesh
     oodx = Lx/dx;oody = Ly/dy        !! one/grid spacings
     write(6,*) "nx,ny",nx,ny
     allocate(rp(nx,ny,dims))
     do i=1,nx
        do j=1,ny
        
           !! Build pointers to transfer from linear to square index
           k=(i-1)*ny + j
           
           !! Build grid
           rp(i,j,1) = dx*dble(i-1)
           rp(i,j,2) = dy*dble(j-1)
        end do
     end do

     !! Allocate arrays for properties
     allocate(u(nx,ny),v(nx,ny),T(nx,ny),w(nx,ny),psi(nx,ny))
     allocate(alpha_out(nx,ny),visc(nx,ny))
     
     return
  end subroutine setup_domain
!! ------------------------------------------------------------------------------------------------
  subroutine initial_solution
     !! Temporary subroutine whilst building mfcomp. Initialises all fields
     integer(ikind) :: i,j,k,n_restart
     real(rkind) :: x,y,tmp,tmpro,tmp2
     character(70) :: fname
     

     !! Values within domain
     !$OMP PARALLEL DO PRIVATE(j,x,y,tmp,tmp2)
     do i=1,nx
        do j=1,ny
           x = rp(i,j,1);y=rp(i,j,2)

           !! Zero initial vorticity/velocity/streamfunction fields
           w(i,j) = zero           
           u(i,j) = zero 
           v(i,j) = zero  
           psi(i,j) = zero 
           
           !! Initial temperature field: hot bottom, cool top, and neutral in the centre. N.B. having
           !! a neutral region in the centre results in a more gentle acceleration of the flow in the 
           !! early stages.
           if(y/Ly.le.0.1d0 + 0.03d0*sin(2.0*pi*x/Lx)) then
              T(i,j) = one
           else if(y/Ly.ge.0.9d0 + 0.03d0*sin(2.0*pi*x/Lx)) then
              T(i,j) = -one
           else
              T(i,j) = zero
           end if



        end do
     end do
     !$OMP END PARALLEL DO
            
     
     return
  end subroutine initial_solution   
!! ------------------------------------------------------------------------------------------------
end module setup
