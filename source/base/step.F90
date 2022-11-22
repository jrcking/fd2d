module step
  !! This module contains time-stepping routines, and a routine to set the time-step.
  !! Time-stepping routines start with "step_" and perform one time step. They are
  !! called from the main loop, and they themselves call routines in the rhs module.

  use kind_parameters
  use common_parameter
  use common_vars
  use rhs
  implicit none

  
contains
!! ------------------------------------------------------------------------------------------------
  subroutine step_euler
     use derivatives
     integer(ikind) :: i,k,j
     real(rkind) :: time0
     real(rkind),dimension(:,:),allocatable :: w_reg1,T_reg1
     real(rkind),dimension(3) :: RKa
     real(rkind),dimension(4) :: RKb
         
     !! Allocate right hand sides
     allocate(rhs_w(nx,ny),rhs_T(nx,ny))  


     !! Temporary storage of time
     time0=time

     !! Build the right-hand sides
     call calc_rhs
     
     !! Set new w,T and _reg1 stores
     !$omp parallel do private(j)
     do i=1,nx
        do j=1,ny
           !! Store next U in register 2
           w(i,j) = w(i,j) + dt*rhs_w(i,j)
           T(i,j) = T(i,j) + dt*rhs_T(i,j)
     

        end do
     end do
     !$omp end parallel do
             
     !! Deallocation    
     deallocate(rhs_w,rhs_T)
  
     
     !! De-alias the vorticity equation via filtering!!
!     call calc_filtered_var(w)

     !! Set the new time   
     time = time0 + dt
     

     return
  end subroutine step_euler
!! ------------------------------------------------------------------------------------------------
  subroutine set_tstep
     integer(ikind) :: i,j
     real(rkind) :: umag,smin,dt_local
    
     !! Find maximum velocity magnitude
     umax = 1.0d-16
!     !$OMP PARALLEL DO PRIVATE(j,umag) REDUCTION(max:umax)
     do i=1,nx
        do j=1,ny
           umag = u(i,j)*u(i,j) + v(i,j)*v(i,j)
           if(umag.ge.umax) umax = umag     
        end do
     end do
!     !$OMP END PARALLEL DO
     umax = sqrt(umax)
     
     !! Find smallest node spacing
     smin = min(dx,dy)

     !! Set time step
     dt = min(min( &
                  0.05*smin*smin*sqrt(Ra/Pr), &      !! Viscous constraint
                  half*half*smin/(umax) ), &          !! Advection constraint
                  0.05*smin*smin*sqrt(Ra*Pr) )        !! Thermal conductivity constraint
   
     write(192,*) time,dt
     flush(192)         

     return
  end subroutine set_tstep
!! ------------------------------------------------------------------------------------------------  
end module step
