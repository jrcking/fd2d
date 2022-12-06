module output
  !! This module contains routines to read in node distributions, modify them with shifting 
  !! pre-process to obtain boundary normal vectors, and write field data out to file.
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  implicit none

contains
!! ------------------------------------------------------------------------------------------------
  subroutine output_to_screen
     !! This routine writes information about the simulation to screen in a fairly easy to read
     !! format. For this to work (nicely) terminal window should be 24 lines tall.
     integer(ikind) :: scr_freq=10
     real(rkind),dimension(6) :: maxphi,minphi
     integer(ikind) :: n_threads_global
     real(rkind) :: t_per_dt_global,t_last_x_global,t_run_global
     real(rkind) :: cput,cput_global
    
     ts_end=omp_get_wtime()
     t_run = t_run + ts_end - ts_start
     t_per_dt = t_run/dble(itime)
     t_last_X = t_last_X + ts_end - ts_start  
     !! Output cpu-time to file.

     write(191,*) itime,ts_end-ts_start
     flush(191)  
  
     ! Some to screen
     if(mod(itime,scr_freq).eq.0)then 
  
        write(6,*)"itime,time,dt=", itime,time,dt
        write(6,*) "nx,ny",nx,ny
        write(6,*) "max velocity = ",umax
        write(6,*) "Number of threads=",n_threads,"Run time=",t_run
        write(6,*) "run-time/dt=",t_per_dt,"Moving avg=",t_last_X/dble(scr_freq)
        write(6,*) "SOR iterations and residual",itercount,residual
        t_last_X = 0.0d0
        
        write(6,'(/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,A)') "  "                           
     
     end if
   
     ts_start=omp_get_wtime()
     return
  end subroutine output_to_screen    
!! ------------------------------------------------------------------------------------------------  
  subroutine output_layer(n_out)
     !! Little subroutine to write out field variables. 
     !! Can be converted to vtk and read into paraview.
     use derivatives
     integer(ikind),intent(in) :: n_out
     integer(ikind) :: i,j,k,np_out_local
     character(70) :: fname
     real(rkind) :: tmpT,xn,yn


     !! set the name of the file...
     !! first number is processor number, second is dump number (allowed up to 9999 processors)


     if( n_out .lt. 10 ) then 
        write(fname,'(A17,I1)') './data_out/layer_',n_out        
     else if( n_out .lt. 100 ) then 
        write(fname,'(A17,I2)') './data_out/layer_',n_out        
     else if( n_out .lt. 1000 ) then
        write(fname,'(A17,I3)') './data_out/layer_',n_out        
     else
        write(fname,'(A17,I4)') './data_out/layer_',n_out        
     end if 
     
     !! Write the main dump files
     open(unit = 20,file=fname)  
     write(20,*) nx*ny
     do i=1,nx
        do j=1,ny
                  
           !! x,y,u,v,vorticity,temperature,alpha_out                     
           write(20,*) rp(i,j,1),rp(i,j,2),u(i,j),v(i,j),w(i,j),T(i,j),alpha_out(i,j)
        end do
     end do

     flush(20)
     close(20)

     !! Write the time,dump number and # nodes to file
     write(21,*) time,nx,ny,n_out,1
     flush(21)
     
     return
  end subroutine output_layer
!! ------------------------------------------------------------------------------------------------
  subroutine statistics
     !! Control subroutine from which to call routines to evaluate the output statistics          
          
     !! Check total kinetic energy
     call check_ke
     
     !! Check the mean temperature
     call check_temp
     
     !! Evaluate the Nusselt number
     call calc_nusselt

     return
  end subroutine statistics 
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_nusselt
     !! Placeholder routine for Nusselt calculation
     integer(ikind) :: i,j
     real(rkind) :: nusselt
       
    

!     write(197,*) time,nusselt           !!! writing to unit 197 will write to the file 'data_out/statistics/nusselt.out'
!     flush(197)
     
     
     return
  end subroutine calc_nusselt    
!! ------------------------------------------------------------------------------------------------  
  subroutine check_ke
     !! Output the L2 of velocity over the domain
     integer(ikind) :: i,j
     real(rkind) :: tot_vel,tmpvel
       
     tot_vel = zero
     !$omp parallel do private(j,tmpvel) reduction(+:tot_vel)
     do i=1,nx
        do j=1,ny             
           tmpvel = (u(i,j)*u(i,j)+v(i,j)*v(i,j))
           tot_vel = tot_vel + tmpvel
        end do
     end do
     !$omp end parallel do
          
     !! Normalise over volume
     tot_vel = sqrt(tot_vel/dble(nx*ny))
     

     write(195,*) time,tot_vel
     flush(195)
     
     
     return
  end subroutine check_ke   
!! ------------------------------------------------------------------------------------------------
  subroutine check_temp
     !! This subroutine calculates the mean temperature
     integer(ikind) :: i,j
     real(rkind) :: tot_T
   
     tot_T = zero
     !$omp parallel do private(j) reduction(+:tot_T)
     do i=1,nx
        do j=1,ny      
           tot_T = tot_T + T(i,j)
        end do        
     end do
     !$omp end parallel do
     
     !! Mean temperature        
     tot_T = tot_T/dble(nx*ny)

     write(196,*) time,tot_T
     flush(196)

     return
  end subroutine check_temp
!! ------------------------------------------------------------------------------------------------ 
end module output
