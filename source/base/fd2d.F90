program fd2d
  use kind_parameters
  use common_parameter
  use common_vars
  use setup
  use output 
  use weights
  use step
  implicit none

  integer(ikind) :: n_out,m_out



  !! Initial conditions
  call initial_setup  
  call setup_domain

  !! Calculate all finite difference weights and link lists we need
  call build_operators
   
  !! Create initial fields for primary variables
  call initial_solution

  !! Initialise the time-step (to something small to be safe...)
  call set_tstep;dt=0.0001*dt  

  !! Initialise time profiling and output counter...
  n_out = 0;ts_start = omp_get_wtime()
  m_out = 0
        
  !! MAIN TIME LOOP ---------------------------------------------------
  do while (time.le.time_end)
    
     !! Output, conditionally: at start, subsequently every dt_out
     if(itime.eq.0.or.time.gt.n_out*dt_out) then 
!     if(itime.eq.0.or.mod(itime,10).eq.0)then
        n_out = n_out + 1
        call output_layer(n_out)
     end if        
    
     !! Set the time step
     call set_tstep     

     !! Perform one time step
     call step_euler

     !! Calculate time-profiling and output to screen
     itime = itime + 1
     call output_to_screen

     !! Call routines to evaluate global statistics and adjust forcing terms if desired
     call statistics
     
  end do
  !! END MAIN TIME LOOP -----------------------------------------------
  
  !! Deallocate particle properties and neighbour lists
  call deallocate_weights
  stop
end program fd2d
!! ------------------------------------------------------------------------------------------------
subroutine deallocate_weights
  use kind_parameters
  use common_parameter
  use common_vars

  deallocate(rp,u,v,w,T)
  deallocate(fd_link_x,fd_link_y)
  deallocate(fd_grad_x,fd_grad2_x)
  deallocate(fd_grad_y,fd_grad2_y)  

  return
end subroutine deallocate_weights
