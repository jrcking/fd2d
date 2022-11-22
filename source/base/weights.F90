module weights
!! This module constructs weights for centred finite differences in the 3rd dimension
!! 8th order, with built in periodicity

  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib  
  implicit none 
  
  private
  public build_operators  
  
  real(rkind),dimension(:),allocatable :: fd_grad_1,fd_grad2_1,fd_hyp_1
  real(rkind),dimension(:),allocatable :: fd_grad_2,fd_grad2_2,fd_hyp_2
  real(rkind),dimension(:),allocatable :: fd_grad_3,fd_grad2_3,fd_hyp_3
  real(rkind),dimension(:),allocatable :: fd_grad_4,fd_grad2_4,fd_hyp_4
  real(rkind),dimension(:),allocatable :: fd_grad,fd_grad2,fd_hyp    
  real(rkind),dimension(:),allocatable :: fd_grad_n,fd_grad2_n,fd_hyp_n
  real(rkind),dimension(:),allocatable :: fd_grad_nm1,fd_grad2_nm1,fd_hyp_nm1
  real(rkind),dimension(:),allocatable :: fd_grad_nm2,fd_grad2_nm2,fd_hyp_nm2
  real(rkind),dimension(:),allocatable :: fd_grad_nm3,fd_grad2_nm3,fd_hyp_nm3      
  
contains
  subroutine build_operators
     integer(ikind) :: i,j
     integer(ikind) :: k,kk,ik,jk,num     
            
     allocate(fd_link_x(nx,9),fd_link_y(ny,9))
     allocate(fd_grad_x(nx,9),fd_grad2_x(nx,9),fd_hyp_x(nx,9))
     allocate(fd_grad_y(ny,9),fd_grad2_y(ny,9),fd_hyp_y(ny,9))     

     !! Build FD neighbours...
     !! These links are built for periodic domains. If the domain is non-periodic, that is
     !! included in the stencil weights, but the links still point to periodic neighbours
     do i=1,nx
        fd_link_x(i,5) = i
        if(i.eq.1) then
          fd_link_x(i,1) = nx-3
          fd_link_x(i,2) = nx-2
          fd_link_x(i,3) = nx-1
          fd_link_x(i,4) = nx          
        else if(i.eq.2) then
          fd_link_x(i,1) = nx-2
          fd_link_x(i,2) = nx-1
          fd_link_x(i,3) = nx
          fd_link_x(i,4) = i-1        
        else if(i.eq.3) then
          fd_link_x(i,1) = nx-1
          fd_link_x(i,2) = nx
          fd_link_x(i,3) = i-2
          fd_link_x(i,4) = i-1
        else if(i.eq.4) then
          fd_link_x(i,1) = nx
          fd_link_x(i,2) = i-3
          fd_link_x(i,3) = i-2
          fd_link_x(i,4) = i-1        
        else 
          fd_link_x(i,1) = i-4
          fd_link_x(i,2) = i-3
          fd_link_x(i,3) = i-2
          fd_link_x(i,4) = i-1        
        end if                        
        if(i.eq.nx) then
          fd_link_x(i,6) = 1
          fd_link_x(i,7) = 2
          fd_link_x(i,8) = 3
          fd_link_x(i,9) = 4
        else if(i.eq.nx-1) then
          fd_link_x(i,6) = i+1
          fd_link_x(i,7) = 1
          fd_link_x(i,8) = 2
          fd_link_x(i,9) = 3
        else if(i.eq.nx-2) then
          fd_link_x(i,6) = i+1
          fd_link_x(i,7) = i+2
          fd_link_x(i,8) = 1
          fd_link_x(i,9) = 2
        else if(i.eq.nx-3) then
          fd_link_x(i,6) = i+1
          fd_link_x(i,7) = i+2
          fd_link_x(i,8) = i+3
          fd_link_x(i,9) = 1        
        else 
          fd_link_x(i,6) = i+1
          fd_link_x(i,7) = i+2
          fd_link_x(i,8) = i+3
          fd_link_x(i,9) = i+4
        end if             
     end do
     
     do i=1,ny
        fd_link_y(i,5) = i
        if(i.eq.1) then
          fd_link_y(i,1) = ny-3
          fd_link_y(i,2) = ny-2
          fd_link_y(i,3) = ny-1
          fd_link_y(i,4) = ny          
        else if(i.eq.2) then
          fd_link_y(i,1) = ny-2
          fd_link_y(i,2) = ny-1
          fd_link_y(i,3) = ny
          fd_link_y(i,4) = i-1        
        else if(i.eq.3) then
          fd_link_y(i,1) = ny-1
          fd_link_y(i,2) = ny
          fd_link_y(i,3) = i-2
          fd_link_y(i,4) = i-1
        else if(i.eq.4) then
          fd_link_y(i,1) = ny
          fd_link_y(i,2) = i-3
          fd_link_y(i,3) = i-2
          fd_link_y(i,4) = i-1        
        else 
          fd_link_y(i,1) = i-4
          fd_link_y(i,2) = i-3
          fd_link_y(i,3) = i-2
          fd_link_y(i,4) = i-1        
        end if                        
        if(i.eq.ny) then
          fd_link_y(i,6) = 1
          fd_link_y(i,7) = 2
          fd_link_y(i,8) = 3
          fd_link_y(i,9) = 4
        else if(i.eq.ny-1) then
          fd_link_y(i,6) = i+1
          fd_link_y(i,7) = 1
          fd_link_y(i,8) = 2
          fd_link_y(i,9) = 3
        else if(i.eq.ny-2) then
          fd_link_y(i,6) = i+1
          fd_link_y(i,7) = i+2
          fd_link_y(i,8) = 1
          fd_link_y(i,9) = 2
        else if(i.eq.ny-3) then
          fd_link_y(i,6) = i+1
          fd_link_y(i,7) = i+2
          fd_link_y(i,8) = i+3
          fd_link_y(i,9) = 1        
        else 
          fd_link_y(i,6) = i+1
          fd_link_y(i,7) = i+2
          fd_link_y(i,8) = i+3
          fd_link_y(i,9) = i+4
        end if             
     end do
     
     !! Initialise FD weights
     call set_fd_weights
     
     !! Derivative weights for X ====================================
     do i=1,nx
        !! First derivative weights    
        fd_grad_x(i,:) =  fd_grad(:)*oodx

        !! Second derivative weights     
        fd_grad2_x(i,:) = fd_grad2(:)*oodx*oodx
        
        !! Filter weights
        fd_hyp_x(i,:) = fd_hyp(:)        

     end do
     
     !! Derivative weights for Y ====================================
     do j=1,ny
        !! First derivative weights      
        fd_grad_y(j,:) =  fd_grad(:)*oody

        !! Second derivative weights     
        fd_grad2_y(j,:) = fd_grad2(:)*oody*oody
        
        !! Filter weights
        fd_hyp_y(j,:) = fd_hyp(:)        

     end do


     !! Derivative weights for Y if Y is non-periodic
     fd_grad_y(1,:) = fd_grad_1(:)*oody
     fd_grad_y(2,:) = fd_grad_2(:)*oody
     fd_grad_y(3,:) = fd_grad_3(:)*oody
     fd_grad_y(4,:) = fd_grad_4(:)*oody               
     fd_grad_y(ny,:) = fd_grad_n(:)*oody
     fd_grad_y(ny-1,:) = fd_grad_nm1(:)*oody
     fd_grad_y(ny-2,:) = fd_grad_nm2(:)*oody
     fd_grad_y(ny-3,:) = fd_grad_nm3(:)*oody               

     fd_grad2_y(1,:) = fd_grad2_1(:)*oody*oody
     fd_grad2_y(2,:) = fd_grad2_2(:)*oody*oody
     fd_grad2_y(3,:) = fd_grad2_3(:)*oody*oody
     fd_grad2_y(4,:) = fd_grad2_4(:)*oody*oody          
     fd_grad2_y(ny,:) = fd_grad2_n(:)*oody*oody
     fd_grad2_y(ny-1,:) = fd_grad2_nm1(:)*oody*oody
     fd_grad2_y(ny-2,:) = fd_grad2_nm2(:)*oody*oody
     fd_grad2_y(ny-3,:) = fd_grad2_nm3(:)*oody*oody   
     
     fd_hyp_y(1,:) = fd_hyp_1(:)
     fd_hyp_y(2,:) = fd_hyp_2(:)
     fd_hyp_y(3,:) = fd_hyp_3(:)
     fd_hyp_y(4,:) = fd_hyp_4(:)
     fd_hyp_y(ny,:) = fd_hyp_n(:)
     fd_hyp_y(ny-1,:) = fd_hyp_nm1(:)
     fd_hyp_y(ny-2,:) = fd_hyp_nm2(:)
     fd_hyp_y(ny-3,:) = fd_hyp_nm3(:)                
     
     !! Derivative weights for X if X is non-periodic
     fd_grad_x(1,:) = fd_grad_1(:)*oodx
     fd_grad_x(2,:) = fd_grad_2(:)*oodx
     fd_grad_x(3,:) = fd_grad_3(:)*oodx
     fd_grad_x(4,:) = fd_grad_4(:)*oodx               
     fd_grad_x(nx,:) = fd_grad_n(:)*oodx
     fd_grad_x(nx-1,:) = fd_grad_nm1(:)*oodx
     fd_grad_x(nx-2,:) = fd_grad_nm2(:)*oodx
     fd_grad_x(nx-3,:) = fd_grad_nm3(:)*oodx               

     fd_grad2_x(1,:) = fd_grad2_1(:)*oodx*oodx
     fd_grad2_x(2,:) = fd_grad2_2(:)*oodx*oodx
     fd_grad2_x(3,:) = fd_grad2_3(:)*oodx*oodx
     fd_grad2_x(4,:) = fd_grad2_4(:)*oodx*oodx          
     fd_grad2_x(nx,:) = fd_grad2_n(:)*oodx*oodx
     fd_grad2_x(nx-1,:) = fd_grad2_nm1(:)*oodx*oodx
     fd_grad2_x(nx-2,:) = fd_grad2_nm2(:)*oodx*oodx
     fd_grad2_x(nx-3,:) = fd_grad2_nm3(:)*oodx*oodx              
     
     fd_hyp_x(1,:) = fd_hyp_1(:)
     fd_hyp_x(2,:) = fd_hyp_2(:)
     fd_hyp_x(3,:) = fd_hyp_3(:)
     fd_hyp_x(4,:) = fd_hyp_4(:)
     fd_hyp_x(nx,:) = fd_hyp_n(:)
     fd_hyp_x(nx-1,:) = fd_hyp_nm1(:)
     fd_hyp_x(nx-2,:) = fd_hyp_nm2(:)
     fd_hyp_x(nx-3,:) = fd_hyp_nm3(:)        
     
              
     deallocate(fd_grad,fd_grad2,fd_hyp)
     return
  end subroutine build_operators
!! ------------------------------------------------------------------------------------------------
  subroutine set_fd_weights
     integer(ikind) :: i,j,k

     !! Central 8th order for internal nodes
     allocate(fd_grad(9),fd_grad2(9),fd_hyp(9))
     !! First derivative weights    
     fd_grad(1) =  3.0d0/840.0d0
     fd_grad(2) = -32.0d0/840.0d0
     fd_grad(3) =  168.0d0/840.0d0
     fd_grad(4) = -672.0d0/840.0d0
     fd_grad(5) =  zero
     fd_grad(6) =  672.0d0/840.0d0
     fd_grad(7) = -168.0d0/840.0d0
     fd_grad(8) =  32.0d0/840.0d0
     fd_grad(9) = -3.0d0/840.0d0

     !! Second derivative weights     
     fd_grad2(1) = -9.0d0/5040.0d0
     fd_grad2(2) =  128.0d0/5040.0d0
     fd_grad2(3) = -1008.0d0/5040.0d0
     fd_grad2(4) =  8064.0d0/5040.0d0
     fd_grad2(5) = -14350.0d0/5040.0d0
     fd_grad2(6) =  8064.0d0/5040.0d0
     fd_grad2(7) = -1008.0d0/5040.0d0
     fd_grad2(8) =  128.0d0/5040.0d0
     fd_grad2(9) = -9.0d0/5040.0d0

     !! Filter weights
     fd_hyp(1) = -one/256.0d0
     fd_hyp(2) =  8.0d0/256.0d0
     fd_hyp(3) = -28.0d0/256.0d0
     fd_hyp(4) =  56.0d0/256.0d0
     fd_hyp(5) = -70.0d0/256.0d0
     fd_hyp(6) =  56.0d0/256.0d0
     fd_hyp(7) = -28.0d0/256.0d0
     fd_hyp(8) =  8.0d0/256.0d0
     fd_hyp(9) = -one/256.0d0   
     
     !! Symmetric 6th order for node 4
     allocate(fd_grad_4(9),fd_grad2_4(9),fd_hyp_4(9))
     !! First derivative weights    
     fd_grad_4(1) =  zero
     fd_grad_4(2) = -one/60.0d0
     fd_grad_4(3) =  9.0d0/60.0d0
     fd_grad_4(4) = -45.0d0/60.0d0
     fd_grad_4(5) =  zero
     fd_grad_4(6) =  45.0d0/60.0d0
     fd_grad_4(7) = -9.0d0/60.0d0
     fd_grad_4(8) =  one/60.0d0
     fd_grad_4(9) =  zero

     !! Second derivative weights     
     fd_grad2_4(1) =  zero
     fd_grad2_4(2) =  2.0d0/180.0d0
     fd_grad2_4(3) = -27.0d0/180.0d0
     fd_grad2_4(4) =  270.0d0/180.0d0
     fd_grad2_4(5) = -490.0d0/180.0d0
     fd_grad2_4(6) =  270.0d0/180.0d0
     fd_grad2_4(7) = -27.0d0/180.0d0
     fd_grad2_4(8) =  2.0d0/180.0d0
     fd_grad2_4(9) =  zero
     
     !! Filter weights
     fd_hyp_4(1) =  zero
     fd_hyp_4(2) =  one/64.0d0
     fd_hyp_4(3) = -6.0d0/64.0d0
     fd_hyp_4(4) =  15.0d0/64.0d0
     fd_hyp_4(5) = -20.0d0/64.0d0
     fd_hyp_4(6) =  15.0d0/64.0d0
     fd_hyp_4(7) = -6.0d0/64.0d0
     fd_hyp_4(8) =  one/64.0d0
     fd_hyp_4(9) =  zero     
     
     !! Symmetric 4th order for node 3
     allocate(fd_grad_3(9),fd_grad2_3(9),fd_hyp_3(9))
     !! First derivative weights    
     fd_grad_3(1) =  zero
     fd_grad_3(2) =  zero
     fd_grad_3(3) =  one/12.0d0
     fd_grad_3(4) = -8.0d0/12.0d0
     fd_grad_3(5) =  zero
     fd_grad_3(6) =  8.0d0/12.0d0
     fd_grad_3(7) = -one/12.0d0
     fd_grad_3(8) =  zero
     fd_grad_3(9) =  zero

     !! Second derivative weights     
     fd_grad2_3(1) =  zero
     fd_grad2_3(2) =  zero
     fd_grad2_3(3) = -one/12.0d0
     fd_grad2_3(4) =  16.0d0/12.0d0
     fd_grad2_3(5) = -30.0d0/12.0d0
     fd_grad2_3(6) =  16.0d0/12.0d0
     fd_grad2_3(7) = -one/12.0d0
     fd_grad2_3(8) =  zero
     fd_grad2_3(9) =  zero
     
     !! Filter weights
     fd_hyp_3(1) =  zero
     fd_hyp_3(2) =  zero
     fd_hyp_3(3) = -one/16.0d0
     fd_hyp_3(4) =  4.0d0/16.0d0
     fd_hyp_3(5) = -6.0d0/16.0d0
     fd_hyp_3(6) =  4.0d0/16.0d0
     fd_hyp_3(7) = -one/16.0d0
     fd_hyp_3(8) =  zero
     fd_hyp_3(9) =  zero      

     
     !! Assymmetric 4th order for node 2
     allocate(fd_grad_2(9),fd_grad2_2(9),fd_hyp_2(9))
     !! First derivative weights    
     fd_grad_2(1) =  zero
     fd_grad_2(2) =  zero
     fd_grad_2(3) =  zero
     fd_grad_2(4) = -3.0d0/12.0d0
     fd_grad_2(5) = -10.0d0/12.0d0
     fd_grad_2(6) =  18.0d0/12.0d0
     fd_grad_2(7) = -6.0d0/12.0d0
     fd_grad_2(8) =  one/12.0d0
     fd_grad_2(9) =  zero

     !! Second derivative weights     
     fd_grad2_2(1) =  zero
     fd_grad2_2(2) =  zero
     fd_grad2_2(3) =  zero
     fd_grad2_2(4) =  11.0d0/12.0d0
     fd_grad2_2(5) = -20.0d0/12.0d0
     fd_grad2_2(6) =  6.0d0/12.0d0
     fd_grad2_2(7) =  4.0d0/12.0d0
     fd_grad2_2(8) = -one/12.0d0
     fd_grad2_2(9) =  zero

     !! Filter weights
     fd_hyp_2(1) =  zero
     fd_hyp_2(2) =  zero
     fd_hyp_2(3) =  zero
     fd_hyp_2(4) = -one/16.0d0
     fd_hyp_2(5) =  4.0d0/16.0d0
     fd_hyp_2(6) = -6.0d0/16.0d0
     fd_hyp_2(7) =  4.0d0/16.0d0
     fd_hyp_2(8) = -one/16.0d0
     fd_hyp_2(9) =  zero      

     
     !! One-sided 4th order for node 1
     allocate(fd_grad_1(9),fd_grad2_1(9),fd_hyp_1(9))
     !! First derivative weights    
     fd_grad_1(1) =  zero
     fd_grad_1(2) =  zero
     fd_grad_1(3) =  zero
     fd_grad_1(4) =  zero
     fd_grad_1(5) = -25.0d0/12.0d0
     fd_grad_1(6) =  48.0d0/12.0d0
     fd_grad_1(7) = -36.0d0/12.0d0
     fd_grad_1(8) =  16.0d0/12.0d0
     fd_grad_1(9) = -3.0d0/12.0d0

     !! Second derivative weights     
     fd_grad2_1(1) =  zero
     fd_grad2_1(2) =  zero
     fd_grad2_1(3) =  zero
     fd_grad2_1(4) =  zero
     fd_grad2_1(5) =  35.0d0/12.0d0
     fd_grad2_1(6) = -104.0d0/12.0d0
     fd_grad2_1(7) =  114.0d0/12.0d0
     fd_grad2_1(8) = -56.0d0/12.0d0
     fd_grad2_1(9) =  11.0d0/12.0d0
     
     !! Filter weights
     fd_hyp_1(1) =  zero
     fd_hyp_1(2) =  zero
     fd_hyp_1(3) =  zero
     fd_hyp_1(4) =  zero
     fd_hyp_1(5) = -one/16.0d0
     fd_hyp_1(6) =  4.0d0/16.0d0
     fd_hyp_1(7) = -6.0d0/16.0d0
     fd_hyp_1(8) =  4.0d0/16.0d0
     fd_hyp_1(9) = -one/16.0d0     

     
     !! Reverse operators for n,n-1,n-2,n-3
     allocate(fd_grad_n(9),fd_grad2_n(9),fd_hyp_n(9))
     allocate(fd_grad_nm1(9),fd_grad2_nm1(9),fd_hyp_nm1(9))             
     allocate(fd_grad_nm2(9),fd_grad2_nm2(9),fd_hyp_nm2(9))             
     allocate(fd_grad_nm3(9),fd_grad2_nm3(9),fd_hyp_nm3(9))        
     
     !! nm3 and 4 are the same
     fd_grad_nm3(:) = fd_grad_4(:)
     fd_grad2_nm3(:) = fd_grad2_4(:)
     fd_hyp_nm3(:) = fd_hyp_4(:)          
     
     !! nm2 and 3 are the same
     fd_grad_nm2(:) = fd_grad_3(:)
     fd_grad2_nm2(:) = fd_grad2_3(:)
     fd_hyp_nm2(:) = fd_hyp_3(:)         
     
     !! Assymmetric 4th order for node nm1 is reversal of node 2
     !! First derivative weights   
     do k=1,9
        j=10-k
        fd_grad_nm1(j) = - fd_grad_2(k)  !! sign reversed
        fd_grad2_nm1(j) = fd_grad2_2(k)
        fd_hyp_nm1(j) = fd_hyp_2(k)        
     end do
     
     !! One-sided 4th order for node n is reversal of node 1
     do k=1,9
        j=10-k
        fd_grad_n(j) = - fd_grad_1(k)  !! sign reversed
        fd_grad2_n(j) = fd_grad2_1(k)
        fd_hyp_n(j) = fd_hyp_1(k)        
     end do
             
     return
  end subroutine set_fd_weights
!! ------------------------------------------------------------------------------------------------
end module weights
