module rhs
  !! This module contains routines to construct the RHS of all evolution equations
  !! These RHS' are built using calls to routines from the derivatives module, and
  !! calls to specific thermodynamics routines.
  
  !! Although separate subroutines, they must be called in correct order, as they rely
  !! on each other (e.g. divvel calculated in calc_rhs_lnro, also used in calc_rhs_vel)
  
  !! We use the lists internal_list and boundary_list to loop through internal and boundary
  !! nodes respectively. The L arrays for boundaries run from 1 to nb, and so use index j
  !! when within a loop over boundary nodes.
  use kind_parameters
  use common_parameter
  use common_vars
  use derivatives
  implicit none
  
 
  private
  public calc_rhs  

  !! Allocatable arrays for 1st and 2nd gradients
  real(rkind),dimension(:,:,:),allocatable :: gradw,gradT,gradvisc
  real(rkind),dimension(:,:),allocatable :: lapw,lapT,lapvisc,d2viscdxdy
  real(rkind),dimension(:,:),allocatable :: dvdymdudx !! dv/dy - du/dx
  
  

contains
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs
     !! Construct the RHS for w-equation
     integer(ikind) :: i,j
     real(rkind) :: coef_T,coef_w
      
     !! Allocate memory for spatial derivatives and stores
     allocate(lapw(nx,ny),lapT(nx,ny))
     allocate(gradw(nx,ny,dims))
     allocate(gradT(nx,ny,dims))

     !! Low order adiabatic thermal BCs on left and right walls
     T(1,1:ny) = T(2,1:ny)
     T(nx,1:ny) = T(nx-1,1:ny)

     !! Calculate spatial derivatives 
     call calc_laplacian(w,lapw)
     call calc_gradient(w,gradw)
     call calc_laplacian(T,lapT)
     call calc_gradient(T,gradT)

     !! Construct and solve Poisson equation to obtain velocity...
     call find_velocity
     
     !! Pre-set the viscosity to unity, then evaluate non-Newtonian viscosity
     visc(:,:) = one
     call evaluate_non_newtonian_viscosity

     !! Evaluate additional derivatives for non-Newtonian terms
     allocate(lapvisc(nx,ny),gradvisc(nx,ny,dims),d2viscdxdy(nx,ny))
     call calc_laplacian(visc,lapvisc)
     call calc_gradient(visc,gradvisc)
     call calc_ddy_of_ddx(gradvisc(:,:,1),d2viscdxdy)


     !! Set viscous and thermal diffusivity coefficients for w and T respectively
     coef_w = sqrt(Pr/Ra)
     coef_T = one/sqrt(Pr*Ra)
         
     !! Build RHS for internal nodes
     !$omp parallel do private(j)
     do i=1,nx
        do j=2,ny-1
         
           !! RHS of vorticity
           rhs_w(i,j) = -u(i,j)*gradw(i,j,1) - v(i,j)*gradw(i,j,2) + coef_w* &
                      ( visc(i,j)*lapw(i,j) + w(i,j)*lapvisc(i,j) &
                      + dot_product(gradw(i,j,:),gradvisc(i,j,:)) &
                      + d2viscdxdy(i,j)*dvdymdudx(i,j) )
                                                 

           !! RHS of Temperature
           rhs_T(i,j) = -u(i,j)*gradT(i,j,1) - v(i,j)*gradT(i,j,2) + lapT(i,j)*coef_T
           
           !! Add the temperature source term
           rhs_w(i,j) = rhs_w(i,j) + gradT(i,j,1)
        end do
     end do
     !$omp end parallel do
              
             
     !! Deallocate any stores no longer required
     deallocate(lapw,lapT)
     deallocate(gradw,gradT)     
     deallocate(dvdymdudx)
     deallocate(lapvisc,gradvisc,d2viscdxdy)
     
     
     !! Reapply vorticity BCs on walls
     rhs_w(:,1) = zero;rhs_w(:,ny) = zero    
     rhs_w(1,:) = zero;rhs_w(nx,:) = zero
     
     !! Reapply temperature BCs on walls
     rhs_T(:,1) = zero;rhs_T(:,ny) = zero
     rhs_T(1,:) = zero;rhs_T(nx,:) = zero

  
 
     return
  end subroutine calc_rhs
!! ------------------------------------------------------------------------------------------------  
  subroutine find_velocity 
     !! Evaluate streamfunction, then get velocity
     !! Currently using second order FD with successive over relaxation method, and low order 
     !! vorticity boundary conditions. 
     integer(ikind) :: i,j,k
     real(rkind),dimension(:,:,:),allocatable :: gradpsi
     real(rkind),dimension(:,:),allocatable :: psi_old
     real(rkind) :: x,y
     integer(ikind),parameter :: max_iters=100
     real(rkind),parameter :: err_tol = 1.0d-6
     integer(ikind) :: ip1,im1
     real(rkind) :: dx2,dy2,dx2pdy2,sor_factor,sor_err,ftn_psi
     logical :: keepgoing
  
     dx2 = dx*dx
     dy2 = dy*dy
     dx2pdy2 = dx2 + dy2
     sor_factor = 1.5d0   !! Over-relaxation factor
  
     allocate(psi_old(nx,ny))
   
     keepgoing = .true.;k=0
     !! Successive over relaxation
     do while(keepgoing)
        !! Store previous iter of psi in psi_old
        psi_old = psi

        !! Loop over x, excluding boundaries
        do i=2,nx-1
           !! Loop over y, excluding boundaries
           do j=2,ny-1
              ftn_psi = dy2*(psi(i+1,j)+psi(i-1,j)) + dx2*(psi(i,j+1)+psi(i,j-1)) + dx2*dy2*w(i,j)
              ftn_psi = half*ftn_psi/dx2pdy2
         
              psi(i,j) = sor_factor*ftn_psi + (one-sor_factor)*psi(i,j)
           end do
           
        end do
        
        !! Measure error
        sor_err=0.0
        do i=1,nx
           do j=1,ny
             sor_err=sor_err + (psi_old(i,j)-psi(i,j))**two
           end do
        end do     
        sor_err = sqrt(sor_err/dble(nx*ny))
        
        k=k+1
        !! Breakout criteria
        if(k.gt.max_iters) keepgoing=.false.
        if(sor_err.le.err_tol) keepgoing=.false.
        
     end do 
     
     !! Store data on SOR performance for screen output    
     itercount = k
     residual = sor_err
     
     !! Apply wall boundary conditions (low order) on vorticity
     w(1:nx,1) = -two*psi(1:nx,2)/dy2
     w(1:nx,ny) = -two*psi(1:nx,ny-1)/dy2
     w(1,1:ny) = -two*psi(2,1:ny)/dx2
     w(nx,1:ny) = -two*psi(nx-1,1:ny)/dx2

     !! Evaluate velocity on internal nodes
     allocate(gradpsi(nx,ny,dims))
     call calc_gradient(psi,gradpsi)
     !$omp parallel do private(j)
     do i=2,nx-1
        do j=2,ny-1
           u(i,j) = gradpsi(i,j,2)
           v(i,j) = -gradpsi(i,j,1)
        end do
     end do
     !$omp end parallel do
     deallocate(gradpsi)    
     
     !! Re-enforce u=v=0 on wall bounds
     u(1,:)=zero;v(1,:)=zero
     u(nx,:)=zero;v(nx,:)=zero
     u(:,1)=zero;v(:,1)=zero
     u(:,ny)=zero;v(:,ny)=zero               
  
     return
  end subroutine find_velocity  
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_non_newtonian_viscosity
     integer(ikind) :: i,j,k
     real(rkind),dimension(:,:,:),allocatable :: gradu,gradv
     real(rkind) :: norm_srt
     
     !! Evaluate velocity gradients
     allocate(gradu(nx,ny,dims),gradv(nx,ny,dims))
     call calc_gradient(u,gradu)
     call calc_gradient(v,gradv)
     
     !! Allocate dv/dy - du/dx which is calculated here and used later for non-Newtonian terms
     allocate(dvdymdudx(nx,ny))
     
     !! Loop over all nodes
     !$omp parallel do private(j,norm_srt)
     do i=1,nx
        do j=1,ny
        
           !! Evaluate the norm of the strain rate tensor
           norm_srt = sqrt( two*gradu(i,j,1)**two &
                               + (gradu(i,j,2)+gradv(i,j,1))**two &
                               + two*gradv(i,j,2)**two )
                               
           !! Set the viscosity according some Non-Newtonian model
           visc(i,j) = non_newtonian_viscosity(norm_srt)
           
           !! Evaluate dv/dy - du/dx
           dvdymdudx(i,j) = gradv(i,j,2) - gradu(i,j,1)
                               
alpha_out(i,j) = visc(i,j)
        end do
     end do
     !$omp end parallel do
     

     !! Release memory
     deallocate(gradu,gradv)

     return
  end subroutine evaluate_non_newtonian_viscosity
!! ------------------------------------------------------------------------------------------------
  function non_newtonian_viscosity(norm_srt) result(viscosity)
     !! This is the function which takes the norm of the strain rate tensor and returns a viscosity

     !! Comment out the models you're not using, and add more as required!!
     real(rkind),intent(in) :: norm_srt
     real(rkind) :: viscosity      
     
     !! Simple Newtonian model
     viscosity = one
     
     !! Power law
!     viscosity = (one + norm_srt)**(power_law_index - one)

     return
  end function non_newtonian_viscosity
!! ------------------------------------------------------------------------------------------------
end module rhs
