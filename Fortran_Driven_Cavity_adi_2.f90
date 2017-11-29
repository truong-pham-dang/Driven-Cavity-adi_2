!****************************************************************************C
!                                                                            C
!     This program solve the Navier-Stokes (NS) equations using the          C
!     Alternating Direction Implicit (ADI) method. The solution is spatially C
!     second order accurate. This code is written by Prof. Ercan Erturk.     C
!     Visit http://www.cavityflow.com                                        C
!                                                                            C
!********************************************C*******************************C
!                                            C
!     s(i,j) ==> streamfunction variable     C
!     v(i,j) ==> vorticity variable          C
!     x(i)   ==> x-coordinate                C
!     y(j)   ==> y-coordinate                C
!     dh     ==> grid spacing                C
!     Re     ==> Renolds Number              C
!                                            C
!********************************************C
    
!     Updated:   27/11/2017
!     Purpose:   Porting from FORTRAN 77 to Fortran 90
!     Developer: DANG Truong
      
      program main
      !$ use omp_lib
      implicit none
      
      integer, parameter :: dp = SELECTED_REAL_KIND(15,307)
      
      integer, parameter :: N=128
      
      real(dp) :: s(0:N,0:N), v(0:N,0:N), s_old(0:N,0:N), v_old(0:N,0:N)
      real(dp) :: x(0:N), y(0:N)
      real(dp) :: aa(0:N),bb(0:N),cc(0:N),dd(0:N)
      real(dp) :: rhs1(0:N,0:N), rhs2(0:N,0:N), sx(0:N,0:N), sy(0:N,0:N)
      
      real(dp) :: Re, dh
      
      real(dp) :: time_start, time_finish
      
      real(dp) :: time_start_omp, time_finish_omp
      
      real(dp) :: alpha, dt
      
      real(dp) :: residual_1_s_B, residual_1_v_B, residual_1_v_A, residual_2_s_A
      real(dp) :: residual_1_s_A, residual_2_s_B, residual_2_v_A, residual_2_v_B
      real(dp) :: residual_3_s_A, residual_3_s_B, residual_3_v_A, residual_3_v_B
      
      integer :: k, i, j
      integer :: iteration
      
      integer, parameter :: nthreads = 2
      
      CALL OMP_SET_NUM_THREADS(nthreads)


       Re=1000.d0

       dh=1.0d0/dble(N)

      do k=0,N
       x(k)=dble(k)/dble(N)
       y(k)=dble(k)/dble(N)
      enddo 

!     Time step
       alpha=0.6d0
       dt=alpha*dh*dh

!     Initial guess
!     Note that when homogeneous initial guess is used, in the first iteration 
!     residual_3 is indeterminate. In order to avoid this, instead of using 
!     exactly zero initial values, at interior points use very very small 
!     numbers that could be considered as zero.
      
      do k=0,N
      do j=0,N
       s(k,j)=1.d-32
       v(k,j)=1.d-32
      enddo
       s(0,k)=0.d0
       v(0,k)=0.d0
       s(N,k)=0.d0
       v(N,k)=0.d0
       s(k,0)=0.d0
       v(k,0)=0.d0
       s(k,N)=0.d0
       v(k,N)=0.d0
      enddo 
      

!     Record the CPU time at start
       call cpu_time(time_start)
       
       time_start_omp = omp_get_wtime()

!     Start the iterations
      do iteration=1,1000000

!     Update old variables
      !$omp parallel 
      do i=1,N-1
      do j=1,N-1
       s_old(i,j)=s(i,j)
       v_old(i,j)=v(i,j)
      enddo 
      enddo 
      !$omp end parallel
      
      

!     SOLVE THE STREAMFUNCTION EQUATION        

!     Implicit in x-direction
!     Calculate the RHS
      
      do i=1,N-1
      do j=1,N-1
       rhs1(i,j)=s(i,j) + 0.5d0*dt*(s(i,j-1)-2.d0*s(i,j)+s(i,j+1))/dh**2. + 0.5d0*dt*v(i,j)        
      enddo 
      enddo 

      do j=1,N-1
      do i=1,N-1
       aa(i)=-0.5d0*dt/dh**2.
       bb(i)=1.d0+0.5d0*dt*2.d0/dh**2.
       cc(i)=-0.5d0*dt/dh**2.
       dd(i)=rhs1(i,j)
      enddo 

!     Note: s(0,j)=0 and s(N,j)=0
!     Forward elimination
      do i=2,N-1
       bb(i)=bb(i)-aa(i)*cc(i-1)/bb(i-1)
       dd(i)=dd(i)-aa(i)*dd(i-1)/bb(i-1)
      enddo 

!     Substitute for the last point
       s(N-1,j)=dd(N-1)/bb(N-1)

!     Backward substitution
      do i=N-2,1,-1
       s(i,j)=(dd(i)-s(i+1,j)*cc(i))/bb(i)
      enddo 
      
      enddo 

!     Implicit in y-direction
!     Calculate the RHS
      do i=1,N-1
      do j=1,N-1
       rhs2(i,j)=s(i,j) + 0.5d0*dt*(s(i-1,j)-2.d0*s(i,j)+s(i+1,j))/dh**2. +0.5d0*dt*v(i,j)        
      enddo 
      enddo 

      do i=1,N-1
      do j=1,N-1
       aa(j)=-0.5d0*dt/dh**2.
       bb(j)=1.d0+0.5d0*dt*2.d0/dh**2.
       cc(j)=-0.5d0*dt/dh**2.
       dd(j)=rhs2(i,j)
      enddo 

!     Note: s(i,0)=0 and s(i,N)=0
!     Forward elimination
      do j=2,N-1
       bb(j)=bb(j)-aa(j)*cc(j-1)/bb(j-1)
       dd(j)=dd(j)-aa(j)*dd(j-1)/bb(j-1)
      enddo 

!     Substitute for the last point
       s(i,N-1)=dd(N-1)/bb(N-1)

!     Backward substitution
      do j=N-2,1,-1
       s(i,j)=(dd(j)-s(i,j+1)*cc(j))/bb(j)
      enddo 
      
      enddo 

!     SOLVE THE VORTICITY EQUATION

!     Calculate vorticity at the wall
!     NOTE:For these boundary conditions please refer to:
!          T. Stortkuhl, C. Zenger, S. ZN-1mer, "An Asymptotic Solution for
!          the Singularity at the Angular Point of the Lid Driven Cavity",
!          International Journal of Numerical Methods for Heat & Fluid Flow
!          1994, Vol 4, pp 47--59
      do k=1,N-1
       v(k,0)=( &
     &       -(s(k-1,1)+s(k,1)+s(k+1,1))/(3.d0*dh**2.) &
     &       -(0.5d0*v(k-1,0)+0.5d0*v(k+1,0) &
     &       +0.25d0*v(k-1,1)+v(k,1)+0.25d0*v(k+1,1))/(9.d0) &
     &        )*9.d0/2.d0
       v(k,N)=(-1.d0/dh &
     &       -(s(k-1,N-1)+s(k,N-1)+s(k+1,N-1))/(3.d0*dh**2.) &
     &       -(0.5d0*v(k-1,N)+0.5d0*v(k+1,N) &
     &       +0.25d0*v(k-1,N-1)+v(k,N-1)+0.25d0*v(k+1,N-1))/(9.d0) &
     &        )*9.d0/2.d0
       v(0,k)=( &
     &       -(s(1,k-1)+s(1,k)+s(1,k+1))/(3.d0*dh**2.) &
     &       -(0.5d0*v(0,k-1)+0.5d0*v(0,k+1) &
     &       +0.25d0*v(1,k-1)+v(1,k)+0.25d0*v(1,k+1))/(9.d0) &
     &        )*9.d0/2.d0
       v(N,k)=( &
     &       -(s(N-1,k-1)+s(N-1,k)+s(N-1,k+1))/(3.d0*dh**2.) &
     &       -(0.5d0*v(N,k-1)+0.5d0*v(N,k+1) &
     &       +0.25d0*v(N-1,k-1)+v(N-1,k)+0.25d0*v(N-1,k+1))/(9.d0) &
     &        )*9.d0/2.d0
     enddo 
       v(0,0)=(-(s(1,1))/(3.d0*dh**2.)-(0.5d0*v(1,0)+0.5d0*v(0,1) &
     &       +0.25d0*v(1,1))/(9.d0))*9.d0
       v(N,0)=(-(s(N-1,1))/(3.d0*dh**2.)-(0.5d0*v(N-1,0)+0.5d0*v(N,1) &
     &       +0.25d0*v(N-1,1))/(9.d0))*9.d0
       v(N,N)=(-0.5d0/dh-(s(N-1,N-1))/(3.d0*dh**2.)-(0.5d0*v(N-1,N) &
     &       +0.5d0*v(N,N-1)+0.25d0*v(N-1,N-1))/(9.d0))*9.d0
       v(0,N)=(-0.5d0/dh-(s(1,N-1))/(3.d0*dh**2.)-(0.5d0*v(1,N) &
     &       +0.5d0*v(0,N-1)+0.25d0*v(1,N-1))/(9.d0))*9.d0

      do i=1,N-1
      do j=1,N-1
       sx(i,j)=Re*(s(i+1,j)-s(i-1,j))/(2.d0*dh)
       sy(i,j)=Re*(s(i,j+1)-s(i,j-1))/(2.d0*dh)
      enddo 
      enddo 

!     Implicit in x-direction
!     Calculate the RHS
      do i=1,N-1
      do j=1,N-1
       rhs1(i,j)=v(i,j) &
     &          +0.5d0*dt*( &
     &               (v(i,j-1)-2.d0*v(i,j)+v(i,j+1))/dh**2. &
     &               +sx(i,j)*(v(i,j+1)-v(i,j-1))/(2.d0*dh) &
     &                    )
      enddo 
      enddo 

      do j=1,N-1
      do i=1,N-1
       aa(i)=-0.5d0*dt/dh**2.-0.5d0*dt*sy(i,j)/(2.d0*dh)
       bb(i)=1.d0+0.5d0*dt*2.d0/dh**2.
       cc(i)=-0.5d0*dt/dh**2.+0.5d0*dt*sy(i,j)/(2.d0*dh)
       dd(i)=rhs1(i,j)
      enddo 

!     Apply BC
       dd(1)=dd(1)-aa(1)*v(0,j)
       dd(N-1)=dd(N-1)-cc(N-1)*v(N,j)

!     Forward elimination
      do i=2,N-1
       bb(i)=bb(i)-aa(i)*cc(i-1)/bb(i-1)
       dd(i)=dd(i)-aa(i)*dd(i-1)/bb(i-1)
      enddo 

!     Substitute for the last point
       v(N-1,j)=dd(N-1)/bb(N-1)

!     Backward substitution
      do i=N-2,1,-1
       v(i,j)=(dd(i)-v(i+1,j)*cc(i))/bb(i)
      enddo 
      enddo

!     Implicit in y-direction
!     Calculate the RHS
      do i=1,N-1
      do j=1,N-1
       rhs2(i,j)=v(i,j) &
     &          +0.5d0*dt*( &
     &               (v(i-1,j)-2.d0*v(i,j)+v(i+1,j))/dh**2. &
     &               -sy(i,j)*(v(i+1,j)-v(i-1,j))/(2.d0*dh) &
     &                    )
      enddo 
      enddo 

      do i=1,N-1
      do j=1,N-1
       aa(j)=-0.5d0*dt/dh**2.+0.5d0*dt*sx(i,j)/(2.d0*dh)
       bb(j)=1.d0+0.5d0*dt*2.d0/dh**2.
       cc(j)=-0.5d0*dt/dh**2.-0.5d0*dt*sx(i,j)/(2.d0*dh)
       dd(j)=rhs2(i,j)
      enddo 

!     Apply BC
       dd(1)=dd(1)-aa(1)*v(i,0)
       dd(N-1)=dd(N-1)-cc(N-1)*v(i,N)

!     Forward elimination
      do j=2,N-1
       bb(j)=bb(j)-aa(j)*cc(j-1)/bb(j-1)
       dd(j)=dd(j)-aa(j)*dd(j-1)/bb(j-1)
      enddo 

!     Substitute for the last point
       v(i,N-1)=dd(N-1)/bb(N-1)

!     Backward substitution
      do j=N-2,1,-1
       v(i,j)=(dd(j)-v(i,j+1)*cc(j))/bb(j)
      enddo 
      enddo 


!     Check the residuals to see if convergence is achieved. 
!     You can comment out all or some, for a faster run.

!     residual_1 is the residual of the governing equations.
       residual_1_s_A=0.d0
       residual_1_v_A=0.d0
!     residual_2 is the change in variables (indicates the significant digit) in a time step.
       residual_2_s_A=0.d0
       residual_2_v_A=0.d0
!     residual_3 is the normalized change in variables (indicates percent change) in a time step.
       residual_3_s_A=0.d0
       residual_3_v_A=0.d0
      
      do i=1,N-1
      do j=1,N-1
       residual_1_s_B=abs( &
     &   (s(i-1,j)-2.d0*s(i,j)+s(i+1,j))/dh**2. &
     &  +(s(i,j-1)-2.d0*s(i,j)+s(i,j+1))/dh**2. &
     &  +v(i,j)          )
       residual_1_v_B=abs( &
     &   (1.d0/Re)*(v(i-1,j)-2.d0*v(i,j)+v(i+1,j))/dh**2. &
     &  +(1.d0/Re)*(v(i,j-1)-2.d0*v(i,j)+v(i,j+1))/dh**2. &
     &  -(s(i,j+1)-s(i,j-1))/(2.d0*dh)*(v(i+1,j)-v(i-1,j))/(2.d0*dh) &
     &  +(s(i+1,j)-s(i-1,j))/(2.d0*dh)*(v(i,j+1)-v(i,j-1))/(2.d0*dh) &
     &                   )

       residual_2_s_B=abs(s(i,j)-s_old(i,j))
       residual_2_v_B=abs(v(i,j)-v_old(i,j))

       residual_3_s_B=abs((s(i,j)-s_old(i,j))/s_old(i,j))
       residual_3_v_B=abs((v(i,j)-v_old(i,j))/v_old(i,j))

       residual_1_s_A=max(residual_1_s_A,residual_1_s_B)
       residual_1_v_A=max(residual_1_v_A,residual_1_v_B)

       residual_2_s_A=max(residual_2_s_A,residual_2_s_B)
       residual_2_v_A=max(residual_2_v_A,residual_2_v_B)

       residual_3_s_A=max(residual_3_s_A,residual_3_s_B)
       residual_3_v_A=max(residual_3_v_A,residual_3_v_B)
      enddo 
      enddo 
      
!     Output the residuals
       write(*,*) ' '
       write(*,*) iteration
       write(*,*) residual_1_s_A,residual_1_v_A
       write(*,*) residual_2_s_A,residual_2_v_A
       write(*,*) residual_3_s_A,residual_3_v_A
!     condition to stop iterations
       if((residual_1_s_A.lt.1.d-6).AND.(residual_1_v_A.lt.1.d-6)) then
        goto 1000
       endif
       
      enddo 
      
!     Record the CPU time at finish
1000 call cpu_time(time_finish)
     time_finish_omp = omp_get_wtime()

       write(*,*) ' '
       write(*,*) 'Convergence is achieved in',iteration,'   iterations'
       write(*,*) 'CPU time=',time_finish-time_start
       write(*,*) 'OpenMP time =',time_finish_omp - time_start_omp
       
       open(2, file ='stat.txt', status = 'replace')
       write(2,*) 'Convergence is achieved in',iteration,'   iterations'
       write(2,*) 'CPU time=',time_finish-time_start
       write(2,*) 'OpenMP time =',time_finish_omp - time_start_omp
       close(2)

!     Output to a file        
      open(unit = 1,file='out.txt')
      do i=0,N
      do j=0,N
       write(1,2222) x(i),y(j),s(i,j),v(i,j)
      enddo 
      enddo 
      close(1)
      
 2222 format(f,x,f,x,es,x,es)

      stop
      end
