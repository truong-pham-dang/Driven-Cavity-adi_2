C****************************************************************************C
C                                                                            C
C     This program solve the Navier-Stokes (NS) equations using the          C
C     Alternating Direction Implicit (ADI) method. The solution is spatially C
C     second order accurate. This code is written by Prof. Ercan Erturk.     C
C     Visit http://www.cavityflow.com                                        C
C                                                                            C
C********************************************C*******************************C
C                                            C
C     s(i,j) ==> streamfunction variable     C
C     v(i,j) ==> vorticity variable          C
C     x(i)   ==> x-coordinate                C
C     y(j)   ==> y-coordinate                C
C     dh     ==> grid spacing                C
C     Re     ==> Renolds Number              C
C                                            C
C********************************************C
      
      program main
      
      implicit double precision (a-h,o-z)
      
      parameter(N=128)
      
      common / flow variables /
     > s(0:N,0:N),v(0:N,0:N),
     > s_old(0:N,0:N),v_old(0:N,0:N)
      common / geometry /
     > x(0:N),y(0:N)
      common / matrix elements /
     > aa(0:N),bb(0:N),cc(0:N),dd(0:N)
      common / other /
     > rhs1(0:N,0:N),rhs2(0:N,0:N),
     > sx(0:N,0:N),sy(0:N,0:N)


       Re=1000.d0

       dh=1.0d0/dble(N)

      do 1 k=0,N
       x(k)=dble(k)/dble(N)
       y(k)=dble(k)/dble(N)
 1    continue

C     Time step
       alpha=0.6d0
       dt=alpha*dh*dh

C     Initial guess
c     Note that when homogeneous initial guess is used, in the first iteration 
c     residual_3 is indeterminate. In order to avoid this, instead of using 
c     exactly zero initial values, at interior points use very very small 
c     numbers that could be considered as zero.
      do 2 k=0,N
      do 102 j=0,N
       s(k,j)=1.d-32
       v(k,j)=1.d-32
 102  continue
       s(0,k)=0.d0
       v(0,k)=0.d0
       s(N,k)=0.d0
       v(N,k)=0.d0
       s(k,0)=0.d0
       v(k,0)=0.d0
       s(k,N)=0.d0
       v(k,N)=0.d0
 2    continue

C     Record the CPU time at start
       call cpu_time(time_start)

C     Start the iterations
      do 999 iteration=1,1000000

C     Update old variables
      do 5 i=1,N-1
      do 5 j=1,N-1
       s_old(i,j)=s(i,j)
       v_old(i,j)=v(i,j)
 5    continue

C     SOLVE THE STREAMFUNCTION EQUATION        

c     Implicit in x-direction
c     Calculate the RHS
      do 10 i=1,N-1
      do 10 j=1,N-1
       rhs1(i,j)=s(i,j)
     >          +0.5d0*dt*(s(i,j-1)-2.d0*s(i,j)+s(i,j+1))/dh**2.
     >          +0.5d0*dt*v(i,j)        
 10   continue

      do 11 j=1,N-1
      do 12 i=1,N-1
       aa(i)=-0.5d0*dt/dh**2.
       bb(i)=1.d0+0.5d0*dt*2.d0/dh**2.
       cc(i)=-0.5d0*dt/dh**2.
       dd(i)=rhs1(i,j)
 12   continue

c     Note: s(0,j)=0 and s(N,j)=0
c     Forward elimination
      do 13 i=2,N-1
       bb(i)=bb(i)-aa(i)*cc(i-1)/bb(i-1)
       dd(i)=dd(i)-aa(i)*dd(i-1)/bb(i-1)
 13   continue

c     Substitute for the last point
       s(N-1,j)=dd(N-1)/bb(N-1)

c     Backward substitution
      do 14 i=N-2,1,-1
       s(i,j)=(dd(i)-s(i+1,j)*cc(i))/bb(i)
 14   continue
 11   continue

c     Implicit in y-direction
c     Calculate the RHS
      do 20 i=1,N-1
      do 20 j=1,N-1
       rhs2(i,j)=s(i,j)
     >          +0.5d0*dt*(s(i-1,j)-2.d0*s(i,j)+s(i+1,j))/dh**2.
     >          +0.5d0*dt*v(i,j)        
 20   continue

      do 21 i=1,N-1
      do 22 j=1,N-1
       aa(j)=-0.5d0*dt/dh**2.
       bb(j)=1.d0+0.5d0*dt*2.d0/dh**2.
       cc(j)=-0.5d0*dt/dh**2.
       dd(j)=rhs2(i,j)
 22   continue

c     Note: s(i,0)=0 and s(i,N)=0
c     Forward elimination
      do 23 j=2,N-1
       bb(j)=bb(j)-aa(j)*cc(j-1)/bb(j-1)
       dd(j)=dd(j)-aa(j)*dd(j-1)/bb(j-1)
 23   continue

c     Substitute for the last point
       s(i,N-1)=dd(N-1)/bb(N-1)

c     Backward substitution
      do 24 j=N-2,1,-1
       s(i,j)=(dd(j)-s(i,j+1)*cc(j))/bb(j)
 24   continue
 21   continue

C     SOLVE THE VORTICITY EQUATION

C     Calculate vorticity at the wall
C     NOTE:For these boundary conditions please refer to:
C          T. Stortkuhl, C. Zenger, S. ZN-1mer, "An Asymptotic Solution for
C          the Singularity at the Angular Point of the Lid Driven Cavity",
C          International Journal of Numerical Methods for Heat & Fluid Flow
C          1994, Vol 4, pp 47--59
      do 30 k=1,N-1
       v(k,0)=(
     >       -(s(k-1,1)+s(k,1)+s(k+1,1))/(3.d0*dh**2.)
     >       -(0.5d0*v(k-1,0)+0.5d0*v(k+1,0)
     >       +0.25d0*v(k-1,1)+v(k,1)+0.25d0*v(k+1,1))/(9.d0)
     >        )*9.d0/2.d0
       v(k,N)=(-1.d0/dh
     >       -(s(k-1,N-1)+s(k,N-1)+s(k+1,N-1))/(3.d0*dh**2.)
     >       -(0.5d0*v(k-1,N)+0.5d0*v(k+1,N)
     >       +0.25d0*v(k-1,N-1)+v(k,N-1)+0.25d0*v(k+1,N-1))/(9.d0)
     >        )*9.d0/2.d0
       v(0,k)=(
     >       -(s(1,k-1)+s(1,k)+s(1,k+1))/(3.d0*dh**2.)
     >       -(0.5d0*v(0,k-1)+0.5d0*v(0,k+1)
     >       +0.25d0*v(1,k-1)+v(1,k)+0.25d0*v(1,k+1))/(9.d0)
     >        )*9.d0/2.d0
       v(N,k)=(
     >       -(s(N-1,k-1)+s(N-1,k)+s(N-1,k+1))/(3.d0*dh**2.)
     >       -(0.5d0*v(N,k-1)+0.5d0*v(N,k+1)
     >       +0.25d0*v(N-1,k-1)+v(N-1,k)+0.25d0*v(N-1,k+1))/(9.d0)
     >        )*9.d0/2.d0
 30   continue
       v(0,0)=(-(s(1,1))/(3.d0*dh**2.)-(0.5d0*v(1,0)+0.5d0*v(0,1)
     >       +0.25d0*v(1,1))/(9.d0))*9.d0
       v(N,0)=(-(s(N-1,1))/(3.d0*dh**2.)-(0.5d0*v(N-1,0)+0.5d0*v(N,1)
     >       +0.25d0*v(N-1,1))/(9.d0))*9.d0
       v(N,N)=(-0.5d0/dh-(s(N-1,N-1))/(3.d0*dh**2.)-(0.5d0*v(N-1,N)
     >       +0.5d0*v(N,N-1)+0.25d0*v(N-1,N-1))/(9.d0))*9.d0
       v(0,N)=(-0.5d0/dh-(s(1,N-1))/(3.d0*dh**2.)-(0.5d0*v(1,N)
     >       +0.5d0*v(0,N-1)+0.25d0*v(1,N-1))/(9.d0))*9.d0

      do 31 i=1,N-1
      do 31 j=1,N-1
       sx(i,j)=Re*(s(i+1,j)-s(i-1,j))/(2.d0*dh)
       sy(i,j)=Re*(s(i,j+1)-s(i,j-1))/(2.d0*dh)
 31   continue

c     Implicit in x-direction
c     Calculate the RHS
      do 40 i=1,N-1
      do 40 j=1,N-1
       rhs1(i,j)=v(i,j)
     >          +0.5d0*dt*(
     >               (v(i,j-1)-2.d0*v(i,j)+v(i,j+1))/dh**2.
     >               +sx(i,j)*(v(i,j+1)-v(i,j-1))/(2.d0*dh)
     >                    )
 40   continue

      do 41 j=1,N-1
      do 42 i=1,N-1
       aa(i)=-0.5d0*dt/dh**2.-0.5d0*dt*sy(i,j)/(2.d0*dh)
       bb(i)=1.d0+0.5d0*dt*2.d0/dh**2.
       cc(i)=-0.5d0*dt/dh**2.+0.5d0*dt*sy(i,j)/(2.d0*dh)
       dd(i)=rhs1(i,j)
 42   continue

c     Apply BC
       dd(1)=dd(1)-aa(1)*v(0,j)
       dd(N-1)=dd(N-1)-cc(N-1)*v(N,j)

c     Forward elimination
      do 43 i=2,N-1
       bb(i)=bb(i)-aa(i)*cc(i-1)/bb(i-1)
       dd(i)=dd(i)-aa(i)*dd(i-1)/bb(i-1)
 43   continue

c     Substitute for the last point
       v(N-1,j)=dd(N-1)/bb(N-1)

c     Backward substitution
      do 44 i=N-2,1,-1
       v(i,j)=(dd(i)-v(i+1,j)*cc(i))/bb(i)
 44   continue
 41   continue

c     Implicit in y-direction
c     Calculate the RHS
      do 50 i=1,N-1
      do 50 j=1,N-1
       rhs2(i,j)=v(i,j)
     >          +0.5d0*dt*(
     >               (v(i-1,j)-2.d0*v(i,j)+v(i+1,j))/dh**2.
     >               -sy(i,j)*(v(i+1,j)-v(i-1,j))/(2.d0*dh)
     >                    )
 50   continue

      do 51 i=1,N-1
      do 52 j=1,N-1
       aa(j)=-0.5d0*dt/dh**2.+0.5d0*dt*sx(i,j)/(2.d0*dh)
       bb(j)=1.d0+0.5d0*dt*2.d0/dh**2.
       cc(j)=-0.5d0*dt/dh**2.-0.5d0*dt*sx(i,j)/(2.d0*dh)
       dd(j)=rhs2(i,j)
 52   continue

c     Apply BC
       dd(1)=dd(1)-aa(1)*v(i,0)
       dd(N-1)=dd(N-1)-cc(N-1)*v(i,N)

c     Forward elimination
      do 53 j=2,N-1
       bb(j)=bb(j)-aa(j)*cc(j-1)/bb(j-1)
       dd(j)=dd(j)-aa(j)*dd(j-1)/bb(j-1)
 53   continue

c     Substitute for the last point
       v(i,N-1)=dd(N-1)/bb(N-1)

c     Backward substitution
      do 54 j=N-2,1,-1
       v(i,j)=(dd(j)-v(i,j+1)*cc(j))/bb(j)
 54   continue
 51   continue


C     Check the residuals to see if convergence is achieved. 
C     You can comment out all or some, for a faster run.

c     residual_1 is the residual of the governing equations.
       residual_1_s_A=0.d0
       residual_1_v_A=0.d0
c     residual_2 is the change in variables (indicates the significant digit) in a time step.
       residual_2_s_A=0.d0
       residual_2_v_A=0.d0
c     residual_3 is the normalized change in variables (indicates percent change) in a time step.
       residual_3_s_A=0.d0
       residual_3_v_A=0.d0
      
      do 60 i=1,N-1
      do 60 j=1,N-1
       residual_1_s_B=abs(
     >   (s(i-1,j)-2.d0*s(i,j)+s(i+1,j))/dh**2.
     >  +(s(i,j-1)-2.d0*s(i,j)+s(i,j+1))/dh**2.
     >  +v(i,j)          )
       residual_1_v_B=abs(
     >   (1.d0/Re)*(v(i-1,j)-2.d0*v(i,j)+v(i+1,j))/dh**2.
     >  +(1.d0/Re)*(v(i,j-1)-2.d0*v(i,j)+v(i,j+1))/dh**2.
     >  -(s(i,j+1)-s(i,j-1))/(2.d0*dh)*(v(i+1,j)-v(i-1,j))/(2.d0*dh)
     >  +(s(i+1,j)-s(i-1,j))/(2.d0*dh)*(v(i,j+1)-v(i,j-1))/(2.d0*dh)
     >                   )

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
 60   continue
      
C     Output the residuals
       write(*,*) ' '
       write(*,*) iteration
       write(*,*) residual_1_s_A,residual_1_v_A
       write(*,*) residual_2_s_A,residual_2_v_A
       write(*,*) residual_3_s_A,residual_3_v_A
c     condition to stop iterations
       if((residual_1_s_A.lt.1.d-6).AND.(residual_1_v_A.lt.1.d-6))
     > goto 1000
      
 999  continue
      
C     Record the CPU time at finish
 1000 call cpu_time(time_finish)

       write(*,*) ' '
       write(*,*) 'Convergence is achieved in',iteration,'   iterations'
       write(*,*) 'CPU time=',time_finish-time_start

c     Output to a file        
       open(1,file='out.txt')
      do i=0,N
      do j=0,N
       write(1,*) x(i),y(j),s(i,j),v(i,j)
      end do
       write(1,*) ' '
      end do
 1112 continue
c 2222 format(f,x,f,x,es,x,es)

      stop
      end
