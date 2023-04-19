module NUMERICAL_METHODS

   implicit none

contains

   subroutine poisson_2d_mixed_bc (nx,ny,dx,dy,f,gb,gt,u)

      !2-D Poisson solver for mixed boundary conditions in a rectangular
      !domain: 1) periodic in the x-direction; 2) inhomogeneous Neumann at
      !bottom and top
   
      use singleton, only: fft
      
      implicit none
   
      integer :: nx,ny
      real*8  :: dx,dy
      
      real*8  :: u(0:nx-1,1:ny)
      real*8  :: f(0:nx-1,1:ny)
      real*8  :: gb(0:nx-1),gt(0:nx-1)
      
      real*8, parameter  :: pi = 4.0d0*datan(1.0d0)
      
      real*8  :: uk_r(0:nx-1,1:ny),uk_i(0:nx-1,1:ny)
      real*8  :: fk_r(0:nx-1,1:ny),fk_i(0:nx-1,1:ny)
      real*8  :: gbk_r(0:nx-1),gbk_i(0:nx-1)
      real*8  :: gtk_r(0:nx-1),gtk_i(0:nx-1) 
   
      real*8  :: alpha_k
      real*8  :: a(1:ny),b(1:ny),c(1:ny),r_r(1:ny),r_i(1:ny),x(1:ny)
   
      complex :: cx(0:nx-1),ck(0:nx-1)
   
      real*8  :: sk1
      
      integer :: i,j,k
      
      !Find the "f" coefficients in Fourier space
      do j = 1,ny
         do i = 0,nx-1
            cx(i) = cmplx(f(i,j),0.0)
         end do
         
         ck = fft(cx)
         
         do i = 0,nx-1
            fk_r(i,j) =  real(ck(i))    !real part
            fk_i(i,j) = aimag(ck(i))    !imaginary part
         end do
      end do
   
      !Find the "gb" coefficients in Fourier space
      do i = 0,nx-1
         cx(i) = cmplx(gb(i),0.0)
      end do
      
      ck = fft(cx)
      
      do i = 0,nx-1
         gbk_r(i) =  real(ck(i))  !real part
         gbk_i(i) = aimag(ck(i))  !imaginary part
      end do
      
      !Find the "gt" coefficients in Fourier space
      do i = 0,nx-1
         cx(i) = cmplx(gt(i),0.0)
      end do
   
      ck = fft(cx)
      
      do i = 0,nx-1
         gtk_r(i) =  real(ck(i))  !real part
         gtk_i(i) = aimag(ck(i))  !imaginary part
      end do
      
      !**********************************************************************
      !Solve the system of equations for "uk_r" and "uk_i"
   
      !----------------------------------------------------------------------
      !For NONZERO k
      
      do k = 1,nx-1
         !Vectors "a", "b", and "c" of the tridiagonal matrix as inputs for "tridag"
         a(1) = 0.d0 !undefined (not used by tridag)
         a(2:ny-1) = 1.d0
         a(ny) = 2.d0
   
         alpha_k = (dy**2/dx**2)*( 2.d0*dcos( 2.d0*pi*dfloat(k)/dfloat(nx) ) - 2.d0 ) - 2.d0
         
         b = alpha_k
         
         c(1) = 2.d0
         c(2:ny-1) = 1.d0
         c(ny) = 0.d0 !undefined (not used by tridag)
   
         !Right-hand side
         r_r(1) = dy**2 * fk_r(k,1) + 2.d0*dy*gbk_r(k)
         r_i(1) = dy**2 * fk_i(k,1) + 2.d0*dy*gbk_i(k)
         do j = 2,ny-1
            r_r(j) = dy**2 * fk_r(k,j)
            r_i(j) = dy**2 * fk_i(k,j)
         end do
         r_r(ny) = dy**2 * fk_r(k,ny) - 2.d0*dy*gtk_r(k)
         r_i(ny) = dy**2 * fk_i(k,ny) - 2.d0*dy*gtk_i(k)
         
         !tridiagonal solver
         call tridag(a,b,c,r_r,x,ny) !real part
         uk_r(k,1:ny) = x(1:ny)
         
         call tridag(a,b,c,r_i,x,ny) !imaginary part
         uk_i(k,1:ny) = x(1:ny)
      end do
   
      !---------------------------------------------------------------------
      !For k = 0
      u(0,1) = 0.d0       !Solvability condition
      
      uk_i(0,1:ny) = 0.d0 !Coefficients for k = 0 are real-valued. The
                        !real parts are evaluated below
      
      !For j = 1
      sk1 = 0.d0
      do k = 1,nx-1
         sk1 = sk1 + uk_r(k,1)
      end do
      uk_r(0,1) = u(0,1) - sk1
      
      !For j = 2
      uk_r(0,2) = uk_r(0,1) + dy*gbk_r(0) + 0.5d0*dy**2 * fk_r(0,1)
   
      !For j = 3, ..., ny-1
      do j = 3,ny-1
         uk_r(0,j) = 2.d0*uk_r(0,j-1) - uk_r(0,j-2) + dy**2 * fk_r(0,j-1)
      end do
   
      !For j = ny
      uk_r(0,ny) = uk_r(0,ny-1) + dy*gtk_r(0) - 0.5d0*dy**2 * fk_r(0,ny)
      
      !Find "u" values on physical space
      !forward fourier transform
      do j = 1,ny
         do k = 0,nx-1
            ck(k) = cmplx( uk_r(k,j), uk_i(k,j) )
         end do
   
         cx = fft(ck,inv=.true.)
   
         do i = 0,nx-1
            u(i,j) = real( cx(i) )
         end do
      end do
      
   end subroutine poisson_2d_mixed_bc
 
   subroutine poisson_2d_mixed_bc_xper_yneu_stagg (nx,ny,dx,dy,f,gb,gt,u)
   
      !2-D Poisson solver for mixed boundary conditions in a rectangular
      !domain: 1) periodic in the x-direction; 2) inhomogeneous Neumann at
      !bottom and top in a STAGGERED grid
   
      use singleton, only: fft
      
      implicit none
   
      integer :: nx,ny
      real*8  :: dx,dy
      
      real*8  :: u(0:nx-1,1:ny)
      real*8  :: f(0:nx-1,1:ny)
      real*8  :: gb(0:nx-1),gt(0:nx-1)
      
      real*8, parameter  :: pi = 4.0d0*datan(1.0d0)
      
      real*8  :: uk_r(0:nx-1,1:ny),uk_i(0:nx-1,1:ny)
      real*8  :: fk_r(0:nx-1,1:ny),fk_i(0:nx-1,1:ny)
      real*8  :: gbk_r(0:nx-1),gbk_i(0:nx-1)
      real*8  :: gtk_r(0:nx-1),gtk_i(0:nx-1) 
   
      real*8  :: alpha_k
      real*8  :: a(1:ny),b(1:ny),c(1:ny),r_r(1:ny),r_i(1:ny),x(1:ny)
   
      complex :: cx(0:nx-1),ck(0:nx-1)
   
      real*8  :: sk1
      
      integer :: i,j,k
      
      !Find the "f" coefficients in Fourier space
      do j = 1,ny
         do i = 0,nx-1
            cx(i) = cmplx(f(i,j),0.0)
         end do
         
         ck = fft(cx)
         
         do i = 0,nx-1
            fk_r(i,j) =  real(ck(i))    !real part
            fk_i(i,j) = aimag(ck(i))    !imaginary part
         end do
      end do
   
      !Find the "gb" coefficients in Fourier space
      do i = 0,nx-1
         cx(i) = cmplx(gb(i),0.0)
      end do
      
      ck = fft(cx)
      
      do i = 0,nx-1
         gbk_r(i) =  real(ck(i))  !real part
         gbk_i(i) = aimag(ck(i))  !imaginary part
      end do
      
      !Find the "gt" coefficients in Fourier space
      do i = 0,nx-1
         cx(i) = cmplx(gt(i),0.0)
      end do
   
      ck = fft(cx)
      
      do i = 0,nx-1
         gtk_r(i) =  real(ck(i))  !real part
         gtk_i(i) = aimag(ck(i))  !imaginary part
      end do
      
      !**********************************************************************
      !Solve the system of equations for "uk_r" and "uk_i"
   
      !For NONZERO k  
      do k = 1,nx-1
         
         !Vectors "a", "b", and "c" of the tridiagonal matrix as inputs for "tridag"
         a(1) = 0.d0 !undefined (not used by tridag)
         a(2:ny-1) = 1.d0
         a(ny) = 1.d0 !Neumann - staggered (NS)
   
         alpha_k = (dy**2/dx**2)*( 2.d0*dcos( 2.d0*pi*dfloat(k)/dfloat(nx) ) - 2.d0 ) - 2.d0
   
         b( 1)     = 1.d0 + alpha_k !NS
         b(ny)     = 1.d0 + alpha_k !NS
         b(2:ny-1) =        alpha_k !NS
         
         c(1) = 1.d0 !NS
         c(2:ny-1) = 1.d0
         c(ny) = 0.d0 !undefined (not used by tridag)
   
         !Right-hand side
         r_r(1) = dy**2 * fk_r(k,1) + dy*gbk_r(k) !NS
         r_i(1) = dy**2 * fk_i(k,1) + dy*gbk_i(k) !NS
         do j = 2,ny-1
            r_r(j) = dy**2 * fk_r(k,j)
            r_i(j) = dy**2 * fk_i(k,j)
         end do
         r_r(ny) = dy**2 * fk_r(k,ny) - dy*gtk_r(k) !NS
         r_i(ny) = dy**2 * fk_i(k,ny) - dy*gtk_i(k) !NS
   
         !tridiagonal solver     
         call tridag(a,b,c,r_r,x,ny) !real part
         uk_r(k,1:ny) = x(1:ny)
         
         call tridag(a,b,c,r_i,x,ny) !imaginary part
         uk_i(k,1:ny) = x(1:ny)
      end do
   
      !For k = 0
      u(0,1) = 0.d0       !Solvability condition
      
      uk_i(0,1:ny) = 0.d0 !Coefficients for k = 0 are real-valued. The
                        !real parts are evaluated below
      
      !For j = 1  
      sk1 = 0.d0
      do k = 1,nx-1
         sk1 = sk1 + uk_r(k,1)
      end do
      uk_r(0,1) = u(0,1) - sk1
      
      !For j = 2
      uk_r(0,2) = uk_r(0,1) + dy*gbk_r(0) + dy**2 * fk_r(0,1) !NS
   
      !For j = 3, ..., ny-1
      do j = 3,ny-1
         uk_r(0,j) = 2.d0*uk_r(0,j-1) - uk_r(0,j-2) + dy**2 * fk_r(0,j-1)
      end do
   
      !For j = ny
      uk_r(0,ny) = uk_r(0,ny-1) + dy*gtk_r(0) - dy**2 * fk_r(0,ny) !NS
      
      !Find "u" values on physical space (inverse fourier transform)
      do j = 1,ny
         do k = 0,nx-1
            ck(k) = cmplx( uk_r(k,j), uk_i(k,j))
         end do
   
         cx = fft(ck,inv=.true.)
   
         do i = 0,nx-1
            u(i,j) = real( cx(i) )
         end do
      end do
   
   end subroutine poisson_2d_mixed_bc_xper_yneu_stagg

   subroutine tridag(a,b,c,r,u,n)

      implicit none

      integer :: n
      real*8  :: a(n),b(n),c(n),r(n),u(n)

      !Solves for a vector u(1:n) of length n the tridiagonal linear
      !ser. a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors and are
      !not modified.

      integer :: j
      real*8  :: bet,gam(n) !One vector of workspace, gam is needed.

      if (b(1).eq.0.d0) then
         write(*,*) "tridag: rewrite equations"
         call exit
         !If this happens then you should rewrite your equations as a set
         !of order N-1, with u_2 trivially eliminated.

      end if

      bet = b(1)
      u(1) = r(1) / bet

      do j = 2,n !Decomposition and forward substitution.
         gam(j) = c(j-1) / bet
         bet = b(j) - a(j)*gam(j)
         if (bet.eq.0.d0) then
            write(*,*) "tridag failed"
            call exit
            !Algorithm fails; see below.
         end if
         u(j) = ( r(j) - a(j)*u(j-1) ) / bet
      end do

      do j = n-1,1,-1 !Back substitution.
         u(j) = u(j) - gam(j+1)*u(j+1)
      end do
   end subroutine tridag
END MODULE NUMERICAL_METHODS