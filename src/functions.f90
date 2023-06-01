module FUNCTIONS
   use CONSTANTS
   implicit none

   interface interpol
      module procedure interpol_scalar, interpol_array
   end interface

   interface SVP
         module procedure SVP_scalar, SVP_array, SVP_matrix
   end interface

   interface SAMPLE_GAUSSIAN
      module procedure SAMPLE_GAUSSIAN_SCALAR, SAMPLE_GAUSSIAN_ARRAY
   end interface

contains

   subroutine progress_bar(current_time,maximum_time)

      use IO_PARAMETERS, only: output_file

      real*8, intent(in) :: current_time, maximum_time
      real*8             :: frac
      integer            :: i, sticks
      character          :: bar(52) !50 points and two brackets
      real*8   , save    :: frac_old = 0.0
      integer  , save    :: start = 0, old = 0, rate = 0
      integer            :: now
      real*8             :: rem_time, total_time, elapsed_time

      if (start == 0) then
         call system_clock(count_rate=rate)
         call system_clock(start)
         old = start
         return
      end if

      frac = current_time/maximum_time

      if (frac > 1.00001D0 .or. frac < 0.D0) then
         print*,'Fraction should be between 0 and 1'
         return
      end if

      bar(1) = '['
      sticks = nint(frac*(size(bar)-2))
      do i = 2,size(bar)-1
         if (i <= sticks+1) then
               bar(i) = '|'
         else
               bar(i) = '.'
         end if
      end do
      bar(size(bar)) = ']'

      call system_clock(now)
      elapsed_time = (real(now - start,8)/real(rate,8))
      rem_time     = (real(now -   old,8)/real(rate,8))/(frac-frac_old)*(1 - frac)
      total_time   = elapsed_time + rem_time

      old = now
      frac_old = frac

      call system('clear')
      print'(A,F0.2,A)', 'Progress: ', frac*100,'%'
      print'(52A)'     ,  bar
      print'(A,A)'     , 'Elapsed time        : ', PRETTY_TIME(elapsed_time)
      print'(A,A)'     , 'Remaining time      : ', PRETTY_TIME(rem_time    )
      print'(A,A)'     , 'Estimated total time: ', PRETTY_TIME(total_time  )
      print'(A,F0.2)'  , 'Performance factor  : ', maximum_time/total_time
      print'(A,A)'     , 'Output file         : ', trim(output_file)

   end subroutine progress_bar

   function PRETTY_TIME(time_input) result(time_string)

      real*8, intent(in)    :: time_input
      integer      :: seconds, minutes, hours
      character(9) :: time_string

      hours   = nint(time_input)/3600
      minutes = nint(time_input)/60 - hours*60
      seconds = mod(nint(time_input),60)

      write(time_string,'(I0.2,A,I2.2,A,I2.2)') hours,':',minutes,':',seconds

   end function PRETTY_TIME

   subroutine BOX_MULLER (G1,G2)
      !Samples up to two standard gaussian variables
      implicit none

      real*8 :: W, G1, X1, X2, RAND(2)
      real*8, optional :: G2

      do
         call RANDOM_NUMBER(RAND)
         X1 = 2.D0*RAND(1)-1.D0
         X2 = 2.D0*RAND(2)-1.D0
         W  = X1**2 + X2**2
         if (W < 1.D0) exit
      end do

      W = sqrt( -2.D0*LOG(W)/W )

      G1 = X1*W
      if (present(G2)) G2 = X2*W

   end subroutine BOX_MULLER

   subroutine random_integer(min, max, r)

      ! Samples a random integer in the interval [min, max]
      ! with uniform distribution

      integer, intent(in)  :: min, max
      integer, intent(out) :: r
      real*8 :: r_float
   
      call RANDOM_NUMBER(r_float)
      r = floor(r_float*(max - min + 1) + min)
   
   end subroutine

   subroutine SAMPLE_GAUSSIAN_SCALAR(X,sigma,mu)
      implicit none
      real*8, intent(in) :: sigma, mu
      real*8  :: X

      call BOX_MULLER(X)
      X = X*sigma + mu

   end subroutine SAMPLE_GAUSSIAN_SCALAR

   subroutine SAMPLE_GAUSSIAN_ARRAY(X,sigma,mu)
      implicit none
      real*8, intent(in) :: sigma, mu
      real*8  :: X(:)
      integer :: i

      do i = 1,size(X)
         call SAMPLE_GAUSSIAN_SCALAR(X(i),sigma,mu)
      end do

   end subroutine SAMPLE_GAUSSIAN_ARRAY

   function OPEN_NEW_FILE(filename) result(ID)
   implicit none
      character(len=*), intent(in) :: filename
      integer, save :: count = 1
      integer :: ID

      open(UNIT = count, FILE = filename)
      ID = count
      count = count + 1;

   end function OPEN_NEW_FILE

   function linspace(a,b,N) result(V)
      real*8   , intent(in)   :: a, b
      integer  , intent(in)   :: N
      real*8                  :: V(N), dv
      integer                 :: i

      if (N .GT. 1) then
         dv = (b-a)/(N-1)
         do i = 1,N
               V(i) = a + dv*(i-1)
         end do
      else
         print*, 'N must be greater than 1.'
         V = 0
      end if
   end function linspace

   function interpol_scalar(coarse_x,coarse_y,fine_x,linear_x) result(fine_y)
      
      real*8 , intent(in) :: coarse_x(:), coarse_y(:), fine_x
      logical, intent(in) :: linear_x
      real*8              :: fine_y
      integer             :: j, N
      real*8              :: DX

      if (size(coarse_x) /= size(coarse_y)) then
         print*, 'Coarse arrays must be of same length.'
         print*, 'Size X', size(coarse_x)
         print*, 'Size Y', size(coarse_y)
         stop
      end if

      N = size(coarse_x)

      if (linear_x) then
         DX = (coarse_x(N) - coarse_x(1))/(N-1)
         j = min(ceiling((fine_x - coarse_x(1)) / DX) + 1, N)
         j = max(j,2)
      else 
         j = 2
         do while (coarse_x(j) <= fine_x .AND. j < N)
            j = j + 1
         end do
      end if

      fine_y = (coarse_y(j) - coarse_y(j-1))/(coarse_x(j) - coarse_x(j-1)) &
               *(fine_x - coarse_x(j-1)) + coarse_y(j-1)
   end function interpol_scalar

   function interpol_array(coarse_x, coarse_y, fine_x, linear_x) result(fine_y)
      real*8 , intent(in) :: coarse_x(:), coarse_y(:), fine_x(:)
      logical, intent(in) :: linear_x
      real*8              :: fine_y(size(fine_x))
      integer             :: i

      do i = 1,size(fine_x)
         fine_y(i) = interpol_scalar(coarse_x,coarse_y,fine_x(i),linear_x)
      end do

   end function interpol_array

   function SVP_scalar(T,phase) result(es)
      use constants
      real*8 :: es
      real*8, intent(in) :: T
      character(1), intent(in) :: phase

      if (phase == 'L') then
         es = ES0*exp(-LV0/R_V*(1/T - 1/T0))
      else if (phase == 'I') then
         es = ES0*exp(-LS0/R_V*(1/T - 1/T0))
      else
         stop 'SVP_SCALAR: Invalid phase tag.'
      end if
   end function SVP_scalar

   function SVP_array(T,phase) result(es)
      character(1), intent(in)  :: phase
      real*8      , intent(in)  :: T(:)
      real*8                    :: es(size(T))
      integer                   :: i

      do i = 1,size(T)
         es(i) = SVP_scalar(T(i),phase)
      end do
   end function SVP_array

   function SVP_matrix(T,phase) result(es)
      character(1), intent(in)  :: phase
      real*8      , intent(in)  :: T(:,:)
      integer                   :: dims(2)
      real*8, allocatable       :: es(:,:)
      integer                   :: i

      dims = shape(T)
      allocate(es(dims(1),dims(2)))
      do i = 1,dims(1)
         es(i,:) = SVP_array(T(i,:),phase)
      end do
   end function SVP_matrix

   function r_vs(T,P,phase) result(r)
      use constants
      real*8 :: r
      real*8, intent(in) :: T,P
      real*8 :: es
      character(1), intent(in) :: phase

      es = SVP(T,phase)
      r = (R_D/R_V)*es/(P - es)

   end function r_vs

   elemental function EXNER(P) result(EX)
      use constants
      real*8, intent(in) :: P
      real*8 :: EX

      EX = (P/P00)**(R_D/CP_D)
   end function EXNER

   function PERN(p,NN) result(pp)
      integer   :: p, pp, NN
      ! Periodic Nodes (different from PER(p,N), which works on cells)
      if (p < 1) then
         pp = p + (NN - 1)
      elseif (p > NN) then
         pp = p - (NN - 1)
      else
         pp = p
      end if
   end function

   subroutine DIV_SYM_TENSOR_2D(DX,DY,Axx,Axy,Ayy,divA_x,divA_y)
      real*8, intent(in ) :: DX, DY, Axx(:,:), Axy(:,:), Ayy(:,:)
      real*8, intent(out) :: divA_x(size(Axx,1),size(Axx,2))
      real*8, intent(out) :: divA_y(size(Axx,1),size(Axx,2))
      real*8              :: DDXx, DDXy, DDYx, DDYy
      integer             :: i, j, N, M
      logical             :: A, B
      !Size check
      A = size(Axx,1)==size(Axy,1).and.size(Axx,1)==size(Ayy,1)
      B = size(Axx,2)==size(Axy,2).and.size(Axx,2)==size(Ayy,2)
      if (A.and.B) then
         N = size(Axx,1)
         M = size(Axx,2)
      else
         stop "DIV_SYM_TENSOR_2D: Input matrices must have the same size."
      end if

      !Divergence calculation
      do i = 1,N
         do j = 1,M
            !Periodic boundaries (sides)
            DDXx = (Axx(PERN(i+1,N),j) - Axx(PERN(i-1,N),j))/(2*DX)
            DDXy = (Axy(PERN(i+1,N),j) - Axy(PERN(i-1,N),j))/(2*DX)
            if (j == 1) then !Lower boundary
               DDYx = (-3*Axy(i,j) + 4*Axy(i,j+1) - Axy(i,j+2))/(2*DY)
               DDYy = (-3*Ayy(i,j) + 4*Ayy(i,j+1) - Ayy(i,j+2))/(2*DY)
            elseif (j == M) then !Upper boundary
               DDYx = ( 3*Axy(i,j) - 4*Axy(i,j-1) + Axy(i,j-2))/(2*DY)
               DDYy = ( 3*Ayy(i,j) - 4*Ayy(i,j-1) + Ayy(i,j-2))/(2*DY)
            else
               DDYx = (Axy(i,j+1) - Axy(i,j-1))/(2*DY)
               DDYy = (Ayy(i,j+1) - Ayy(i,j-1))/(2*DY)
            end if
            divA_x(i,j) = DDXx + DDYx
            divA_y(i,j) = DDXy + DDYy
         end do
      end do
   end subroutine DIV_SYM_TENSOR_2D

   subroutine GRAD_VECTOR_2D(DX,DY,Ux,Uy,grad)
      real*8, intent(in ) :: DX, DY, Ux(:,:), Uy(:,:)
      real*8, intent(out) :: grad(size(Ux,1),size(Ux,2),2,2)
      integer             :: i, j, N, M

      !Size check
      if (size(Ux,1)==size(Uy,1).and.size(Ux,2)==size(Uy,2)) then
         N = size(Ux,1)
         M = size(Ux,2)
      else
         stop "GRAD_VECTOR_2D: Input matrices must have the same size."
      end if

      !Divergence calculation
      do i = 1,N
         do j = 1,M
            !Periodic boundaries (sides)
            grad(i,j,1,1) = (Ux(PERN(i+1,N),j) - Ux(PERN(i-1,N),j))/(2*DX)
            grad(i,j,1,2) = (Uy(PERN(i+1,N),j) - Uy(PERN(i-1,N),j))/(2*DX)
            if (j == 1) then !Lower boundary
               grad(i,j,2,1) = (-3*Ux(i,j) + 4*Ux(i,j+1) - Ux(i,j+2))/(2*DY)
               grad(i,j,2,2) = (-3*Uy(i,j) + 4*Uy(i,j+1) - Uy(i,j+2))/(2*DY)
            elseif (j == M) then !Upper boundary
               grad(i,j,2,1) = ( 3*Ux(i,j) - 4*Ux(i,j-1) + Ux(i,j-2))/(2*DY)
               grad(i,j,2,2) = ( 3*Uy(i,j) - 4*Uy(i,j-1) + Uy(i,j-2))/(2*DY)
            else
               grad(i,j,2,1) = (Ux(i,j+1) - Ux(i,j-1))/(2*DY)
               grad(i,j,2,2) = (Uy(i,j+1) - Uy(i,j-1))/(2*DY)
            end if

         end do
      end do
   end subroutine GRAD_VECTOR_2D

   subroutine DIV_VECTOR_2D(DX,DY,U,div)
      real*8, intent(in ) :: DX, DY, U(:,:,:)
      real*8, intent(out) :: div(size(U,1),size(U,2))
      real*8              :: DUx, DUy
      integer             :: i, j, N, M

      !Size check
      if (size(U,3)==2) then
         N = size(U,1)
         M = size(U,2)
      else
         stop "DIV_VECTOR_2D: The third dimension of U must be 2."
      end if

      !Divergence calculation
      do i = 1,N
         do j = 1,M
            !Periodic boundaries (sides)
            DUx = (U(PERN(i+1,N),j,1) - U(PERN(i-1,N),j,1))/(2*DX)
            if (j == 1) then !Lower boundary
               DUy = (-3*U(i,j,2) + 4*U(i,j+1,2) - U(i,j+2,2))/(2*DY)
            elseif (j == M) then !Upper boundary
               DUy = ( 3*U(i,j,2) - 4*U(i,j-1,2) + U(i,j-2,2))/(2*DY)
            else
               DUy = (U(i,j+1,2) - U(i,j-1,2))/(2*DY)
            end if
            div(i,j) = DUx + DUy
         end do
      end do
   end subroutine DIV_VECTOR_2D

end module FUNCTIONS

module CCN
   implicit none
   real*8, allocatable :: CDF_TF(:,:)  !CDF for the freezing temperature
   real*8, allocatable :: CDF_DR(:,:)  !CDF for the Dry radii
   private :: SET_CDF_DRY, SET_CDF_TFREEZING, INIT_RANDOM_SEED

contains

   subroutine SET_CDF_DRY
      use FUNCTIONS

      integer   :: N_CDF
      real*8    :: r_cutoff, mu, sig
      !Building the CDF table
      N_CDF = 1000
      r_cutoff = 5.D-2    !Cut off radius in micrometers
      mu = log(7.5D-2)    !Mean
      sig = log(1.6D0)    !Standard deviation

      allocate(CDF_DR(2,N_CDF))
      CDF_DR(1,:) = linspace(r_cutoff,5.D-1,N_CDF)
      CDF_DR(2,:) = 0.5D0 + 0.5D0*erf((log(CDF_DR(1,:)) - mu)/(sqrt(2.D0)*sig))
   end subroutine SET_CDF_DRY

   subroutine SAMPLE_DRY(R_DRY)
      use FUNCTIONS

      integer :: i
      real*8, intent(out)  :: R_DRY(:)

      call SET_CDF_DRY
      call INIT_RANDOM_SEED

      do i = 1,size(R_DRY)
         call random_number(R_DRY(i))
         do while (R_DRY(i) < CDF_DR(2,1))
            call random_number(R_DRY(i))
         end do
      end do
      R_DRY = interpol(CDF_DR(2,:),CDF_DR(1,:),R_DRY,linear_x = .false.)
      R_DRY = R_DRY*1D-6 !Conversion from micro meters to meters
   end subroutine SAMPLE_DRY


   subroutine SET_CDF_TFREEZING

      !Sets the cumulative distribution function (CDF)
      !for the freezing temperature.
      use constants, only: pi

      real*8  :: T,TMIN,TMAX,DT
      real*8  :: TC,TCMIN,TCMAX
      real*8  :: A,NS,DNS_DT
      real*8  :: PDF
      real*8  :: NSF,X
      integer :: I,M

      NSF(X) = EXP( - 0.517D0*X + 8.934D0 ) ![m^{-2}]
      !---------------------------------------------------------------
      !Assume the radius of the insolube IN equal to 0.5 micron
      !---------------------------------------------------------------
      A = 4.D0*PI*(5.D-7)**2
      !---------------------------------------------------------------

      TCMIN = -36.D0
      TCMAX = -12.D0

      TMIN = -40.D0 + 273.15D0
      TMAX =   0.D0 + 273.15D0

      M = 1000

      DT = (TMAX-TMIN)/REAL(M)

      allocate( CDF_TF(0:1,0:M) )

      CDF_TF(0,0) = TMIN
      CDF_TF(1,0) = 0.D0
      PDF = 0.D0

      do I = 1,M
         T = TMIN + DT*I
         TC = T - 273.15D0
         if (TC.GT.TCMAX) then
            NS = 0.D0
            DNS_DT = 0.D0
         else if ( (TC.LE.TCMAX).AND.(TC.GT.TCMIN) ) then
            NS = NSF(TC)
            DNS_DT = -0.517D0*NS
         else if ( TC.LT.TCMIN ) then
            NS = NSF(TCMIN)
            DNS_DT = 0.D0
         end if

         CDF_TF(0,I) = T
         CDF_TF(1,I) = EXP(-A*NS)

         PDF = -A*DNS_DT*EXP(-A*NS)
      end do
   end subroutine SET_CDF_TFREEZING


   subroutine SAMPLE_TFREEZING (TF, N)

      !Samples "N" values of freezing temperatures (stored in TF(:)) using
      !the inverse transform sampling method. It uses the pre-calculated
      !cumulative distribution function (CDF_TF).

      use FUNCTIONS, only: interpol, random_integer
      

      real*8, allocatable, intent(out) :: TF(:)
      integer            , intent(in ) :: N

      real*8  :: R, P
      integer :: k, N_insol
      integer :: current_size

      print'(A)', 'Randomizing Freezing temperatures. This might take a while...'
      call SET_CDF_TFREEZING
      call INIT_RANDOM_SEED

      P = 0.1 !Percentage of SDs with rd_insol
      N_insol = ceiling(P*N)
      allocate(TF(N))

      current_size = N - N_insol
      TF(1:current_size) = -40.D0 + 273.15D0
      
      do while (current_size < N)
         call random_integer(1, current_size, k)
         call RANDOM_NUMBER(R)
         R = interpol(CDF_TF(1,:), CDF_TF(0,:), R, linear_x=.false.)

         TF(current_size+1) = TF(k)
         TF(k) = R
         current_size = current_size + 1
      end do

   end subroutine SAMPLE_TFREEZING

   subroutine INIT_RANDOM_SEED

      implicit none

      integer :: I,N,CLOCK
      integer, allocatable :: SEED(:)

      call RANDOM_SEED(SIZE = N)
      allocate(SEED(N))
      call SYSTEM_CLOCK(COUNT=CLOCK)
      SEED = CLOCK + 37 * (/ (I - 1, I = 1, N) /)
      call RANDOM_SEED(PUT = SEED)
      deallocate(SEED)
   end subroutine INIT_RANDOM_SEED

end module CCN

module CIC_ROUTINES
   use CUSTOM_TYPES
   use GRID


   implicit none

contains

   function bilinear_weight(r) result(b)
      real*8, intent(in) :: r(2)
      real*8 :: b

      b = max(0.D0 , 1.D0 - abs(r(1))/DX) * max(0.D0 , 1.D0 - abs(r(1))/DZ)
   end function

   function periodic_box(i,j) result(b)
      integer   :: i, j
      type(box_address) :: b

      b = box_address(i,j)

      if (b%j < 1 .or. b%j > NZ) then
         ! Null indeces for boxes below ground or above ceiling
         b%i = 0
         b%j = 0
      else
         ! Periodic B.C. on side walls
         if (b%i < 1) then
            b%i = b%i + NX
         elseif (b%i > NX) then
            b%i = b%i - NX
         end if
      end if
   end function

   subroutine INIT_CIC_NODES

      integer   :: i, j

      if (.not.allocated(cell_nodes)) allocate(cell_nodes(NXP,NZP))

      do i = 1,NXP
         do j = 1,NZP
            !Defining the 4 neighbouring boxes around each node
            cell_nodes(i,j)%srd_box(1) = periodic_box(i  ,j  )
            cell_nodes(i,j)%srd_box(2) = periodic_box(i-1,j  )
            cell_nodes(i,j)%srd_box(3) = periodic_box(i-1,j-1)
            cell_nodes(i,j)%srd_box(4) = periodic_box(i  ,j-1)

            !Defining the position of each node in space (m)
            cell_nodes(i,j)%pos = [DX*(i-1), DZ*(j-1)]
         end do
      end do
   end subroutine

   function CIC_SCALAR_AT_NODES(LAGR_SCALAR) result(SCALAR_FIELD)
      use SD_VARIABLES, only: N_sd, X
      use CUSTOM_TYPES
      use OMP_LIB
      use STEPPER, only: TIME
      !Debugging only
      use THERMODYNAMIC_VAR, only: DPB
      !use SD_VARIABLES, only: SD_box

      real*8, intent(in)  :: LAGR_SCALAR(N_sd)
      real*8                 :: SCALAR_FIELD(NXP,NZP)
      real*8                 :: NN, r(2)
      integer,   allocatable :: srd_particles(:)
      real*8 ,   allocatable :: bw(:)

      integer   :: i, j, k, p, q, m, n, length

      !Debugging only
      character(len=30) :: DPB_file_name


      !$OMP parallel do private(j, k, p, q, m, n, length, srd_particles, r, bw, NN)
      do i = 1,NXP
         do j = 1,NZP
            !Get the amount of particles in surrounding boxes
            length = 0;
            do k = 1,4
               m = cell_nodes(i,j)%srd_box(k)%i
               n = cell_nodes(i,j)%srd_box(k)%j
               if (m/=0 .and. n/=0) length = length + size(boxes(m,n)%p_list)
            end do

            !Reset the surrounding particles list
            if (allocated(srd_particles)) deallocate(srd_particles)
            allocate(srd_particles(length))
            
            !Build a list of all surrounding particles
            p = 1
            do k = 1,4
               m = cell_nodes(i,j)%srd_box(k)%i
               n = cell_nodes(i,j)%srd_box(k)%j
               if (m/=0 .and. n/=0) then
                  do q = 1,size(boxes(m,n)%p_list)
                     srd_particles(p) = boxes(m,n)%p_list(q)
                     p = p + 1
                  end do
               end if
            end do
            
            !Reset the weights list
            if (allocated(bw)) deallocate(bw)
            allocate(bw(size(srd_particles)))

            !Calculate bw and NN
            do k = 1,size(bw)
               r = X(srd_particles(k),:) - cell_nodes(i,j)%pos
               if (abs(r(1)) > 10*DX) then
                  r(1) = r(1) - sign(LX,r(1)) !Correction for periodic boxes
               end if
               bw(k) = bilinear_weight(r)
            end do
            NN = sum(bw)

            if (NN == 0.D0) then
               print'(A,I0,A,I0,A)', "No surrounding particles at node (",i,',',j,")" 
               write(DPB_file_name,'(A,F0.2,A)') 'DPB_TIME=',TIME,'.dat'
               open(unit=42,file=DPB_file_name)
               do p = 1,size(DPB,2)
                  write(42,*) DPB(:,p)
               end do
               close(42)
               NN = 1.D0
            end if

            SCALAR_FIELD(i,j) = 0.D0
            do k = 1,size(srd_particles)
               SCALAR_FIELD(i,j) = SCALAR_FIELD(i,j) + LAGR_SCALAR(srd_particles(k)) * bw(k)
            end do
            SCALAR_FIELD(i,j) = SCALAR_FIELD(i,j)/NN
         end do
      end do
      !$OMP end parallel do

   end function CIC_SCALAR_AT_NODES

   function CIC_SCALAR_AT_PARTICLE(SCALAR_FIELD, k) result(LAGR_SCALAR)
      use SD_VARIABLES, only: X, SD_box
      real*8   , intent(in) :: SCALAR_FIELD(NXP,NZP)
      integer  , intent(in) :: k
      real*8                :: LAGR_SCALAR
      integer   :: i, j, m, n
      
      LAGR_SCALAR = 0.D0
      i = SD_box(k)%i
      j = SD_box(k)%j
      do m = i,i+1
         do n = j,j+1
            LAGR_SCALAR = LAGR_SCALAR + SCALAR_FIELD(m,n) * bilinear_weight(X(k,:) - cell_nodes(m,n)%pos)
         end do
      end do
   end function

   ! ----------- EXPERIMENTAL ---------------
   function CIC_SCALAR_AT_PARTICLE_2(SCALAR_FIELD, k) result(LAGR_SCALAR)
      use SD_VARIABLES, only: X
      real*8   , intent(in) :: SCALAR_FIELD(2,2)
      integer  , intent(in) :: k
      real*8                :: LAGR_SCALAR
      integer   :: m, n
      
      LAGR_SCALAR = 0.D0
      do m = 1,2
         do n = 1,2
            LAGR_SCALAR = LAGR_SCALAR + SCALAR_FIELD(m,n) * bilinear_weight(X(k,:) - cell_nodes(m,n)%pos)
         end do
      end do
   end function
end module CIC_ROUTINES

module SD_FUNCTIONS
    implicit none

    private :: FREEZE

contains

   subroutine stats(status)

      use SD_VARIABLES, only: R, R_mean, R_mean_ice, R_mean_water,sigma, N_sd, Frozen

      character(len = *), optional, intent(in) :: status
      integer   :: j, ice_count, water_count
      real*8    :: mean

      R_mean = (sum(R(2,:)**3)/N_sd)**(1.D0/3)
      mean  = sum(R(2,:))/N_sd
      sigma = sqrt(sum((R(2,:) - mean)**2)/N_sd)
      !Mean radius by phase
      R_mean_ice   = 0.D0
      ice_count    = 0
      R_mean_water = 0.D0
      water_count  = 0
      do j = 1,N_sd
         if (Frozen(j)) then
            R_mean_ice = R_mean_ice + R(2,j)**3
            ice_count = ice_count + 1;
         else
            R_mean_water = R_mean_water + R(2,j)**3
            water_count = water_count + 1;
         end if
      end do
      R_mean_ice   = (R_mean_ice/max(1,ice_count))**(1.D0/3)
      R_mean_water = (R_mean_water/max(1,water_count))**(1.D0/3)

      if (present(status) .and. status == 'open') then
         open(unit=1,file='stats.dat')
         write(1,*) '# R_mean   R_mean_ice   R_mean_water   sigma'
      end if
      write(1,*) R_mean, R_mean_ice, R_mean_water, sigma
      if (present(status) .and. status == 'close') then
         close(1)
         print*, 'Saved stats.dat'
      end if
   end subroutine stats

   subroutine WHICH_BOX(k)
      use GRID, only: DX,DZ
      use SD_VARIABLES, only: SD_box, X

      integer, intent(in) :: k
      ! OBS: Only works for rectangular and uniform grid.
      SD_box(k)%i = floor(X(k,1)/DX) + 1
      SD_box(k)%j = floor(X(k,2)/DZ) + 1
   end subroutine WHICH_BOX

   subroutine UPDATE_PARTICLE_BOX_MAP
      use SD_VARIABLES     , only: SD_box, DPB_input, N_sd, xi, X
      use THERMODYNAMIC_VAR, only: DPB, N_DROPS_BOX
      use GRID, only: NX, NZ, boxes
      

      integer   :: i, j, k
      integer  , allocatable :: temp_list(:)
      integer   :: prelocated_list_size
      integer, parameter :: prelocation_safety_factor = 2

      !Resetting particle lists
      prelocated_list_size = prelocation_safety_factor*DPB_input
      if (allocated(boxes)) deallocate(boxes)
      allocate(boxes(NX,NZ))
      do i=1,NX
         do j=1,NZ
            allocate(boxes(i,j)%p_list(prelocated_list_size))
         end do
      end do

      !Counting particles in each box and giving them an address.
      DPB = 0.D0;
      N_DROPS_BOX = 0.D0;
      do k = 1,N_sd
         call WHICH_BOX(k)
         i = SD_box(k)%i
         j = SD_box(k)%j
         !Debugging only
         if (i < 0 .or. i > NX .or. j < 0 .or. j > NZ) then
            print*, 'Fatal error: particle out of bounds'
            print*, 'Particle index: ', k
            print*, 'Particle position: ', X(k,:)
            print*, 'Box indeces:', i, j
            stop 
         end if
         !Debugging only ^^^
         DPB(i,j) = DPB(i,j) + 1
         N_DROPS_BOX(i,j) = N_DROPS_BOX(i,j) + xi(k)

         if (DPB(i,j) > prelocated_list_size) then
            !Append to the end of the list
            boxes(i,j)%p_list = [boxes(i,j)%p_list, k]
         else  
            !Assign to position 
            boxes(i,j)%p_list(DPB(i,j)) = k
         end if
      end do

      ! "Trimming" lists
      do i = 1,NX
         do j = 1,NZ
            allocate(temp_list(DPB(i,j)))
            temp_list = boxes(i,j)%p_list(1:DPB(i,j))
            call move_alloc(from=temp_list,to=boxes(i,j)%p_list)
         end do
      end do
   end subroutine UPDATE_PARTICLE_BOX_MAP

   function WET_RADIUS(dry) result(R_wet)
      use SD_VARIABLES     , only: SD_box, C_gamma, kappa, X
      use THERMODYNAMIC_VAR, only: TEMP, RV
      use CONSTANTS        , only: R_D, R_V
      use FUNCTIONS        , only: interpol, svp
      use ENVIRONMENT

      integer   :: k
      real*8    :: A, B, Rc, Sc
      real*8    :: error, tol, next
      real*8    :: F, F_prime, SAT, T, Q, e
      real*8, intent(in)  :: dry(:)
      real*8 :: R_wet(size(dry))


      A = C_gamma
      tol = 1.D-8
      !Initial super droplet radii - Newton-Rhapson Scheme
      do k = 1,size(R_wet)
         B = kappa*dry(k)**3
         error = 1.D0
         R_wet(k) = dry(k)   !Initial guess
         T   = TEMP(SD_box(k)%i,SD_box(k)%j)
         Q   = RV  (SD_box(k)%i,SD_box(k)%j)
         e   = Q/(Q + R_D/R_V)*interpol(Z_E,PR_E,X(k,2),linear_x=.true.)
         SAT = e/SVP(T,'L') - 1
         Rc  = sqrt(3*B/A)
         Sc  = A/Rc - B/(Rc**3)
         if (SAT > Sc) then
            R_wet(k) = Rc
         else
            do while (error > tol)
               F = SAT*R_wet(k)**3 - A*R_wet(k)**2 + B
               F_prime = 3*SAT*R_wet(k)**2 - 2*A*R_wet(k)
               next = R_wet(k) - F/F_prime
               error = abs(F)
               R_wet(k) = next
            end do
            if (isnan(R_wet(k))) then
               stop "WET_RADIUS: Failed to converge. Encountered value is NaN."
            end if
         end if
      end do
   end function WET_RADIUS

   subroutine INIT_SD_POSITIONS
      use SD_VARIABLES, only: DPB_input, X
      use GRID, only: NX, NZ, DX, DZ
      implicit none

      real*8    :: rn(2)
      integer   :: i, j, k,m

      call INIT_RANDOM_SEED
      k = 0
      do i = 1,NX
         do j = 1,NZ
            do m = 1,DPB_input
               k = k + 1
               call random_number(rn)
               X(k,1) = (i - rn(1))*DX
               X(k,2) = (j - rn(2))*DZ
            end do
         end do
      end do
   end subroutine INIT_SD_POSITIONS

   subroutine INIT_SD
      use CCN
      use SD_VARIABLES
      use THERMODYNAMIC_VAR, only: RV,THETA!, DPB
      use GRID             , only: NX, NZ
      use CUSTOM_TYPES
      implicit none
      integer   :: k

      N_sd = DPB_input*NX*NZ
      allocate(X(N_sd,2),XP(N_sd,2))
      allocate(SD_box(N_sd))
      allocate(R_dry(N_sd), R_crit(N_sd))
      allocate(R(2,N_sd))
      allocate(xi(N_sd), S_crit(N_sd))
      allocate(Activated(N_sd))
      allocate(phase_change(N_sd))
      allocate(Q_k(N_sd), TH_k(N_sd))
      allocate(u_prime(N_sd,2),u_prime_old(N_sd,2))
      allocate(Sk_log(N_sd))

      xi  = ceiling(n_drops/N_sd)*NX*NZ

      call INIT_SD_POSITIONS
      call UPDATE_PARTICLE_BOX_MAP
      call SAMPLE_DRY(R_dry)

      !Calculate initial wet radii
      do k = 1,N_sd
         Q_k(k)  = RV   (SD_box(k)%i , SD_box(k)%j)
         TH_k(k) = THETA(SD_box(k)%i , SD_box(k)%j)
      end do
      R(2,:) = WET_RADIUS(R_dry)
      R(1,:) = R(2,:)

      if (im_freezing) then
         call SAMPLE_TFREEZING(T_Freeze, N_sd)
      end if

      !Initially frozen droplets
      call INIT_FROZEN_ARRAY(Frozen, N_sd, init_ice_frac)

      !Initialize velocity fluctuations
      u_prime     = 0.D0
      u_prime_old = 0.D0
   end subroutine INIT_SD

   subroutine FREEZE(S,T,k)

      use SD_VARIABLES, only: R, R_dry, T_Freeze, Frozen, phase_change

      implicit none

      real*8   , intent(in) :: S, T
      integer  , intent(in) :: k
      logical :: A, B, C

      A = R(1,k) .GT. R_dry(k)
      B = S      .GT. 0.D0
      C = T      .LT. T_Freeze(k)

      if (Frozen(k)) then
         Frozen(k) = T < 273.15D0
         if (Frozen(k)) then
            phase_change(k) = 0
         else
            phase_change(k) = -1
         end if
      else
         Frozen(k) = A .AND. B .AND. C
         if (Frozen(k)) then
            phase_change(k) = 1
         else
            phase_change(k) = 0
         end if
      end if
   end subroutine FREEZE

   subroutine INIT_FROZEN_ARRAY(F_ARRAY, N, P)

      use FUNCTIONS, only: random_integer
      
      logical, allocatable, intent(out) :: F_ARRAY(:)
      integer             , intent(in ) :: N
      real*8              , intent(in ) :: P

      integer :: current_size, k

      current_size = floor( (1-P) * N )
      allocate(F_ARRAY(N))

      F_ARRAY(1:current_size) = .false.

      do while (current_size < N)
         call random_integer(1, current_size, k)

         F_ARRAY(current_size+1) = F_ARRAY(k)
         F_ARRAY(k) = .true.
         current_size = current_size + 1
      end do

   end subroutine


   subroutine GROW_DROPLETS
      use CONSTANTS
      use FUNCTIONS
      use THERMODYNAMIC_VAR
      use STEPPER           , only: dt
      use ENVIRONMENT       , only: Z_E, PR_E
      use ADVECTION         , only: stochastic_micro_phys, C_one, omega, z_eps
      use SD_VARIABLES      , only: SD_box, R, kappa, N_sd, phase_change, Frozen, X
      use SD_VARIABLES      , only: Q_k, TH_k, im_freezing, C_gamma, R_dry
      use OMP_LIB

      implicit none

      integer   :: k, m, n, count
      real*8, parameter  :: tol = 1.D-14
      real*8 :: error, next, F, F_prime, A
      real*8 :: D, EX, P, T, omega_local
      real*8 :: e_k, S_k, dq_k, dm_F, dTH_k ! Field localization parameters
      real*8 :: aw, Daw, Y, rd

      aw(y,rd)  = (y**(3.0/2) - rd**3)/(y**(3.0/2) - (rd**3)*(1 - kappa))
      Daw(y,rd) = 3.0/2*y**(1.0/2)/(y**(3.0/2) - (rd**3)*(1 - kappa))*(1 - aw(y,rd))

      !Solving the growth equation
      !$OMP parallel private(e_k,S_k,A,D,error,Y,count,F,F_prime,next,dq_k,dm_F,dTH_k,m,n,EX,P)
      !$OMP do
      do k = 1,N_sd
         !P = interpol(Z_E,PR_E,X(j,2)) !Pressure at the position of SD
         m  = SD_box(k)%i; n = SD_box(k)%j
         P  = interpol(Z_E,PR_E,X(k,2),linear_x=.true.)
         EX = EXNER(P)
         

         !Updating local mixing ratio and potential temperature
         if (stochastic_micro_phys) then
            dq_k     = -RHO_LIQ*N_DROPS_BOX(m,n)*4.0/3.0*pi*(R(2,k)**3 - R(1,k)**3)
            dm_F = phase_change(k)*N_DROPS_BOX(m,n)*RHO_LIQ*(4.D0/3.D0*pi)*R(1,k)**3
            if (Frozen(k)) then
               dTH_k = (-LS0*dq_k + (LS0 - LV0)*dm_F )/(CP_D * EX)
            else
               dTH_k = (-LV0*dq_k + (LS0 - LV0)*dm_F )/(CP_D * EX)
            end if

            
            omega_local = interpol(z_eps, omega, X(k,2), linear_x=.true.)

            Q_k(k)  = Q_k (k) + (   RV(m,n) -  Q_k(k)) * (C_one*omega_local*dt) + dq_k
            TH_k(k) = TH_k(k) + (THETA(m,n) - TH_k(k)) * (C_one*omega_local*dt) + dTH_k
         else
            !Conventional case: Homogeneous box with no fluctuations
            Q_k(k)  =    RV(m,n)
            TH_k(k) = THETA(m,n)
         end if

         !Saving old radius
         R(1,k) = R(2,k)

         e_k = P*Q_k(k)/(Q_k(k) + R_D/R_V)
         T   = TH_k(k)*EX

         if (im_freezing) then
            call FREEZE(S_k,T,k)
         end if

         if (Frozen(k)) then
            S_k = e_k/SVP(T,'I') - 1
            A = 0.D0
            D = (R_V*T*RHO_LIQ/D0/svp(T,'I') + LS0*RHO_LIQ/KT/T*(LS0/R_V/T - 1))**(-1)
         else
            S_k = e_k/SVP(TH_k(k)*EX,'L') - 1
            A = C_gamma
            D = (R_V*T*RHO_LIQ/D0/svp(T,'L') + LV0*RHO_LIQ/KT/T*(LV0/R_V/T - 1))**(-1)
         end if

         error = 1.D0
         Y = R(2,k)**2 !Initial guess for the Y = R² parameter
         count = 0
         do while (error > tol)
            F = 2*D*(S_k + 1 - aw(Y,R_dry(k))*exp(A/(Y**(1.D0/2))))*dt + R(2,k)**2 - Y
            F_prime = 2*D*exp(A/(Y**(1.D0/2)))*(aw(Y,R_dry(k))*A/2/(Y**(3.D0/2)) - Daw(Y,R_dry(k)))*dt - 1
            next = max(Y - F/F_prime,R_dry(k)**2)
            error = abs(F)
            Y = next
            count = count + 1
            if (count > 100) then
               stop 'GROW_DROPLETS: Growth equation did not converge'
            end if
         end do
         R(2,k) = sqrt(Y)
      end do
      !$OMP end do
      !$OMP end parallel
   end subroutine GROW_DROPLETS

   subroutine GET_ANALITICAL_K(KK,GRAD_KK,POS_Z)
      use ENVIRONMENT , only: EDDY_MAX, Z_INV, TT
      use CONSTANTS   , only: pi
      
      real*8, intent(in ) :: POS_Z
      real*8, intent(out) :: KK, GRAD_KK(2)

      KK         = 0.5D0*EDDY_MAX*ERFC((POS_Z - Z_INV)/(TT*sqrt(2.D0)))
      GRAD_KK(1) = 0.D0
      GRAD_KK(2) = -EDDY_MAX/(TT*sqrt(2.D0*pi))*exp(-(POS_Z - Z_INV)**2/(2*TT**2))
   end subroutine GET_ANALITICAL_K

   function B_COEFF(grad, k) result(b)
      use CIC_ROUTINES
      use ADVECTION   , only: C_one, omega, z_eps
      use FUNCTIONS   , only: interpol
      use SD_VARIABLES, only: X

      integer   :: i, j
      integer  , parameter  :: eye(2,2) = reshape([1,0,0,1],shape(eye))
      integer  , intent(in) :: k
      real*8   , intent(in) :: grad(NXP,NZP,2,2)
      real*8                :: b(2,2)
      real*8                :: omega_local
      
      ! Obtaining the negative transposed gradient
      do i = 1,2
         do j = 1,2
            b(i,j) = - CIC_SCALAR_AT_PARTICLE(grad(:,:,j,i), k)
         end do
      end do

      ! Adding the C1*w*I term
      omega_local = interpol(z_eps, omega, X(k,2), linear_x=.true.)
      b(:,:)      = b(:,:) - C_one * omega_local * eye(:,:)
   end function

   subroutine PARTICLE_BC(pos,L_x,L_z,vel)

      implicit none
      real*8          , intent(inout) :: pos(2), L_x, L_z
      real*8, optional, intent(inout) :: vel(2)
      real*8 ::  pos_correction(2), vel_correction(2)
      

      pos_correction = 0.D0; vel_correction = 0.D0
      ! Apply boundary conditions - Periodic(sides) / Reflection(top and bottom)
      if (pos(1) < 0.D0) then
         pos_correction(1) =    L_x
      elseif (pos(1) >  L_x) then
         pos_correction(1) =  - L_x
      end if

      if (pos(2) < 0.D0) then
         pos_correction(2) = - 2*pos(2)
         if (present(vel)) then
            vel_correction(2) = -2*vel(2)
         end if
      elseif (pos(2) >  L_z) then
         pos_correction(2) = - 2*(pos(2) - L_z)
         if (present(vel)) then
            vel_correction(2) = -2*vel(2)
         end if
      end if

      pos(:) = pos(:) + pos_correction(:)
      if (present(vel)) then
         vel(:) = vel(:) + vel_correction(:)
      end if
   end subroutine

   subroutine UPDATE_U_PRIME
      use CIC_ROUTINES
      use FUNCTIONS
      use ADVECTION, only: UXN, UZN, C_zero, z_eps, eps, stochastic_micro_phys
      use STEPPER  , only: DT_ADV
      use SD_VARIABLES, only: u_prime, u_prime_old, X, XP, N_sd
      use omp_lib

      real*8    :: uu(NXP,NZP), uw(NXP,NZP), ww(NXP,NZP)
      real*8    :: div_uu_x(NXP,NZP), div_uu_z(NXP,NZP)
      real*8    :: grad(NXP,NZP,2,2)
      real*8    :: a(2), b(2,2), c, ksi(2)
      real*8    :: eps_local, u_local(2), Du_prime(2), u_avg(2)
      integer   :: i, j, k


      if (.not.stochastic_micro_phys) then
         u_prime_old = u_prime
         u_prime = 0.D0
         XP = X
         return
      end if

      !Calculate the components of <u'u'> on the nodes
      uu = CIC_SCALAR_AT_NODES(u_prime(:,1)*u_prime(:,1))
      uw = CIC_SCALAR_AT_NODES(u_prime(:,1)*u_prime(:,2))
      ww = CIC_SCALAR_AT_NODES(u_prime(:,2)*u_prime(:,2))

      !Calculate the divergence of <u'u'>
      call DIV_SYM_TENSOR_2D(DX,DZ,uu,uw,ww,div_uu_x,div_uu_z)

      !Calculate the gradient of <u> at the nodes
      call GRAD_VECTOR_2D(DX,DZ,UXN,UZN,grad)

      !!$OMP parallel do private(u_local, a, b, c, eps_local, ksi, Du_prime)
      do k = 1,N_sd
         !Obtain mean flow velocity at particle position
         u_local(1) = CIC_SCALAR_AT_PARTICLE(UXN,k)
         u_local(2) = CIC_SCALAR_AT_PARTICLE(UZN,k)
      
         !Advance particle to midpoint position
         XP(k,:) = X(k,:) !Save previous position
         X (k,:) = XP(k,:) + 0.5*(u_local + u_prime(k,:))*DT_ADV
         
         !Apply boundary conditions - Periodic(sides) / Reflection(top and bottom)
         call PARTICLE_BC(X(k,:), LX, LZ)
         
         !Update box addresses (some particles may have changed boxes)
         call WHICH_BOX(k)

         !Obtain coefficient a (div(<u'u'>) at particle position)
         a(1) = CIC_SCALAR_AT_PARTICLE(div_uu_x, k)
         a(2) = CIC_SCALAR_AT_PARTICLE(div_uu_z, k)

         !Obtain b = -(grad <u>)' + C_one*omega*eye(2) at particle position
         b = B_COEFF(grad, k)
         
         !Obtain c = C_zero*eps at particle position
         eps_local = interpol(z_eps,eps,X(k,2), linear_x=.true.)
         c = C_zero*eps_local
         
         !Get random number with standard gaussian distribution
         call SAMPLE_GAUSSIAN(ksi,1.D0,0.D0)
         
         !Runge-Kutta of 2nd order to advance u_prime in time
         Du_prime = (a + matmul(b,u_prime_old(k,:)))*DT_ADV + sqrt(c*DT_ADV)*ksi
         u_prime_old(k,:) = u_prime(k,:)
         u_prime(k,:) = u_prime_old(k,:) + Du_prime + 0.5*matmul(b,Du_prime)*DT_ADV
      end do
      !!$OMP end parallel do
    
      !Correction for zero average in each cell
      call UPDATE_PARTICLE_BOX_MAP

      do i = 1,NX
         do j = 1,NZ
            u_avg(1) = sum(u_prime(boxes(i,j)%p_list,1)) / size(boxes(i,j)%p_list)
            u_avg(2) = sum(u_prime(boxes(i,j)%p_list,2)) / size(boxes(i,j)%p_list)
            u_prime(boxes(i,j)%p_list,1) = u_prime(boxes(i,j)%p_list,1) - u_avg(1)
            u_prime(boxes(i,j)%p_list,2) = u_prime(boxes(i,j)%p_list,2) - u_avg(2)
         end do
      end do
      
   end subroutine UPDATE_U_PRIME

   SUBROUTINE POSITION_CORRECTION (POS)
  
      !The positions are corrected to force consistency
     
      use SD_VARIABLES, only: N_sd
      USE GRID        , only: NX, NZ, DX, DZ, LX, LZ
      use OMP_LIB
     
      IMPLICIT NONE
     
      REAL*8, intent(inout) :: POS(N_sd,2)
      REAL*8 :: DPX(NX+1,NZ),DPZ(NX,NZ+1) !Gradient at the cell faces
      REAL*8 :: DPX_K,DPZ_K
    
      REAL*8 :: X0,Z0,C,D !Variables for linear interpolation
     
      INTEGER :: IX,IZ, k

      !Poisson problem for mean particle kinematic pressure (position correction)
      CALL EVAL_KINEMATIC_PARTICLE_PRESSURE_POS (PX=DPX, PZ=DPZ)
      
      !$OMP parallel do private(IX, IZ, X0, Z0, C, D, DPX_K, DPZ_K)
      do k = 1,N_sd
         !Linear interpolation of the pressure gradient at the particle
         !position - the pressure gradient is stored at the cell faces
         IX = INT(POS(k,1)/DX) + 1
         IZ = INT(POS(k,2)/DZ) + 1
         
         !NOTA: Alterar para usar as posições dos nós (já calculadas) - vale a pena?
         X0 = DX*(IX-1)
         Z0 = DZ*(IZ-1)
   
         C = (POS(k,1)-X0)/DX
         D = (POS(k,2)-Z0)/DZ
         
         DPX_K = (1.D0-C)*DPX(IX,IZ) + C*DPX(IX+1,IZ  )
         DPZ_K = (1.D0-D)*DPZ(IX,IZ) + D*DPZ(IX  ,IZ+1)
         
         !Correction
         POS(k,:) = POS(k,:) - [DPX_K, DPZ_K]
         call PARTICLE_BC(POS(k,:), LX, LZ)    
      end do
      !$OMP end parallel do

      call UPDATE_PARTICLE_BOX_MAP
   END SUBROUTINE POSITION_CORRECTION

   SUBROUTINE EVAL_KINEMATIC_PARTICLE_PRESSURE_POS (P,PX,PZ)

      USE GRID             , only: NX, NZ, NXP, NZP, DX, DZ
      use SD_VARIABLES     , only: N_sd
      use THERMODYNAMIC_VAR, only: DPB
      use NUMERICAL_METHODS, only: poisson_2d_mixed_bc_xper_yneu_stagg
   
      IMPLICIT NONE
   
      REAL*8, optional, intent(out) :: P(NX,NZ)
      real*8, intent(out) :: PX(NXP,NZ),PZ(NX,NZP)
   
      REAL*8 :: DENS_0
   
      REAL*8 :: S(0:NX-1,1:NZ),F(0:NX-1,1:NZ)
      REAL*8 :: GB(0:NX-1),GT(0:NX-1)
   
      INTEGER :: i , j
      INTEGER :: IX, IZ
   
      !real*8  :: lapl,s1,s2
   
      !*************************************************************
      !Right-hand side (RHS) of the Poisson equation for "P"
      !-------------------------------------------------------------
   
      DENS_0 = REAL(N_sd)/(NX*NZ)
   
      do i = 1,NX
         do j = 1,NZ
            F(i-1,j) = 1.D0 - DPB(i,j)/DENS_0 !Right-hand side (source)
         end do
      end do
   
      !Bottom and top Neumann boundary conditions 
      DO IX = 1,NX
         GB(IX-1) = 0.D0 
         GT(IX-1) = 0.D0 
      END DO
   
      !********************************************************************
      !Poisson solver
      call poisson_2d_mixed_bc_xper_yneu_stagg (NX,NZ,DX,DZ,F,GB,GT,S)

      !********************************************************************
      !Pressure gradient in the staggered grid
   
      !1) Vertical gradient
      DO IX = 1,NX
         !Interior points
         DO IZ = 2,NZ
            PZ(IX,IZ) = ( S(IX-1,IZ) - S(IX-1,IZ-1) ) / DZ
         END DO
         
         !Boundary points
         PZ(IX,  1) = GB(IX-1) !According to the Neumann b. cond. (bottom)
         PZ(IX,NZP) = GT(IX-1) !According to the Neumann b. cond. (top)
      END DO
      
      !2) Horizontal gradient
      DO IZ = 1,NZ
         !Interior points
         DO IX = 2,NX
            PX(IX,IZ) = ( S(IX-1,IZ) - S(IX-2,IZ) ) / DX
         END DO
         
         !Boundary points - periodicity
         PX(  1,IZ) = ( S(0,IZ) - S(NX-1,IZ) ) / DX
         PX(NXP,IZ) = PX(1,IZ)
      END DO
      
      if (present(P)) then
         DO IX = 1,NX
            DO IZ = 1,NZ
               P(IX,IZ) = S(IX-1,IZ)
            END DO
         END DO
      end if

   END SUBROUTINE EVAL_KINEMATIC_PARTICLE_PRESSURE_POS

   subroutine ADVECTION_SD
      use STEPPER     , only: DT_ADV
      use GRID        , only: LX , LZ!, cell_nodes
      use ADVECTION   , only: UXN, UZN, ADVECT
      use CIC_ROUTINES, only: CIC_SCALAR_AT_PARTICLE, CIC_SCALAR_AT_PARTICLE_2
      use SD_VARIABLES, only: N_sd, u_prime, u_prime_old, X, XP
      use OMP_LIB

      integer   :: k
      real*8    :: u_local(2)

      if (.NOT.ADVECT) return
      
      !Update velocity fluctuations and move particles to midpoint
      call UPDATE_U_PRIME 
      
      !!$OMP parallel do private(u_local, i,j,UX_local,UZ_local) firstprivate(p_X)
      do k = 1,N_sd
         !Obtain local velocity at midpoint (u_prime is temporarily overriden to zero)
         u_local(1) = CIC_SCALAR_AT_PARTICLE(UXN,k) + 0.5*(u_prime_old(k,1) + u_prime(k,1))
         u_local(2) = CIC_SCALAR_AT_PARTICLE(UZN,k) + 0.5*(u_prime_old(k,2) + u_prime(k,2))

         !Calculate the provisional position (without compressibility correction)
         X(k,:) = XP(k,:) + u_local * DT_ADV
      
         ! Apply boundary conditions - Periodic(sides) / Reflection(top and bottom)
         call PARTICLE_BC(X(k,:), LX, LZ, u_prime(k,:))         
      end do
      !!$OMP end parallel do
      
      !Update box-particle maps before applying position correction based on particle density
      call UPDATE_PARTICLE_BOX_MAP
      
      !Correct droplet positions in order to mantain constant droplet density(second map update inside)
      CALL POSITION_CORRECTION(X)
      
   end subroutine ADVECTION_SD

end module SD_FUNCTIONS

module BOX
   implicit none

contains

   subroutine UPDATE_BOXES
      use CONSTANTS
      use FUNCTIONS        , only: EXNER, interpol
      use ENVIRONMENT      , only: Exn
      use GRID             , only: NX, NZ
      use SD_VARIABLES     , only: R, xi, N_sd, Frozen, phase_change, R_dry, SD_box
      use SD_VARIABLES     , only: Q_k, TH_k, u_prime
      use THERMODYNAMIC_VAR, only: RL, RI, RV, THETA, TEMP
      use ADVECTION        , only: ADVECT
      use OMP_LIB


      integer   :: k, i, j
      real*8    :: DELTA_L(NX,NZ), DELTA_I(NX,NZ), DELTA_F(NX,NZ)
      real*8, allocatable, save    :: FTH_TURB_DT(:,:), FRV_TURB_DT(:,:)

     
      if (.NOT.allocated(FTH_TURB_DT)) allocate(FTH_TURB_DT(NX,NZ), FRV_TURB_DT(NX,NZ))

      DELTA_L = 0.D0
      DELTA_I = 0.D0
      DELTA_F = 0.D0

      RL      = 0.D0
      RI      = 0.D0

      !$OMP parallel do private(i,j) reduction(+:DELTA_I,DELTA_L,DELTA_F,RI,RL)
      do k = 1,N_sd
         i = SD_box(k)%i; j = SD_box(k)%j
         if (Frozen(k)) then
            DELTA_I(i,j) = DELTA_I(i,j) + xi(k)*(R(2,k)**3 - R(1,k)  **3)
            RI(i,j)      = RI(i,j)      + xi(k)*(R(2,k)**3 - R_dry(k)**3)
         else
            DELTA_L(i,j) = DELTA_L(i,j) + xi(k)*(R(2,k)**3 - R(1,k)  **3)
            RL(i,j)      = RL(i,j)      + xi(k)*(R(2,k)**3 - R_dry(k)**3)
         end if
         DELTA_F(i,j)     = DELTA_F(i,j) + xi(k)*(R(1,k)**3 - R_dry(k)**3)*phase_change(k)
      end do
      !$OMP end parallel do

      DELTA_I(:,:) = DELTA_I(:,:) * RHO_LIQ*(4.D0/3*pi)
      DELTA_L(:,:) = DELTA_L(:,:) * RHO_LIQ*(4.D0/3*pi)
      DELTA_F(:,:) = DELTA_F(:,:) * RHO_LIQ*(4.D0/3*pi)
      RI(:,:)      = RI(:,:)      * RHO_LIQ*(4.D0/3*pi)
      RL(:,:)      = RL(:,:)      * RHO_LIQ*(4.D0/3*pi)

      
      !Obtain turbulent fluxes
      if (ADVECT) then
         call GET_TURB_FLUXES_CIC(u_prime, Q_k, FRV_TURB_DT)
         call GET_TURB_FLUXES_CIC(u_prime, TH_k, FTH_TURB_DT) 
      end if

      !Updating potential temperature and vapor mixing ratio
      do j = 1,NZ
         THETA(:,j) = THETA(:,j) + (LV0*DELTA_L(:,j) + LS0*DELTA_I(:,j) + DELTA_F(:,j)*(LS0 - LV0)) / (CP_D*Exn(j)) &
                    + FTH_TURB_DT(:,j)
         TEMP (:,j)  = Exn(j)*THETA(:,j)
      end do
      RV(:,:) = RV(:,:) - DELTA_L(:,:) - DELTA_I(:,:) + FRV_TURB_DT(:,:)

   end subroutine UPDATE_BOXES

   subroutine SAT_FIELD
      use ENVIRONMENT      , only: PR_E
      use THERMODYNAMIC_VAR, only: TEMP, RV, SAT
      use CONSTANTS        , only: R_D, R_V
      use FUNCTIONS        , only: SVP, interpol
      use GRID             , only: NX, NZ
      real*8     :: e(NX,NZ),P(NX,NZ)
      integer    :: j

      do j = 1,NZ
         P(:,j) = PR_E(j)
      end do

      e   = P*RV/(RV + R_D/R_V)
      SAT = e/SVP(TEMP,'L') - 1
   end subroutine SAT_FIELD

   subroutine GET_SCALAR_FORCES_TURB(FTH_TURB_DT, FRV_TURB_DT)

      !Get the "forces" (i.e., the RHS of THE Eulerian (scalar) advection
      !equations) accounting for turbulent fluxes

      USE ENVIRONMENT, only: EDDY_DIFF
      USE THERMODYNAMIC_VAR
      USE STEPPER
      USE GRID

      real*8, allocatable ::  JX(:,:), JZ(:,:), DIV(:,:)
      real*8, intent(out) :: FTH_TURB_DT(NX,NZ), FRV_TURB_DT(NX,NZ)
      integer   :: i

      allocate ( JX(NX,NZ), JZ(NX,NZ), DIV(NX,NZ) )

      !Turbulent flux of potential temperature THETA(:,:)
      call EVAL_GRADIENT_SCALAR ( THETA, JX, JZ )

      do i = 1,NX
         JX(i,:) = EDDY_DIFF(:)*JX(i,:)
         JZ(i,:) = EDDY_DIFF(:)*JZ(i,:)
      end do

      JZ(1:NX, 1) = 0.D0  !ZERO vertical turbulent flux boundary condition
      JZ(1:NX,NZ) = 0.D0  !ZERO vertical turbulent flux boundary condition

      call EVAL_DIVERGENCE ( DIV, JX, JZ )

      FTH_TURB_DT = DIV*DT

      !Turbulent flux of vapor mixing ratio QV(:,:)
      call EVAL_GRADIENT_SCALAR ( RV, JX, JZ )

      do i = 1,NX
         JX(i,:) = EDDY_DIFF(:)*JX(i,:)
         JZ(i,:) = EDDY_DIFF(:)*JZ(i,:)
      end do

      JZ(1:NX, 1) = 0.D0  !ZERO vertical turbulent flux boundary condition
      JZ(1:NX,NZ) = 0.D0  !ZERO vertical turbulent flux boundary condition

      call EVAL_DIVERGENCE ( DIV, JX, JZ )

      FRV_TURB_DT = DIV*DT

      deallocate (JX,JZ,DIV)
   end subroutine GET_SCALAR_FORCES_TURB

   subroutine GET_TURB_FLUXES_CIC(u_prim, LAGR_SCALAR, F_TURB)
      !Get the "forces" (i.e., the RHS of THE Eulerian (scalar) advection
      !equations) accounting for turbulent fluxes
      use GRID        , only: NXP, NZP, NX, NZ, DX, DZ
      use SD_VARIABLES, only: N_sd
      use FUNCTIONS   , only: DIV_VECTOR_2D
      use CIC_ROUTINES, only: CIC_SCALAR_AT_NODES, CIC_SCALAR_AT_PARTICLE
      use STEPPER     , only: DT
      use OMP_LIB

      implicit none


      real*8, intent(in ) :: u_prim(N_sd,2), LAGR_SCALAR(N_sd)
      real*8, intent(out) :: F_TURB(NX,NZ)
      real*8 :: EUL_AVG_GRID(NXP,NZP), COV_GRID(NXP,NZP,2)
      real*8 :: LAG_AVG_PARTICLE(N_sd), COV_PARTICLE(N_sd,2)
      integer :: i, j, k

      !Obtain a grid with average values
      EUL_AVG_GRID = CIC_SCALAR_AT_NODES(LAGR_SCALAR)
      
      ! OBS: Parallelizing this section makes it worse for some reason
      !!$OMP parallel do private(LAG_AVG_PARTICLE)
      do k = 1,N_sd
         !Get the average values at the particle position
         LAG_AVG_PARTICLE(k)   = CIC_SCALAR_AT_PARTICLE(EUL_AVG_GRID,k)
         !calculate u'phi'
         COV_PARTICLE (k,:) = (LAGR_SCALAR(k) - LAG_AVG_PARTICLE(k)) * u_prim(k,:)
      end do
      !!$OMP end parallel do
      
      !Obtain a grid with <u'phi'>
      COV_GRID(:,:,1) = CIC_SCALAR_AT_NODES(COV_PARTICLE(:,1)) ! X component
      COV_GRID(:,:,2) = CIC_SCALAR_AT_NODES(COV_PARTICLE(:,2)) ! Z component
      
      !Calculate div(<u'phi'>)
      call DIV_VECTOR_2D(DX,DZ,COV_GRID,EUL_AVG_GRID)

      !Back to primary grid
      do i = 1,NX
         do j = 1,NZ
            F_TURB(i,j) = sum(EUL_AVG_GRID(i:i+1,j:j+1)) * DT/4.D0
         end do
      end do
   end subroutine GET_TURB_FLUXES_CIC

   function PER(I,N)

      integer   :: PER, I, N

      if (I.LT.1) then
         PER = I + N
      elseif (I.GT.N) then
         PER = I - N
      else
         PER = I
      end if
   end function PER

   subroutine EVAL_GRADIENT_SCALAR (S,DSX,DSZ)

      !Evaluates the gradient of the scalar "S" at every point of the
      !primary grid. Periodic boundary conditions in the horizontal (x)
      !direction and no-flux BC in the vertical (z) direction.

      USE GRID

      real*8    :: S(NX,NZ), DSX(NX,NZ), DSZ(NX,NZ)
      integer   :: IX, IZ

      do IX = 1,NX
         do IZ = 1,NZ
            !Horizontal gradient (periodic boundary conditions)
            DSX(IX,IZ) = 0.5D0/DX*( S(PER(IX+1,NX),IZ) - S(PER(IX-1,NX),IZ) )

            !Vertical gradient
            if (IZ.EQ.1) then
               DSZ(IX,IZ) = 0.5D0/DZ*( -3.D0*S(IX,IZ) + 4.D0*S(IX,IZ+1) - S(IX,IZ+2) )
            else if (IZ.EQ.NZ) then
               DSZ(IX,IZ) = 0.5D0/DZ*(  3.D0*S(IX,IZ) - 4.D0*S(IX,IZ-1) + S(IX,IZ-2) )
            else
               !Interior points
               DSZ(IX,IZ) = 0.5D0/DZ*( S(IX,IZ+1) - S(IX,IZ-1) )
            end if
         end do
      end do
   end subroutine EVAL_GRADIENT_SCALAR

   subroutine EVAL_DIVERGENCE (DIV,JX,JZ)

      !Evaluates the divergence of the vector (JX,JZ) at every point of
      !the primary grid. Periodic boundary conditions in the horizontal
      !(x) direction and no-flux in the vertical (z) direction.

      USE GRID

      real*8    :: DIV(NX,NZ), JX(NX,NZ), JZ(NX,NZ)
      real*8    :: AX, AZ
      integer   :: IX, IZ

      do IX = 1,NX
         do IZ = 1,NZ
            !Horizontal contribution (periodic boundary conditions)
            AX = 0.5D0/DX*( JX(PER(IX+1,NX),IZ) - JX(PER(IX-1,NX),IZ) )

            !Vertical contribution
            if (IZ.EQ.1) then
               AZ = 0.5D0/DZ*( -3.D0*JZ(IX,IZ) + 4.D0*JZ(IX,IZ+1) - JZ(IX,IZ+2) )
            else if (IZ.EQ.NZ) then
               AZ = 0.5D0/DZ*(  3.D0*JZ(IX,IZ) - 4.D0*JZ(IX,IZ-1) + JZ(IX,IZ-2) )
            else
               !Interior points
               AZ = 0.5D0/DZ*( JZ(IX,IZ+1) - JZ(IX,IZ-1) )
            end if

            DIV(IX,IZ) = AX + AZ
         end do
      end do
   end subroutine EVAL_DIVERGENCE

end module BOX

module ADV_FUNCTIONS

   implicit none

contains

   subroutine CHECK_ADV
      use ADVECTION, only: ADVECT
      use STEPPER  , only: TIME, TIME_ADV, DT_ADV, DT
      integer  , save :: IT = 0, IT_OLD = 0

      ! Workarround for float precision
      if (abs(TIME/DT_ADV - nint(TIME/DT_ADV)) < DT/2) then
         IT = nint(TIME/DT_ADV)
      else
         IT = floor(TIME/DT_ADV)
      end if

      if (IT == IT_OLD) then
         ADVECT = .FALSE.
      else
         IT_OLD   = IT
         ADVECT   = .TRUE.
         TIME_ADV = [TIME, TIME + DT_ADV]
      end if
   end subroutine CHECK_ADV

   function COURANT()
      use ADVECTION, only: UXT, UZT
      use STEPPER  , only: DT_ADV
      use GRID     , only: DX, DZ
      real*8 :: CX, CZ, COURANT
      real*8, save :: C_MAX = 0.D0

      CX = maxval(UXT)*DT_ADV/DX
      CZ = maxval(UZT)*DT_ADV/DZ
      COURANT = max(CX,CZ)

      if (COURANT > C_MAX) then
         C_MAX = COURANT
      else
         COURANT = C_MAX
      end if
   end function COURANT

   subroutine GET_VELOCITIES

      !Gets the velocities "UXT" and "UZT" in the dual grid and the flow
      !velocities "UX" and "UZ" (in the PRIMARY grid) by interpolation
      !Derived numerically from the prescribed INCOMPRESSIBLE streamfunction
      use ADVECTION
      use GRID
      use IO_PARAMETERS, only: input_file
      use HDF5_functions, only: HDF5_READ_PSI

      integer         :: I,K
      logical, save   :: first_call = .true.

      if (first_call) then
         allocate (UX  (NX ,NZ ), UZ  (NX ,NZ ))
         allocate (UXN (NXP,NZP), UZN (NXP,NZP))
         allocate (UXT (NXP,NZ ), UZT (NX ,NZP))
         allocate (UXTP(NXP,NZ ), UZTP(NX ,NZP))
         allocate (PSI (NXP,NZP), PSIP(NXP,NZP))
         UXT  = 0.D0; UZT  = 0.D0

         call HDF5_READ_PSI(input_file)
         call GET_VEL_FROM_PSI_STAGGERED_GRID
      end if
      
      ! Get transport velocities UXT and UZT, as well as
      ! dual grid velocities UXN and UZN
      if (ADVECT.or.first_call) then
         ! Save Previous values of velocity and stream function
         UXTP = UXT
         UZTP = UZT
         PSIP = PSI

         call HDF5_READ_PSI(input_file)
         call GET_VEL_FROM_PSI_STAGGERED_GRID
         !Flow velocity "UX" and "UZ" in the primary grid by interpolation
         do I = 1,NX
            do K = 1,NZ
               UX(I,K) = 0.5D0*( UXTP(I,K) + UXTP(I+1,K  ) )
               UZ(I,K) = 0.5D0*( UZTP(I,K) + UZTP(I  ,K+1) )
            end do
         end do
         call GET_ADVECTIVE_VELOCITIES
         if (first_call) first_call = .false.
      end if
      
      call GET_VEL_FROM_PSI_DUAL_GRID

   end subroutine GET_VELOCITIES

   subroutine GET_VEL_FROM_PSI_STAGGERED_GRID

      use ADVECTION, only: UXT, UZT, PSI
      use GRID     , only: NX, NZ, NXP, NZP, DX, DZ
      use OMP_LIB

      integer       :: i, k

      !Spatial derivation of the streamfunction

      ! Obtaining transport velocities at the faces of cells
      !$OMP parallel private(k)
      !UX velocity
      !$OMP do
      do i = 1,NXP
         do k = 1,NZ
            UXT(I,K) = - ( PSI(I,K+1) - PSI(I,K) ) / DZ
         end do
      end do
      !$OMP end do

      !UZ velocity
      !$OMP do
      do i = 1,NX
         do k = 1,NZP
            UZT(I,K) = ( PSI(I+1,K) - PSI(I,K) ) / DX
         end do
      end do
      !$OMP end do
      !$OMP end parallel
   end subroutine GET_VEL_FROM_PSI_STAGGERED_GRID

   subroutine GET_VEL_FROM_PSI_DUAL_GRID

      use GRID     , only: NXP, NZP, DX, DZ, NZ
      use STEPPER  , only: TIME_ADV, TIME
      use ADVECTION, only: PSIP, PSI, UXN, UZN
      use FUNCTIONS, only: interpol, PERN
      use OMP_LIB

      integer   :: i, k
      real*8    :: PSI_now(NXP,NZP)
      
      ! Interpolate the stream function to the present moment
      !$OMP parallel private(k)
      !$OMP do
      do i = 1,NXP
         do k = 1,NZP
            PSI_now(i,k) = interpol(TIME_ADV,[PSIP(i,k), PSI(i,k)],TIME, linear_x=.true.)
         end do
      end do
      !$OMP end do
      !Obtaining mean velocity at nodes (2nd order expressions)
      !$OMP do
      do i = 1,NXP
         UXN(i,1) = ( 3*PSI_now(i,1) - 4*PSI_now(i,2) + PSI_now(i,3))/(2*DZ)
         UZN(i,1) = (PSI_now(PERN(i+1,NXP),1) - PSI_now(PERN(i-1,NXP),1))/(2*DX)
         do k = 2,NZ
            UXN(i,k) = (-PSI_now(i,k+1) + PSI_now(i,k-1))/(2*DZ)
            UZN(i,k) = (PSI_now(PERN(i+1,NXP),k) - PSI_now(PERN(i-1,NXP),k))/(2*DX)
         end do
         UXN(i,NZP) = (-3*PSI_now(i,NZP) + 4*PSI_now(i,NZP-1) - PSI_now(i,NZP-2))/(2*DZ)
         UZN(i,NZP) = (PSI_now(PERN(i+1,NXP),NZP) - PSI_now(PERN(i-1,NXP),NZP))/(2*DX)

      end do
      !$OMP end do
      !$OMP end parallel
   end subroutine GET_VEL_FROM_PSI_DUAL_GRID

   subroutine GET_ADVECTIVE_VELOCITIES

      !Gets the velocities "UXA" and "UZA" in the dual grid (at time level
      !"n+1/2") for the advection by MPDATA.

      !*** Uses "UXT" and "UZT" and their values "UXTP" and "UZTP" in the previous time level ***

      use ADVECTION
      use GRID
      use STEPPER, only: DT_ADV
      use ENVIRONMENT, only: RHO

      real*8    :: RHO_BAR
      integer   :: IZ
      logical, save :: first_call = .true.

      if (first_call) then
         allocate (UXA(NXP,NZ ))
         allocate (UZA(NX ,NZP))
         first_call = .false.
      end if

      !Horizontal velocities (applies periodic boundary conditions)

      !Interpolating horizontal velocities
      UXA = 0.5D0*(UXT + UXTP)*DT_ADV/DX

      !Interpolating vertical velocities
      UZA = 0.5D0*(UZT + UZTP)*DT_ADV/DZ

      !Multiplying by the base-state dry density RHO
      do IZ = 1,NZ
            UXA(:,IZ) = RHO(IZ)*UXA(:,IZ)
      end do

      RHO_BAR = 0.5D0*( 3.D0*RHO(1) - RHO(2) )
      UZA(:,1) = RHO_BAR*UZA(:,1)
      do IZ = 2,NZ
         RHO_BAR = 0.5D0*( RHO(IZ-1) + RHO(IZ) )
         UZA(:,IZ) = RHO_BAR*UZA(:,IZ)
      end do
      RHO_BAR = 0.5D0*( 3.D0*RHO(NZ) - RHO(NZ-1) )
      UZA(:,NZP) = RHO_BAR*UZA(:,NZP)
   end subroutine GET_ADVECTIVE_VELOCITIES

   subroutine ADVECTION_MPDATA

      !2-D advection of "theta" and "RV" fields
      !To be used in the Super-droplet scheme

      use THERMODYNAMIC_VAR
      use ADVECTION
      use GRID

      integer :: IORD,ISOR,NONOS,IDIV !MPDATA parameters

      if (.NOT.ADVECT) return
      
      IORD  = 2
      ISOR  = 0
      NONOS = 1
      IDIV  = 0

      call MPDATA2D (UXA,UZA,THETA,GAC,NX,NZ,IORD,ISOR,NONOS,IDIV,1)
      call MPDATA2D (UXA,UZA,RV   ,GAC,NX,NZ,IORD,ISOR,NONOS,IDIV,2)

   end subroutine ADVECTION_MPDATA

   function INTERMITTENCY_FACTOR(z, z_inv, DZ) result(gamma)

      real*8, intent(in) :: z, z_inv, DZ
      real*8             :: gamma

      gamma = 0.5 - 0.5 * erf( (z - z_inv) / (DZ) )

   end function

   subroutine INIT_SG_PARAMETERS
      use ADVECTION    , only: z_eps, eps, omega, L
      use GRID         , only: NZP
      use IO_PARAMETERS, only: eps_file_name
      use STEPPER      , only: DT_ADV

      integer           :: i
      integer           :: eof               !End-of-file indicator
      real*8            :: z_read, eps_read  !Values read form file
      real*8, parameter :: C = 1.D0          !Proportionality constant in k ~ (eps*L)^(2/3)
      real*8            :: omega_laminar, gamma, eps_transition

      !Read eps from data file
      allocate(z_eps(0), eps(0))
      open(unit=1,file=eps_file_name, status='old')
      read(unit=1,fmt=*,iostat=eof) z_read, eps_read
      do while (eof==0)
         z_eps = [z_eps, z_read  ]
         eps   = [eps  , eps_read]
         read(unit=1,fmt=*,iostat=eof) z_read, eps_read
      end do
      close(1)

      !Size check
      if (size(eps) < NZP) then
         call system('clear')
         print*, "WARNING: Resolution for TKE profile is too low! Check input file"
         stop
      end if

      !Calculate omega as a function of eps (lookup table)
      allocate(omega(size(eps)))
      omega_laminar = 1 / DT_ADV
      eps_transition = 0.05 * maxval(eps)
      do i = 1,size(omega)
         gamma = INTERMITTENCY_FACTOR(eps(i), eps_transition, -eps_transition/2)
         omega(i) = (1/C)*(eps(i)/L**2)**(1.D0/3) * gamma + (1 - gamma) * omega_laminar
      end do

   end subroutine INIT_SG_PARAMETERS
end module ADV_FUNCTIONS