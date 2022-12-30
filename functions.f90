MODULE functions
    use CONSTANTS
    use OUTPUT, only: FILE_NAME
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
!==========================DEFINIÇÃO DE FUNÇÕES=================================

    subroutine progress_bar(frac)
        real*8, intent(in) :: frac
        integer :: i, bar_size,sticks
        character, allocatable :: bar(:)
        real*8   , save :: frac_old = 0.0
        integer*8, save :: start = 0, old = 0, rate = 0
        integer*8       :: now
        real*8          :: rem_time, total_time, elapsed_time

        if (start == 0) then
            call system_clock(count_rate=rate)
            call system_clock(start)
            old = start
            return
        end if

        if (frac > 1.00001D0 .or. frac < 0.D0) then
            print*,'Fraction should be between 0 and 1'
            return
        end if

        bar_size = 50
        allocate(bar(bar_size+2))
        bar(1) = '['
        sticks = nint(frac*bar_size)
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
        print '(A,F6.2,A)', 'Progresso: ', frac*100,'%'
        print*, bar
        print'(A,A)','Tempo Decorrido     :  ', PRETTY_TIME(elapsed_time)
        print'(A,A)','Tempo Restante      :  ', PRETTY_TIME(rem_time    )
        print'(A,A)','Tempo Total Estimado:  ', PRETTY_TIME(total_time  )
        print*, 'Output file: ', trim(FILE_NAME)
    end subroutine progress_bar

    function PRETTY_TIME(time_input) result(time_string)

        real*8, intent(in)    :: time_input
        integer*8    :: seconds, minutes, hours
        character(9) :: time_string

        
        hours   = nint(time_input)/3600
        minutes = nint(time_input)/60 - hours*60
        seconds = mod(nint(time_input),60)

        write(time_string,'(I3,A,I2.2,A,I2.2)') hours,':',minutes,':',seconds

    end function PRETTY_TIME

    SUBROUTINE BOX_MULLER (G1,G2)

        IMPLICIT NONE
  
        REAL*8 :: G1
        real*8, optional :: G2
        REAL*8 :: RAND(2)
        REAL*8 :: X1,X2,W
  
        DO
            CALL RANDOM_NUMBER(RAND)
  
            X1 = 2.D0*RAND(1)-1.D0
            X2 = 2.D0*RAND(2)-1.D0
            W  = X1*X1 + X2*X2
  
            IF (W.LT.1.D0) THEN
                EXIT
            END IF
        END DO
  
        W = SQRT(-2.D0*LOG(W)/W )
  
        G1 = X1*W  
        if (present(G2)) then
            G2 = X2*W
        end if
    END SUBROUTINE BOX_MULLER

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
        implicit NONE
        character(len=*), intent(in) :: filename
        integer, save :: count = 1
        integer :: ID

        open(UNIT = count, FILE = filename)
        ID = count
        count = count + 1;

    end function OPEN_NEW_FILE


    function linspace(a,b,N) result(V)
        real*8, intent(in) :: a, b
        integer, intent(in) :: N
        real*8, dimension(N) :: V
        integer :: i
        real*8 :: dv

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

    function interpol_scalar(coarse_x,coarse_y,fine_x) result(fine_y)
        real*8, intent(in) :: coarse_x(:), coarse_y(:), fine_x
        real*8 :: fine_y
        integer :: j, N

        if (size(coarse_x) /= size(coarse_y)) then
            print*, 'Coarse arrays must be of same length.'
        end if

        N = size(coarse_x)
        j = 1

        do while (coarse_x(j) <= fine_x .AND. j < N)
            j = j + 1
        end do
        
        if (j == 1) then
            fine_y = (coarse_y(2) - coarse_y(1))/(coarse_x(2) - coarse_x(1)) &
                 *(fine_x - coarse_x(1)) + coarse_y(1)
        else
            fine_y = (coarse_y(j) - coarse_y(j-1))/(coarse_x(j) - coarse_x(j-1)) &
                     *(fine_x - coarse_x(j-1)) + coarse_y(j-1)
        end if
    end function interpol_scalar

    function interpol_array(coarse_x,coarse_y,fine_x) result(fine_y)
        real*8, intent(in) :: coarse_x(:), coarse_y(:), fine_x(:)
        real*8  :: fine_y(size(fine_x))
        integer :: i 

        do i = 1,size(fine_x)
            fine_y(i) = interpol_scalar(coarse_x,coarse_y,fine_x(i))
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
            print*,'Invalid phase tag.'
            stop
        end if
    end function SVP_scalar

    function SVP_array(T,phase) result(es)
        character(1), intent(in)  :: phase
        real*8      , intent(in)  :: T(:)
        real*8                    :: es(size(T))
        integer*8                 :: i

        do i = 1,size(T)
            es(i) = SVP_scalar(T(i),phase)
        end do
    end function SVP_array

    function SVP_matrix(T,phase) result(es)
        character(1), intent(in)  :: phase
        real*8      , intent(in)  :: T(:,:)
        integer*8                 :: dims(2)
        real*8, allocatable       :: es(:,:)
        integer*8                 :: i

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

    subroutine sort(x)
        real*8 :: x(:)
        real*8 :: a
        integer :: i, check

        check = 1
        do while (check .ne. 0)
            check = 0
            do i = 1, size(x)-1
                if (x(i) > x(i+1)) then
                    a = x(i+1)
                    x(i+1) = x(i)
                    x(i) = a
                    check = 1
                end if
            end do
        end do
    end subroutine sort

end module functions

MODULE CCN
    implicit none
    real*8, allocatable :: CDF_TF(:,:)  !CDF for the freezing temperature
    real*8, allocatable :: CDF_DR(:,:)  !CDF for the Dry radii
    private :: SET_CDF_DRY, SET_CDF_TFREEZING, INIT_RANDOM_SEED

    contains

    subroutine SET_CDF_DRY
        use FUNCTIONS

        implicit none

        integer :: N_CDF
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

        implicit none
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
        R_DRY = interpol(CDF_DR(2,:),CDF_DR(1,:),R_DRY)
        R_DRY = R_DRY*1D-6 !Conversion from micro meters to meters
    end subroutine SAMPLE_DRY


    SUBROUTINE SET_CDF_TFREEZING

        !Sets the cumulative distribution function (CDF)
        !for the freezing temperature.
        use constants, only: pi
        IMPLICIT NONE

        REAL*8  :: T,TMIN,TMAX,DT
        REAL*8  :: TC,TCMIN,TCMAX 
        REAL*8  :: A,NS,DNS_DT
        REAL*8  :: PDF
        REAL*8  :: NSF,X
        INTEGER :: I,M

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

        ALLOCATE( CDF_TF(0:1,0:M) )

        CDF_TF(0,0) = TMIN
        CDF_TF(1,0) = 0.D0
        PDF = 0.D0

        DO I = 1,M
            T = TMIN + DT*I
            TC = T - 273.15D0
            IF (TC.GT.TCMAX) THEN
                NS = 0.D0
                DNS_DT = 0.D0
            ELSE IF ( (TC.LE.TCMAX).AND.(TC.GT.TCMIN) ) THEN
                NS = NSF(TC)
                DNS_DT = -0.517D0*NS
            ELSE IF ( TC.LT.TCMIN ) THEN
                NS = NSF(TCMIN)
                DNS_DT = 0.D0
            END IF

            CDF_TF(0,I) = T
            CDF_TF(1,I) = EXP(-A*NS)

            PDF = -A*DNS_DT*EXP(-A*NS)
        END DO
    END SUBROUTINE SET_CDF_TFREEZING


    SUBROUTINE SAMPLE_TFREEZING (TF)

        !Samples "N" values of freezing temperatures (stored in TF(:)) using
        !the inverse transform sampling method. It uses the pre-calculated  
        !cumulative distribution function (CDF_TF).


        IMPLICIT NONE

        REAL*8  :: TF(:), P
        REAL*8  :: RAND(1),DELTA
        INTEGER :: M
        INTEGER :: I,J,k

        call SET_CDF_TFREEZING
        call INIT_RANDOM_SEED

        P = 0.1 !Percentage of SDs with rd_insol
        k = ceiling(size(TF)*P)

        M = SIZE( CDF_TF, DIM = 2 ) - 1

        DO I = 1,size(TF)
            CALL RANDOM_NUMBER(RAND)
            J = 0
            DO WHILE (.NOT.((RAND(1).GE.CDF_TF(1,J)).AND.(RAND(1).LT.CDF_TF(1,J+1))))
                J = J + 1
            END DO
            if (I > k) then
                TF(I) = -40.D0 + 273.15D0
            else
                !Evaluates the temperature corresponding to CDF_TF = RAND(1) by linear interpolation
                DELTA = ( RAND(1) - CDF_TF(1,J) )/( CDF_TF(1,J+1) - CDF_TF(1,J) )
                TF(I) = CDF_TF(0,J) + ( CDF_TF(0,J+1) - CDF_TF(0,J) )*DELTA
            end if
        END DO
    END SUBROUTINE SAMPLE_TFREEZING

    SUBROUTINE INIT_RANDOM_SEED

        IMPLICIT NONE

        INTEGER :: I,N,CLOCK
        INTEGER, ALLOCATABLE :: SEED(:)

        CALL RANDOM_SEED(SIZE = N)
        ALLOCATE(SEED(N))
        CALL SYSTEM_CLOCK(COUNT=CLOCK)
        SEED = CLOCK + 37 * (/ (I - 1, I = 1, N) /)
        CALL RANDOM_SEED(PUT = SEED)
        DEALLOCATE(SEED)
    END SUBROUTINE INIT_RANDOM_SEED

END MODULE CCN

module SD_VARIABLES
    implicit NONE

    integer*8           :: N_sd, NPB
    real*8              :: n_drops =1.D9        !Droplets per unit mass of dry air[1/kg] (pristine: 1e8, poluted: 1e9)
    real*8, allocatable :: X(:,:),XP(:,:)       !Positions of droplets (X: current, XP: previous)
    real*8, allocatable :: R_dry(:), R_crit(:)      !Dry and Activation radius
    real*8, allocatable :: R(:,:)                   !Instantaneous radius
    real*8              :: init_ice_radius
    real*8, allocatable :: Q_k(:), T_k(:)         !Local mixing ratio and potentitial temperature
    real*8, allocatable :: w_prime(:)             !Local velocity fluctuations
    real*8, allocatable :: xi(:), S_crit(:)         !multiplicity and critical supersaturation
    real*8, allocatable :: T_Freeze(:)
    real*8, parameter   :: C_gamma = 1.5D-9, kappa = 0.61D0
    !Statistical Parameters
    real*8 :: R_mean,R_mean_ice,R_mean_water, sigma  !Average R and std deviation
    real*8 :: Activated_frac, init_ice_frac
    !Logical Attributes
    logical, allocatable :: Activated(:), Frozen(:)
    integer, allocatable :: phase_change(:) !0 for no change, 1 for freezing and -1 for melting
    logical :: im_freezing
    !Log da saturação local
    real*8, allocatable :: Sk_log(:)
end module SD_VARIABLES

module SD_FUNCTIONS
    use SD_VARIABLES
    implicit NONE

    private :: FREEZE

contains

    subroutine stats(status)
        character(len = *), optional, intent(in) :: status
        integer*8 :: j, ice_count, water_count
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
            write(1,*) '# R_mean	R_mean_ice	R_mean_water	sigma'
        end if
        write(1,*) R_mean, R_mean_ice, R_mean_water, sigma
        if (present(status) .and. status == 'close') then
               close(1)
               print*, 'Saved stats.dat'
        end if
    end subroutine stats


    subroutine WHICH_BOX(k,i,j)
        use GRID, only: DX,DZ

        integer*8, intent( in) :: k
        integer*8, intent(out) :: i,j
        ! OBS: Only works for rectangular and uniform grid.
        i = floor(X(k,1)/DX) + 1
        j = floor(X(k,2)/DZ) + 1
 
    end subroutine WHICH_BOX

    function WET_RADIUS(dry) result(R_wet)
        use THERMODYNAMIC_VAR, only: TEMP, RV
        use CONSTANTS, only: R_D, R_V
        use FUNCTIONS, only: interpol, svp
        use ENVIRONMENT

        integer*8 :: i, j, k
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
            call WHICH_BOX(k,i,j)
            T   = TEMP(i,j)
            Q   = RV(i,j)
            e   = Q/(Q + R_D/R_V)*interpol(Z_E,PR_E,X(k,2))
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
            end if
        end do
    end function WET_RADIUS

    subroutine INIT_SD_POSITIONS
        use PERIOD3D
        implicit none
        
        real*8    :: rn(2)
        integer*8 :: i

        call INIT_RANDOM_SEED
        do i = 1,N_sd
            call random_number(rn)
            X(i,1) = rn(1)*LX
            X(i,2) = rn(2)*LZ
        end do
    end subroutine INIT_SD_POSITIONS

    subroutine INIT_SD
        use CCN
        use THERMODYNAMIC_VAR, only: RV,TEMP
        use GRID, only: NX, NZ
        implicit none
        integer*8 :: k, i, j

        N_sd = NPB*NX*NZ
        allocate(X(N_sd,2),XP(N_sd,2))
        allocate(R_dry(N_sd), R_crit(N_sd))
        allocate(R(2,N_sd))
        allocate(xi(N_sd), S_crit(N_sd))
        allocate(Activated(N_sd), Frozen(N_sd))
        allocate(phase_change(N_sd))
        allocate(T_Freeze(N_sd))
        allocate(Q_k(N_sd), T_k(N_sd),w_prime(N_sd))
        allocate(Sk_log(N_sd))

        xi  = ceiling(n_drops/N_sd)
        w_prime = 0.D0
        
        call INIT_SD_POSITIONS 
        call SAMPLE_DRY(R_dry)
        R(1,:) = R_dry;
        !Calculate initial wet radii
        do k = 1,N_sd
            call which_box(k,i,j)
            Q_k(k) = RV(i,j)
            T_k(k) = TEMP(i,j)
        end do
        R(2,:) = WET_RADIUS(R_dry)        
        k = NINT(N_sd*init_ice_frac)
        call SAMPLE_TFREEZING(T_Freeze)
        Frozen(1:k)      = .true.
        Frozen(k+1:N_sd) = .false.
    end subroutine INIT_SD

    subroutine FREEZE(S,T,k)

        implicit none

        real*8   , intent(in) :: S, T
        integer*8, intent(in) :: k
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


    subroutine GROW_DROPLETS
        use CONSTANTS
        use FUNCTIONS
        use THERMODYNAMIC_VAR
        use STEPPER           , only: dt
        use ENVIRONMENT       , only: PR_E
        use ADVECTION         , only: L, eps
        use ENVIRONMENT       , only: Z_INV
        use OMP_LIB
        implicit none

        integer*8 :: j, m, n, count
        real*8, parameter  :: tol = 1.D-8
        real*8, parameter  :: Ct = 0.63D0 !Lasher-Trapp et al. 2005
        real*8, parameter  :: Ce = 0.89D0 !Schumann 1991
        real*8 :: error, next, F, F_prime, A
        real*8 :: D, tau, var!, psi
        real*8 :: e_k, S_k, dq_k, dq_k_F, dT_k, P, T ! Field localization parameters
        real*8 :: aw, Daw, Y, rd

        ! Variables for benchmarking. Remove later
        integer*8, save:: start, end, rate

        aw(y,rd)  = (y**(3.0/2) - rd**3)/(y**(3.0/2) - (rd**3)*(1 - kappa))
        Daw(y,rd) = 3.0/2*y**(1.0/2)/(y**(3.0/2) - (rd**3)*(1 - kappa))*(1 - aw(y,rd))
        
        call system_clock(count=start,count_rate=rate)


        if (L /= 0.D0) then
            tau = Ct*(L**2/eps)**(1.D0/3)         ! Turbulent relaxation time
            var = 2.D0/3*Ce*(L*eps)**(2.D0/3)     ! Variance of w
        end if
        
        !Solving the growth equation
        !$OMP PARALLEL private(P,T,e_k,S_k,A,D,error,Y,count,F,F_prime,next,dq_k,dq_k_F,dT_k,m,n) !Removed psi
        !$OMP DO
        do j = 1,N_sd

            !P = interpol(Z_E,PR_E,X(j,2)) !Pressure at the position of SD
            
            call WHICH_BOX(j,m,n)
            P = PR_E(n)!Pressure at the center of the box
            T = THETA(m,n)*EXNER(P)


            !Updating local mixing ratio and potential temperature
            if (L /= 0.D0) then
                dq_k     = -RHO_LIQ*N_DROPS_BOX(m,n)*4.0/3.0*pi*(R(2,j)**3 - R(1,j)**3)
                dq_k_F = phase_change(j)*N_DROPS_BOX(m,n)*RHO_LIQ*(4.D0/3.D0*pi)*R(1,j)**3
                if (Frozen(j)) then
                    dT_k = (-LS0*dq_k + (LS0 - LV0)*dq_k_F - G*(X(j,2)-XP(j,2)) )/CP_D !Removed "+ w_prime(j)*dt"
                else
                    dT_k = (-LV0*dq_k + (LS0 - LV0)*dq_k_F - G*(X(j,2)-XP(j,2)) )/CP_D !Removed "+ w_prime(j)*dt"
                end if
                !call SAMPLE_GAUSSIAN(psi,1.D0,0.D0)

                !Setting eps = 0 above inversion
                if (X(j,2) < Z_INV) then
                    Q_k(j) = Q_k(j) + (RV(m,n) - Q_k(j))*(dt/tau) + dq_k
                    T_k(j) = T_k(j) + (T - T_k(j))*(dt/tau) + dT_k
                    !w_prime(j) = w_prime(j)*exp(-dt/tau) + sqrt(var*(1 - exp(-2*dt/tau)))*psi
                else
                    Q_k(j)     = Q_k(j) + dq_k
                    T_k(j)     = T_k(j) + dT_k
                    !w_prime(j) = 0.D0
                end if
            else 
                !Conventional case: Homogeneous box with no fluctuation
                Q_k(j) = RV(m,n)
                T_k(j) = T
            end if

            !Saving old radii
            R(1,j) = R(2,j)
            
            e_k = P*Q_k(j)/(Q_k(j) + R_D/R_V)
            S_k = e_k/SVP(T_k(j),'L') - 1
            
            if (im_freezing) then
                call FREEZE(S_k,T_k(j),j)
            end if

            if (Frozen(j)) then
                S_k = e_k/SVP(T_k(j),'I') - 1
                A = 0.D0
                D = (R_V*T_k(j)*RHO_LIQ/D0/svp(T_k(j),'I') + LS0*RHO_LIQ/KT/T_k(j)*(LS0/R_V/T_k(j) - 1))**(-1)
            else
                A = C_gamma
                D = (R_V*T_k(j)*RHO_LIQ/D0/svp(T_k(j),'L') + LV0*RHO_LIQ/KT/T_k(j)*(LV0/R_V/T_k(j) - 1))**(-1)
            end if

            error = 1.D0
            Y = R(2,j)**2 !Initial guess for the Y = R² parameter
            count = 0
            do while (error > tol)
                F = 2*D*(S_k + 1 - aw(Y,R_dry(j))*exp(A/(Y**(1.D0/2))))*dt + R(2,j)**2 - Y
                F_prime = 2*D*exp(A/(Y**(1.D0/2)))*(aw(Y,R_dry(j))*A/2/(Y**(3.D0/2)) - Daw(Y,R_dry(j)))*dt - 1
                next = max(Y - F/F_prime,R_dry(j)**2)
                error = abs(F)
                Y = next
                count = count + 1
                if (count > 100) then
                    print*, 'Growth equation did not converge'
                    stop
                end if 
            end do
            R(2,j) = sqrt(Y)              
            !Sk_log(j) = S_k
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        call system_clock(end)
        print'(A,F6.2,A)', 'Updating the radii of all SDs took ', real(end-start)/rate, ' s.'
        call sleep(1)
    end subroutine GROW_DROPLETS

    subroutine GET_ANALITICAL_K(KK,GRAD_KK,k)
        use ENVIRONMENT, only: EDDY_MAX, Z_INV, TT
        use CONSTANTS  , only: pi
        integer*8, intent(in ) :: k
        real*8   , intent(out) :: KK, GRAD_KK(2)

        KK         = 0.5D0*EDDY_MAX*ERFC((X(k,2) - Z_INV)/(TT*sqrt(2.D0)))
        GRAD_KK(1) = 0.D0
        GRAD_KK(2) = -EDDY_MAX/(TT*sqrt(2.D0*pi))*exp(-(X(k,2) - Z_INV)**2/(2*TT**2))
    end subroutine GET_ANALITICAL_K

    subroutine ADVECTION_SD
        use ADVECTION, only: UXT, UXTP, UZT, UZTP
        use PERIOD3D,  only: LX, LZ
        use STEPPER,   only: DT, TIME, TIME_ADV
        use GRID,      only: DX, DZ
        use FUNCTIONS, only: interpol, BOX_MULLER
        use OMP_LIB
        integer*8 :: i, j, k, l
        real*8    :: U_SD(N_sd,2),UP_SD(N_sd,2)
        real*8    :: X_WALLS(2), Z_WALLS(2)
        real*8    :: UX_NOW(2), UZ_NOW(2)
        real*8    :: KKP, KK, GRAD_KK(2), ksi(2)

        XP = X
        
        !$OMP parallel private(i,j,l,X_WALLS,Z_WALLS,UX_NOW,UZ_NOW,KKP,KK,GRAD_KK,ksi)
        !$OMP DO
        do k = 1,N_sd
            call BOX_MULLER(ksi(1),ksi(2))
            do l = 0,1
                call WHICH_BOX(k,i,j)
                X_WALLS = (/(i-1)*DX, i*DX/)
                Z_WALLS = (/(j-1)*DZ, j*DZ/)
                
                UX_NOW(1)  = interpol(TIME_ADV,(/UXTP(i  ,j),UXT(i  ,j)/),TIME + l*DT)
                UX_NOW(2)  = interpol(TIME_ADV,(/UXTP(i+1,j),UXT(i+1,j)/),TIME + l*DT)
                UZ_NOW(1)  = interpol(TIME_ADV,(/UZTP(i,j  ),UZT(i,j  )/),TIME + l*DT)
                UZ_NOW(2)  = interpol(TIME_ADV,(/UZTP(i,j+1),UZT(i,j+1)/),TIME + l*DT)

                if (l==0) then
                    !First step of predictor-corrector
                    call GET_ANALITICAL_K(KKP,GRAD_KK,k)
                    UP_SD(k,1) = interpol(X_WALLS,UX_NOW,XP(k,1)) + GRAD_KK(1)
                    UP_SD(k,2) = interpol(Z_WALLS,UZ_NOW,XP(k,2)) + GRAD_KK(2)
                    X(k,:) = XP(k,:) + UP_SD(k,:)*DT + sqrt(2*KKP*DT)*ksi
                else
                    !Second step of predictor-corrector
                    call GET_ANALITICAL_K(KK,GRAD_KK,k)
                    U_SD(k,1) = interpol(X_WALLS,UX_NOW,X(k,1)) + GRAD_KK(1)
                    U_SD(k,2) = interpol(Z_WALLS,UZ_NOW,X(k,2)) + GRAD_KK(2)
                    X(k,:) = XP(k,:) + 0.5*(U_SD(k,:) + UP_SD(k,:))*DT + 0.5*(sqrt(2*KKP) + sqrt(2*KK))*sqrt(DT)*ksi
                end if                            

                ! Periodic boundary conditions on the sides
                if (X(k,1) > LX) then
                    X(k,1) = X(k,1) - LX
                else if (X(k,1) < 0) then
                    X(k,1) = X(k,1) + LX
                end if

                ! Reflective boundary conditions on top and bottom
                if (X(k,2) > LZ) then
                    X(k,2) = 2*LZ - X(k,2)
                else if (X(k,2) < 0) then
                    X(k,2) = -X(k,2)
                end if
            end do
        end do
        !$OMP END DO        
        !$OMP end parallel
    end subroutine ADVECTION_SD

end module SD_FUNCTIONS

module BOX
    use ADVECTION
    use THERMODYNAMIC_VAR
    use SD_FUNCTIONS, only: which_box
    implicit none

contains

    subroutine UPDATE_BOXES
        use CONSTANTS
        use FUNCTIONS
        use ENVIRONMENT
        use SD_VARIABLES, only: R, xi, N_sd, Frozen, phase_change, R_dry
        use GRID        , only: NX, NZ, DZ 
        use OMP_LIB

        integer*8 :: k, i, j
        real*8    :: DELTA_L(NX,NZ), DELTA_I(NX,NZ), DELTA_F(NX,NZ)
        real*8    :: P, Exn(NX,NZ)
        real*8    :: N_boxes

        DELTA_L = 0.D0
        DELTA_I = 0.D0
        DELTA_F = 0.D0

        RL      = 0.D0
        RI      = 0.D0

        N_boxes = NX*NZ
        do k = 1,N_sd
            call which_box(k,i,j)
            if (Frozen(k)) then
                DELTA_I(i,j) = DELTA_I(i,j) + xi(k)*(R(2,k)**3 - R(1,k)  **3)
                RI(i,j)      = RI(i,j)      + xi(k)*(R(2,k)**3 - R_dry(k)**3)
            else
                DELTA_L(i,j) = DELTA_L(i,j) + xi(k)*(R(2,k)**3 - R(1,k)  **3)
                RL(i,j)      = RL(i,j)      + xi(k)*(R(2,k)**3 - R_dry(k)**3)
            end if
            DELTA_F(i,j)     = DELTA_F(i,j) + xi(k)*(R(1,k)**3 - R_dry(k)**3)*phase_change(k)
        end do
        
        DELTA_I = DELTA_I * RHO_LIQ*(4.D0/3*pi)*N_boxes
        DELTA_L = DELTA_L * RHO_LIQ*(4.D0/3*pi)*N_boxes
        DELTA_F = DELTA_F * RHO_LIQ*(4.D0/3*pi)*N_boxes
        RI      = RI      * RHO_LIQ*(4.D0/3*pi)*N_boxes
        RL      = RL      * RHO_LIQ*(4.D0/3*pi)*N_boxes
    
        !Gathering value of Pressure and Exner function at grid box center 
        !from environmental pressure profile
        do j = 1,NZ
            P        = interpol(Z_E,PR_E,(real(j)-0.5)*DZ)
            Exn(:,j) = (P/P00)**(R_D/CP_D)
        end do
    
        !Updating potential temperature and vapor mixing ratio
        call GET_SCALAR_FORCES_TURB
        THETA = THETA + (LV0*DELTA_L + LS0*DELTA_I + DELTA_F*(LS0 - LV0))/CP_D/Exn + FTH_TURB_DT
        RV    = RV - DELTA_L - DELTA_I + FRV_TURB_DT
    
        !Updating  temperature
        TEMP = Exn*THETA
    end subroutine UPDATE_BOXES

    subroutine DROPS_PER_BOX(normalized)
        use SD_VARIABLES     , only: N_sd, xi
        use THERMODYNAMIC_VAR, only: DPB, N_DROPS_BOX
        use GRID             , only: NX, NZ

        integer*8           :: i, j, k
        logical, optional   :: normalized

        DPB = 0
        N_DROPS_BOX = 0;
        do k = 1,N_sd
            call which_box(k,i,j)
            DPB(i,j)         = DPB(i,j) + 1
            N_DROPS_BOX(i,j) = N_DROPS_BOX(i,j) + xi(k)*NX*NZ
        end do

        if (present(normalized) .and. normalized) then
            DPB = DPB/N_sd
        end if

    end subroutine DROPS_PER_BOX

    subroutine SAT_FIELD
        use GRID             , only: NX, NZ
        use ENVIRONMENT      , only: PR_E
        use THERMODYNAMIC_VAR, only: TEMP, RV, SAT
        use CONSTANTS        , only: R_D, R_V
        use FUNCTIONS        , only: SVP, interpol
        real*8     :: e(NX,NZ),P(NX,NZ)
        integer*8  :: j
 
        do j = 1,NZ
            P(:,j) = PR_E(j)
        end do         

        e   = P*RV/(RV + R_D/R_V)
        SAT = e/SVP(TEMP,'L') - 1        
    end subroutine SAT_FIELD

    SUBROUTINE GET_SCALAR_FORCES_TURB
  
        !Get the "forces" (i.e., the RHS of THE Eulerian (scalar) advection
        !equations) accounting for turbulent fluxes
       
        USE ENVIRONMENT, only: EDDY_DIFF
        USE THERMODYNAMIC_VAR
        USE STEPPER 
        USE GRID
     
        REAL*8, ALLOCATABLE ::  JX(:,:), JZ(:,:), DIV(:,:)
        integer*8 :: i
     
        IF ( .NOT.ALLOCATED(FTH_TURB_DT) ) THEN
           ALLOCATE ( FTH_TURB_DT(NX,NZ), FRV_TURB_DT(NX,NZ) )
        END IF
       
        ALLOCATE ( JX(NX,NZ), JZ(NX,NZ), DIV(NX,NZ) )
     
        !Turbulent flux of potential temperature THETA(:,:)
        CALL EVAL_GRADIENT_SCALAR ( THETA, JX, JZ )
        
        do i = 1,NX
            JX(i,:) = EDDY_DIFF(:)*JX(i,:)
            JZ(i,:) = EDDY_DIFF(:)*JZ(i,:)
        end do
        
        JZ(1:NX, 1) = 0.D0  !ZERO vertical turbulent flux boundary condition
        JZ(1:NX,NZ) = 0.D0  !ZERO vertical turbulent flux boundary condition
        
        CALL EVAL_DIVERGENCE ( DIV, JX, JZ )
        
        FTH_TURB_DT = DIV*DT

        !Turbulent flux of vapor mixing ratio QV(:,:)       
        CALL EVAL_GRADIENT_SCALAR ( RV, JX, JZ )
        
        do i = 1,NX
            JX(i,:) = EDDY_DIFF(:)*JX(i,:)
            JZ(i,:) = EDDY_DIFF(:)*JZ(i,:)
        end do
       
        JZ(1:NX, 1) = 0.D0  !ZERO vertical turbulent flux boundary condition
        JZ(1:NX,NZ) = 0.D0  !ZERO vertical turbulent flux boundary condition
        
        CALL EVAL_DIVERGENCE ( DIV, JX, JZ )
       
        FRV_TURB_DT = DIV*DT
     
        DEALLOCATE (JX,JZ,DIV)
    END SUBROUTINE GET_SCALAR_FORCES_TURB

    FUNCTION PER(I,N)
            
        INTEGER*8 :: PER, I, N
       
        IF (I.LT.1) THEN
           PER = I + N
        ELSEIF (I.GT.N) THEN
           PER = I - N
        ELSE
           PER = I
        END IF  
    END FUNCTION PER
     
    SUBROUTINE EVAL_GRADIENT_SCALAR (S,DSX,DSZ)
       
        !Evaluates the gradient of the scalar "S" at every point of the
        !primary grid. Periodic boundary conditions in the horizontal (x)
        !direction and no-flux BC in the vertical (z) direction.
     
        USE GRID

        REAL*8    :: S(NX,NZ), DSX(NX,NZ), DSZ(NX,NZ)
        INTEGER*8 :: IX, IZ
     
        DO IX = 1,NX
           DO IZ = 1,NZ
              !Horizontal gradient (periodic boundary conditions)
              DSX(IX,IZ) = 0.5D0/DX*( S(PER(IX+1,NX),IZ) - S(PER(IX-1,NX),IZ) )
             
              !Vertical gradient
              IF (IZ.EQ.1) THEN
                 DSZ(IX,IZ) = 0.5D0/DZ*( -3.D0*S(IX,IZ) + 4.D0*S(IX,IZ+1) - S(IX,IZ+2) )
              ELSE IF (IZ.EQ.NZ) THEN
                 DSZ(IX,IZ) = 0.5D0/DZ*(  3.D0*S(IX,IZ) - 4.D0*S(IX,IZ-1) + S(IX,IZ-2) )   
              ELSE
                 !Interior points
                 DSZ(IX,IZ) = 0.5D0/DZ*( S(IX,IZ+1) - S(IX,IZ-1) )
              END IF
           END DO
        END DO
    END SUBROUTINE EVAL_GRADIENT_SCALAR
     
    SUBROUTINE EVAL_DIVERGENCE (DIV,JX,JZ)
       
        !Evaluates the divergence of the vector (JX,JZ) at every point of
        !the primary grid. Periodic boundary conditions in the horizontal
        !(x) direction and no-flux in the vertical (z) direction.
       
        USE GRID

        REAL*8    :: DIV(NX,NZ), JX(NX,NZ), JZ(NX,NZ)
        REAL*8    :: AX, AZ
        INTEGER*8 :: IX, IZ
       
        DO IX = 1,NX
           DO IZ = 1,NZ
              !Horizontal contribution (periodic boundary conditions)
              AX = 0.5D0/DX*( JX(PER(IX+1,NX),IZ) - JX(PER(IX-1,NX),IZ) )
             
              !Vertical contribution
              IF (IZ.EQ.1) THEN
                 AZ = 0.5D0/DZ*( -3.D0*JZ(IX,IZ) + 4.D0*JZ(IX,IZ+1) - JZ(IX,IZ+2) )
              ELSE IF (IZ.EQ.NZ) THEN
                 AZ = 0.5D0/DZ*(  3.D0*JZ(IX,IZ) - 4.D0*JZ(IX,IZ-1) + JZ(IX,IZ-2) )
              ELSE
                 !Interior points
                 AZ = 0.5D0/DZ*( JZ(IX,IZ+1) - JZ(IX,IZ-1) )
              END IF
     
              DIV(IX,IZ) = AX + AZ
           END DO
        END DO
    END SUBROUTINE EVAL_DIVERGENCE
end module BOX

module ADV_FUNCTIONS

    implicit none

    contains

    subroutine CHECK_ADV
        use ADVECTION, only: ADVECT
        use STEPPER  , only: TIME, DT_ADV, DT
        INTEGER*8, SAVE :: IT = 0, IT_OLD = 0
        
        ! Workarround for float precision
        if (abs(TIME/DT_ADV - nint(TIME/DT_ADV)) < DT/2) then 
            IT = nint(TIME/DT_ADV)
        else
            IT = floor(TIME/DT_ADV)
        end if
  
        if (IT == IT_OLD) then
            ADVECT = .FALSE.
        else    
            IT_OLD = IT
            ADVECT = .TRUE.
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

    SUBROUTINE GET_VELOCITIES

        !Gets the velocities "UXT" and "UZT" in the dual grid and the flow
        !velocities "UX" and "UZ" (in the PRIMARY grid) by interpolation
        
        !Derived numerically from the prescribed INCOMPRESSIBLE streamfunction
  
        use ADVECTION
        use GRID
        use STEPPER, only: TIME, DT_ADV, TIME_ADV
  
        INTEGER*8       :: I,K
  
        IF (.NOT.ALLOCATED(UXT)) THEN
            ALLOCATE (UXT(NXP,NZ ))
            ALLOCATE (UZT(NX ,NZP))
  
            UXT  = 0.D0
            UZT  = 0.D0
  
            ALLOCATE (UX(NX,NZ))
            ALLOCATE (UZ(NX,NZ))
  
            ALLOCATE (UXTP(NXP,NZ ))
            ALLOCATE (UZTP(NX ,NZP))

            CALL HDF5_READ_VELOCITIES(velocities_file)
        END IF
  
        UXTP = UXT
        UZTP = UZT
  
        !Get transport velocities UXT and UZT from HDF5 file
        CALL HDF5_READ_VELOCITIES(velocities_file)
  
        !Flow velocity "UX" and "UZ" in the primary grid by interpolation
  
        DO I = 1,NX
            DO K = 1,NZ
                UX(I,K) = 0.5D0*( UXT(I,K) + UXT(I+1,K  ) )
                UZ(I,K) = 0.5D0*( UZT(I,K) + UZT(I  ,K+1) )
            END DO
        END DO
  
        CALL GET_ADVECTIVE_VELOCITIES

        TIME_ADV = (/TIME, TIME + DT_ADV/)
  
    END SUBROUTINE GET_VELOCITIES

    SUBROUTINE GET_ADVECTIVE_VELOCITIES

        !Gets the velocities "UXA" and "UZA" in the dual grid (at time level
        !"n+1/2") for the advection by MPDATA.
        
        !*** Uses "UXT" and "UZT" and their values "UXTP" and "UZTP" in the previous time level ***
  
        use ADVECTION
        use GRID
        use STEPPER, ONLY: DT_ADV
        use ENVIRONMENT, ONLY: RHO
  
        REAL*8  :: RHO_BAR
        INTEGER*8 :: IX,IZ
  
        IF (.NOT.ALLOCATED(UXA)) THEN
            ALLOCATE (UXA(NXP,NZ ))
            ALLOCATE (UZA(NX ,NZP))
        END IF
  
        !Horizontal velocities (applies periodic boundary conditions)
        
        !Extrapolating (Abade)
        ! DO IZ = 1,NZ
        !     DO IX = 1,NXP
        !         UXA(IX,IZ) = ( 1.5D0*UXT(IX,IZ) - 0.5D0*UXTP(IX,IZ) )*DT_ADV/DX
        !     END DO
        ! END DO
        
        !Interpolating
        UXA = 0.5D0*(UXT + UXTP)*DT_ADV/DX

        !Vertical velocities
        
        !Extrapolating (Abade)
        ! DO IX = 1,NX
        !     DO IZ = 1,NZP
        !         UZA(IX,IZ) = ( 1.5D0*UZT(IX,IZ) - 0.5D0*UZTP(IX,IZ) )*DT_ADV/DZ
        !     END DO
        ! END DO
  
        !Interpolating
        UZA = 0.5D0*(UZT + UZTP)*DT_ADV/DZ
        ! * * *
  
        !Multiplying by the base-state dry density RHO
  
        DO IZ = 1,NZ
            DO IX = 1,NXP
                UXA(IX,IZ) = RHO(IZ)*UXA(IX,IZ)
            END DO
        END DO
  
  
        DO IZ = 1,NZP
  
            IF (IZ.EQ.1) THEN
                RHO_BAR = 0.5D0*( 3.D0*RHO(1) - RHO(2) )
            ELSE IF (IZ.EQ.NZP) THEN
                RHO_BAR = 0.5D0*( 3.D0*RHO(NZ) - RHO(NZ-1) )
            ELSE
                RHO_BAR = 0.5D0*( RHO(IZ-1) + RHO(IZ) )
            END IF
  
            DO IX = 1,NX
                UZA(IX,IZ) = RHO_BAR*UZA(IX,IZ)
            END DO
  
        END DO
  
        ! * * *
  
    END SUBROUTINE GET_ADVECTIVE_VELOCITIES

    SUBROUTINE ADVECTION_MPDATA

        !2-D advection of "theta" and "RV" fields
        !To be used in the Super-droplet scheme
    
        use THERMODYNAMIC_VAR
        use ADVECTION
        use GRID
    
        INTEGER :: IORD,ISOR,NONOS,IDIV !MPDATA parameters
  
        IORD  = 2
        ISOR  = 0
        NONOS = 1
        IDIV  = 0
  
        CALL MPDATA2D (UXA,UZA,THETA,GAC,NX,NZ,IORD,ISOR,NONOS,IDIV,1)
        CALL MPDATA2D (UXA,UZA,   RV,GAC,NX,NZ,IORD,ISOR,NONOS,IDIV,2)
  
    END SUBROUTINE ADVECTION_MPDATA

end module ADV_FUNCTIONS