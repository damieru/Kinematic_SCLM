
!********************************************************************************
!********************************************************************************
!********************************************************************************

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
  IMPLICIT NONE

  REAL*8  :: T,TMIN,TMAX,DT
  REAL*8  :: TC,TCMIN,TCMAX
  REAL*8  :: A,NS,DNS_DT
  REAL*8  :: PI
  REAL*8  :: PDF
  REAL*8  :: NSF,X
  INTEGER :: I,M

  NSF(X) = EXP( - 0.517D0*X + 8.934D0 ) ![m^{-2}]
  PI = ACOS(-1.D0)
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

!********************************************************************************
!********************************************************************************
!********************************************************************************

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

!********************************************************************************
!********************************************************************************
!********************************************************************************

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

!********************************************************************************
!********************************************************************************
!********************************************************************************

END MODULE CCN
