module HDF5_VARIABLES
    use HDF5
    integer(HID_T) :: file_id, dset_id, dspace_id !Identifiers
    integer(HID_T) :: memspace
    integer :: error !Error flag
    integer :: rnk = 3
    integer(HSIZE_T), dimension(3) :: data_dims !Data dimensions
    integer(HSIZE_T), dimension(3) :: dimsm     !Subset Data dimensions
    integer(HSIZE_T), dimension(3) :: count, offset
    integer(HSIZE_T), dimension(2) :: count2, offset2
    integer(HSIZE_T), dimension(2) :: dimsm2
    integer(HSIZE_T), dimension(3) :: stride  = (/1,1,1/)
    integer(HSIZE_T), dimension(3) :: block   = (/1,1,1/)
    integer(HSIZE_T), dimension(2) :: stride2 = (/1,1/)
    integer(HSIZE_T), dimension(2) :: block2  = (/1,1/)

end module HDF5_VARIABLES

MODULE VELOCITY_MODES

  INTEGER :: MK,NK                         !Number of modes in "x" and
                                           !"z" directions

  REAL*8, ALLOCATABLE :: KM(:),KN(:)       !Wavenumbers in "x" and "z" directions

  REAL*8, ALLOCATABLE :: DM(:),CN(:)       !Fourier coefficients

  REAL*8, ALLOCATABLE :: AMN(:,:),BMN(:,:) !Random coefficients

  REAL*8 :: LLX,LLZ                        !Periods in "x" and "z" directions

  REAL*8 :: DZ_INV_W                       !Vertical velocity
                                           !overshoot above inversion

END MODULE VELOCITY_MODES

MODULE ENVIRONMENT

  REAL*8, ALLOCATABLE :: Z_E(:)

  REAL*8, ALLOCATABLE :: PR_E(:),TH_E(:),RV_E(:),T_E(:)

  REAL*8, ALLOCATABLE :: RHO(:) !Base state density profile

  REAL*8, ALLOCATABLE :: EDDY_DIFF(:) !Eddy diffusivity profile

  REAL*8 :: TT, EDDY_MAX !Transition thickness and maximum value of the eddy diffusivity profile

  REAL*8 :: Z_INV = 1.D6 !Vertical position of the inversion: here
                         !initialized with an unrealistically large
                         !value

END MODULE ENVIRONMENT

MODULE PERIOD3D

  REAL*8 :: LX,LY,LZ

END MODULE PERIOD3D

MODULE GRID

  INTEGER*8 :: NX,NZ
  INTEGER*8 :: NXP,NXM,NZP,NZM

  REAL*8    :: DX,DY,DZ
END MODULE GRID

MODULE THERMODYNAMIC_VAR

  REAL*8, ALLOCATABLE :: THETA(:,:),TEMP(:,:)

  REAL*8, ALLOCATABLE :: FTH_TURB_DT(:,:), FRV_TURB_DT(:,:)

  REAL*8, ALLOCATABLE :: RV(:,:),RL(:,:),RI(:,:),SAT(:,:)

  REAL*8, ALLOCATABLE :: GAC(:,:)

  REAL*8, ALLOCATABLE :: DPB(:,:)

  REAL*8, ALLOCATABLE :: N_DROPS_BOX(:,:)

END MODULE THERMODYNAMIC_VAR

MODULE CONSTANTS

  REAL*8, parameter :: pi = acos(-1.D0)
  REAL*8, parameter :: G  = 9.8067D0 !m s^{-2}

  REAL*8, parameter :: R_D = 287.1D0
  REAL*8, parameter :: R_V = 461.4D0

  REAL*8, parameter :: CP_D = 1005.D0
  REAL*8, parameter :: CP_V = 1850.D0

  REAL*8, parameter :: C_ICE = 1960.D0
  REAL*8, parameter :: C_LIQ = 4218.D0

  REAL*8, parameter :: RHO_ICE = 917.D0
  REAL*8, parameter :: RHO_LIQ = 1000.D0

  REAL*8, parameter :: NU_AIR = 14.D-6  !Air kinematic viscosity

  !Properties at the triple point

  REAL*8, parameter :: T0  = 273.16D0
  REAL*8, parameter :: LV0 = 2.5D6
  REAL*8, parameter :: LS0 = 2.83D6
  REAL*8, parameter :: LM0 = 3.3D5
  REAL*8, parameter :: ES0 = 611.73D0 ![Pa]

  !Reference pressure for potential temperature

  REAL*8, parameter :: P00 = 1.D5 ! = 1000 mb

  !Water vapor and heat diffusivities

  REAL*8, parameter :: D0 = 2.21D-5   ![m^{2} s^{-1}]
  REAL*8, parameter :: KT = 2.4D-2    ![m^{2} s^{-1}]

  !Surface tension (Kelvin) coefficient

  REAL*8, parameter :: A_SIGMA = 4.33D-7

END MODULE CONSTANTS

MODULE ADVECTION

  REAL*8, ALLOCATABLE :: UXA(:,:),UZA(:,:)

  REAL*8, ALLOCATABLE :: PSI(:,:)      !Stream function

  !Flow velocity in the primary grid

  REAL*8, ALLOCATABLE :: UX(:,:), UZ(:,:)

  !Flow velocities at the interfaces of the primary grid

  REAL*8, ALLOCATABLE :: UXT (:,:), UZT (:,:)
  REAL*8, ALLOCATABLE :: UXTP(:,:), UZTP(:,:)

  !Typical eddy size and TKE dissipation rate
  real*8 :: L, eps

  !Path to velocities file
  character(100) :: velocities_file

  ! Tells if advected variables will be updated
  logical :: ADVECT

  ! * * *

  INTEGER :: FLAG_TRANSIENT_FLOW      ! = 1 -> the prescribed flow
                                      ! changes in time. Otherwise,
                                      ! the flow is frozen

END MODULE ADVECTION

MODULE STEPPER

  REAL*8  :: TIME, T_MAX, OUT_PER
  REAL*8  :: DT

  REAL*8  :: DT_ADV
  real*8 :: TIME_ADV(2)
  INTEGER :: N_STEPS

END MODULE STEPPER

MODULE OUTPUT

  CHARACTER(100) :: FILE_NAME
  CHARACTER(100) :: FILE_PLACE

END MODULE OUTPUT


