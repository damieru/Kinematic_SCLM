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
   integer(HSIZE_T), dimension(3) :: stride  = [1,1,1]
   integer(HSIZE_T), dimension(3) :: block   = [1,1,1]
   integer(HSIZE_T), dimension(2) :: stride2 = [1,1]
   integer(HSIZE_T), dimension(2) :: block2  = [1,1]
end module

module CUSTOM_TYPES

   type :: box_address
      integer   :: i, j
   end type

   type :: grid_node
      type(box_address) :: srd_box(4)
      real*8            :: pos(2)
   end type

   type :: particle_list
      integer  , allocatable :: p_list(:)
   end type

end module CUSTOM_TYPES

MODULE VELOCITY_MODES

  INTEGER :: MK,NK                         !Number of modes in "x" and
                                           !"z" directions

  REAL*8, ALLOCATABLE :: KM(:),KN(:)       !Wavenumbers in "x" and "z" directions

  REAL*8, ALLOCATABLE :: DM(:),CN(:)       !Fourier coefficients

  REAL*8, ALLOCATABLE :: AMN(:,:),BMN(:,:) !Random coefficients

  REAL*8 :: LLX,LLZ                        !Periods in "x" and "z" directions

END MODULE

MODULE ENVIRONMENT

  REAL*8, ALLOCATABLE :: Z_E(:), Exn(:)

  REAL*8, ALLOCATABLE :: PR_E(:),TH_E(:),RV_E(:),T_E(:)

  REAL*8, ALLOCATABLE :: RHO(:) !Base state density profile

  REAL*8, ALLOCATABLE :: EDDY_DIFF(:) !Eddy diffusivity profile

  REAL*8 :: TT, EDDY_MAX !Transition thickness and maximum value of the eddy diffusivity profile

  REAL*8 :: Z_INV = 1.D6 !Vertical position of the inversion: here
                         !initialized with an unrealistically large
                         !value

END MODULE

MODULE GRID
  
  use CUSTOM_TYPES

  REAL*8    :: LX , LZ
  REAL*8    :: DX , DZ
  integer   :: NX , NZ
  integer   :: NXP, NZP

  real*8             , allocatable :: Z_nodes(:)
  type(grid_node)    , allocatable :: cell_nodes(:,:)
  type(particle_list), allocatable :: boxes(:,:)

END MODULE

module SD_VARIABLES
   use CUSTOM_TYPES
   implicit NONE

   integer                        :: N_sd, DPB_input               ! Total number of superdroplets and number of SDs per box (input)
   real*8                         :: init_ice_radius               ! Initial ice radius for ice particles
   real*8                         :: n_drops =1.D9                 ! Droplets per unit mass of dry air[1/kg] (pristine: 1e8, poluted: 1e9)
   real*8, target   , allocatable :: X(:,:),XP(:,:)                ! Positions of droplets (X: current, XP: previous)
   type(box_address), allocatable :: SD_box(:)                     ! Box addresses for each SD
   real*8           , allocatable :: R_dry(:), R_crit(:)           ! Dry and Activation radius
   real*8           , allocatable :: R(:,:)                        ! Instantaneous radius
   real*8           , allocatable :: Q_k(:), TH_k(:)                ! Local mixing ratio and potentitial temperature
   real*8           , allocatable :: u_prime(:,:),u_prime_old(:,:) ! Local velocity fluctuations
   real*8           , allocatable :: xi(:), S_crit(:)              ! multiplicity and critical supersaturation
   real*8           , allocatable :: T_Freeze(:)                   ! Freezing temperature
   real*8           , parameter   :: C_gamma=1.5D-9, kappa=0.61D0  ! Droplet growth constants
  
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

MODULE THERMODYNAMIC_VAR

  REAL*8   , ALLOCATABLE :: THETA(:,:), TEMP(:,:)
  REAL*8   , ALLOCATABLE :: RV(:,:), RL(:,:), RI(:,:), SAT(:,:)
  REAL*8   , ALLOCATABLE :: GAC(:,:)
  REAL*8   , ALLOCATABLE :: N_DROPS_BOX(:,:)
  integer  , ALLOCATABLE :: DPB(:,:)

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

  REAL*8, ALLOCATABLE :: PSI(:,:), PSIP(:,:)   !Stream function

  !Flow velocity in the primary grid

  REAL*8, ALLOCATABLE :: UX(:,:), UZ(:,:)

  !Flow velocities at the interfaces of the primary grid
  REAL*8, ALLOCATABLE :: UXT (:,:), UZT (:,:)
  REAL*8, ALLOCATABLE :: UXTP(:,:), UZTP(:,:)

  !Flow velocities at the nodes
  real*8, allocatable :: UXN (:,:), UZN (:,:)

  !Characteristic eddy size and TKE dissipation rate
  real*8 :: L, eps_input
  real*8, allocatable :: eps(:) ! A value of eps for each height.
  real*8, allocatable :: z_eps(:)

  !Stochastic boolean to determine which microphysics model to be used
  logical :: stochastic_micro_phys

  ! Stochastic model constants
  real*8, parameter   :: C_zero = 2.D0, C_one = 0.5D0 + 0.75D0*C_zero 

  !Relaxation frequency omega
  real*8, allocatable :: omega(:)
  

  ! Tells if advected variables will be updated
  logical :: ADVECT

END MODULE ADVECTION

MODULE STEPPER

  REAL*8  :: TIME, T_MAX, OUT_PER
  REAL*8  :: DT

  REAL*8  :: DT_ADV
  real*8 :: TIME_ADV(2)
  INTEGER :: N_STEPS

END MODULE STEPPER

MODULE IO_PARAMETERS

  character(len=100) :: input_file, output_file, output_folder, eps_file_name
  character(len=  3) :: field ! Input variable. 'VEL' - velocities, 'PSI' - stream function
   
END MODULE IO_PARAMETERS


