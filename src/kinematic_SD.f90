PROGRAM MAIN
   use STEPPER
   use FUNCTIONS
   use ADV_FUNCTIONS 
   use SD_FUNCTIONS
   use BOX
   use CIC_ROUTINES
   !For Optimization and debugging only
   use profiler
   
   IMPLICIT NONE

   integer :: i
   !Debugging
   logical :: acabou

   !Initialization -----------------------------------------------------
   TIME = 0.D0
   CALL PRINT_REMARKS
   CALL READIN                              !Reads the input parameters
   CALL INIT_SG_PARAMETERS                  !Builds the profiles of epsilon and omega
   CALL GET_ENVIRONMENT_PROFILES_DYCOMS_II  !Sets up the environmental base state with a strong temperature inversion at z_inv
   CALL INIT_THERMODYNAMIC_VAR_SD           !Initializes the thermodynamic fields using environmental profiles
   CALL INIT_RANDOM_SEED                    !For the random sub-grid velocity fluctuations
   CALL INIT_SD
   CALL INIT_CIC_NODES
   CALL UPDATE_BOXES
   CALL SAT_FIELD
   CALL GET_VELOCITIES
   CALL HDF5_CREATE_FILE
   CALL HDF5_SAVE_RESULTS(OUT_PER)
   call progress_bar(TIME, T_MAX)
   
   call tictoc(time_unit='ms'); acabou = .false.
   DO i = 2,N_STEPS
      acabou = (i == N_STEPS)
      TIME = TIME + DT
      CALL CHECK_ADV
      !call tictoc(section_name='CHECK_ADV')
      CALL ADVECTION_MPDATA   ! Advects thermodynamic fields using Euler scheme
      !call tictoc(section_name='MPDATA')
      CALL GET_VELOCITIES     ! Obtains velocity field from stream function file
      !call tictoc(section_name='GET_VEL')
      CALL ADVECTION_SD       ! Transport of super-droplets by the flow
      !call tictoc(section_name='ADV_DROPLETS')
      CALL GROW_DROPLETS      ! Solves growth equations
      !call tictoc(section_name='GROWTH_EQ')
      CALL UPDATE_BOXES       ! Updates grid boxes mean properties after condensation
      !call tictoc(section_name='UPDATE_BOXES')
      CALL SAT_FIELD          ! Diagnostic only
      !call tictoc(section_name='SAT_FIELD')
      call HDF5_SAVE_RESULTS(OUT_PER)
      !call tictoc(section_name='SAVE_RES',finish=acabou)

      if (mod(nint(TIME/DT),10) == 0) call progress_bar(TIME, T_MAX)
   END DO
   !call system('cat tictoc.txt')
END PROGRAM MAIN

SUBROUTINE PRINT_REMARKS

  IMPLICIT NONE

  WRITE(*,*)
  WRITE(*,*) "* PLEASE NOTE **************************************************"
  WRITE(*,*)
  WRITE(*,*) "Uses the PINSKY ET AL. (JAS,2008) turbulent-like ABL flow"
  WRITE(*,*)
  WRITE(*,*) "Sets a constant value for the dry air density RHO in the spirit"
  WRITE(*,*) "of the Boussinesq approx."
  WRITE(*,*)
  WRITE(*,*) "****************************************************************"
  WRITE(*,*)

  CALL SLEEP(1) !Waits 1 second

END SUBROUTINE PRINT_REMARKS

SUBROUTINE READIN

   use GRID
   use STEPPER
   use ADVECTION
   use IO_PARAMETERS
   use VELOCITY_MODES, only: MK,NK,LLX,LLZ
   use SD_VARIABLES  , only: im_freezing, init_ice_frac,DPB_input
   use ENVIRONMENT   , only: Z_INV
   use FUNCTIONS     , only: linspace
   

   IMPLICIT NONE   

   character(len=100) :: file_line

   ! GRID ------------------------------------------------------------------------
   file_line = get_line(); READ(file_line,*) NX, NZ
   file_line = get_line(); READ(file_line,*) DX, DZ

   NXP = NX + 1
   NZP = NZ + 1

   LX = DX*NX
   LZ = DZ*NZ

   Z_nodes = linspace(0.D0,LZ,NZP)


   ! VELOCITY MODES --------------------------------------------------------------
   file_line = get_line(); READ(file_line,*) MK, NK

   LLX = LX
   LLZ = 2.D0*LZ

   WRITE(*,'(/A,F6.2)') "Turbulence resolution length scale: ", SQRT( LLX*LLZ/REAL(MK*NK) )

   ! STEPPER ---------------------------------------------------------------------
   file_line = get_line(); READ(file_line,*) DT
   file_line = get_line(); READ(file_line,*) DT_ADV
   file_line = get_line(); READ(file_line,*) T_MAX
   file_line = get_line(); READ(file_line,*) OUT_PER

   N_STEPS = NINT(T_MAX/DT) + 1
   TIME_ADV = [0.D0, DT_ADV]

   ! INVERSION -------------------------------------------------------------------
   file_line = get_line(); READ(file_line,*) Z_INV    !Thermodynamic initial inversion

   ! Number of super droplets per box
   file_line = get_line(); READ(file_line,*) DPB_input

   !Eddy size and TKE dissipation rate--------------------------------------------
   file_line = get_line(); READ(file_line,*) L

   ! Uses the deterministic model if L = 0
   if (L == 0.D0) then
      stochastic_micro_phys = .false.
   else
      stochastic_micro_phys = .true.
   end if

   ! Immersion freezing ON / OFF
   file_line = get_line(); READ(file_line,*) im_freezing

   ! Initial ice fraction
   file_line = get_line(); READ(file_line,*) init_ice_frac

   ! INPUT PARAMETERS ------------------------------------------------------------
   file_line = get_line(); READ(file_line,*) field !REMOVE THIS
   file_line = get_line(); READ(file_line,*) input_file
   file_line = get_line(); READ(file_line,*) eps_file_name

   ! OUTPUT PARAMETERS -----------------------------------------------------------
   file_line = get_line(); READ(file_line,*) output_folder
   file_line = get_line(); READ(file_line,*) output_file

   WRITE(*,'(/A,I6.3)') "NUMBER OF GRID BOXES: ", NX*NZ

   contains

   function get_line() result(line)
      ! Returns only the lines that don't start with '#'
      ! and are not empty. Ignores text after ":".
      character(len=100) :: probe, line
      integer :: error, i
      read(*,'(A)',iostat=error) probe

      do while ((probe(1:1) == '#' .or. len_trim(probe) == 0) .and. error == 0) 
         read(*,'(A)',iostat=error) probe 
      end do   

      do i = 1,len_trim(probe)
         if (probe(i:i)==':') exit
      end do

      line = trim(probe(1:i-1))
   end function get_line

END SUBROUTINE READIN

SUBROUTINE INIT_THERMODYNAMIC_VAR_SD

  use GRID
  use ENVIRONMENT
  use THERMODYNAMIC_VAR

  IMPLICIT NONE

  REAL*8  :: RHO_BAR
  integer :: IZ

  IF ( .NOT.(ALLOCATED(THETA)) ) THEN
     ALLOCATE (      THETA(NX,NZ) )
     ALLOCATE (       TEMP(NX,NZ) )
     ALLOCATE (         RV(NX,NZ) )
     ALLOCATE (         RL(NX,NZ) )
     ALLOCATE (         RI(NX,NZ) )
     ALLOCATE (        GAC(NX,NZ) )
     ALLOCATE (        DPB(NX,NZ) )
     ALLOCATE (N_DROPS_BOX(NX,NZ) )
     ALLOCATE (        SAT(NX,NZ) )     
  END IF

  !-----------------------------------------------------------------
  !Initializes the thermodynamic fields using the environmental
  !profile. The initial fields are assumed to be horizontally
  !homogeneous.

  DO IZ = 1,NZ
     THETA  (1:NX,IZ) = TH_E(IZ)
     TEMP   (1:NX,IZ) =  T_E(IZ)
     RV     (1:NX,IZ) = RV_E(IZ)
     RL     (1:NX,IZ) = 0.D0
     RI     (1:NX,IZ) = 0.D0
  END DO

  !******************** A T T E N T I O N **************************
  !*****************************************************************
  !For Boussinesq approximation as in Pinsky et al.
  !Sets a constant value for dry air density distribution

  RHO_BAR = sum(RHO)/NZ !Or simply: RHO_BAR = RHO(1)
  RHO     = RHO_BAR

  !*****************************************************************
  !*****************************************************************

  DO IZ = 1,NZ
     GAC(1:NX,IZ) = RHO(IZ)
  END DO

END SUBROUTINE INIT_THERMODYNAMIC_VAR_SD

SUBROUTINE INIT_RANDOM_SEED

  IMPLICIT NONE

  INTEGER :: I,N,CLOCK
  INTEGER, ALLOCATABLE :: SEED(:)

  CALL RANDOM_SEED(SIZE = N)
  ALLOCATE(SEED(N))

  CALL SYSTEM_CLOCK(COUNT=CLOCK)

  SEED = CLOCK + 37 * [ (I - 1, I = 1, N) ]

  CALL RANDOM_SEED(PUT = SEED)

  DEALLOCATE(SEED)

END SUBROUTINE INIT_RANDOM_SEED

SUBROUTINE HDF5_CREATE_FILE

   use HDF5_VARIABLES
   use GRID
   use STEPPER, ONLY: T_MAX, OUT_PER
   use SD_VARIABLES, only: N_sd
   use IO_PARAMETERS
   
   integer :: NT
   character(100) :: file_path
   logical :: f_exists

   inquire(file=output_folder, exist=f_exists)
   if (.NOT.f_exists) then 
       call system('mkdir '//output_folder)
   end if

   file_path = trim(output_folder)//'/'//trim(output_file)

   !inquire(file=file_path, exist=f_exists)
   !do while (f_exists)
   !    write(f_index,'(i2.2)') i
   !    file_path = trim(output_folder)//'/'//'('//trim(f_index)//')'//trim(output_file)
   !    inquire(file=file_path, exist=f_exists)
   !end do

   NT = NINT(T_MAX/OUT_PER) + 1 ! Number of saved instants
   call h5open_f(error) !Opens HDF5 Fortran interface
   call h5fcreate_f(file_path, H5F_ACC_TRUNC_F, file_id, error) !Creates file

   call h5screate_simple_f(1, int([NZ],8), dspace_id, error) !Creates dataspace
   !Creating RHO dataset
   call h5dcreate_f(file_id, 'RHO', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)

   call h5screate_simple_f(rnk, int([NX,NZ,NT],8), dspace_id, error) !Creates dataspace
   !Creating THETA dataset
   call h5dcreate_f(file_id, 'THETA', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)
   !Creating RV dataspace and dataset
   call h5dcreate_f(file_id, 'RV', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)
   !Creating RL dataspace and dataset
   call h5dcreate_f(file_id, 'RL', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)
   !Creating RI dataspace and dataset
   call h5dcreate_f(file_id, 'RI', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)
   !Creating DPB dataspace and dataset
   call h5dcreate_f(file_id, 'DPB', H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)
   !Creating SAT dataspace and dataset
   call h5dcreate_f(file_id, 'SAT', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)

   !Creating dataspace for UXN and UXZ (NXP by NZP)
   call h5screate_simple_f(rnk, int([NXP,NZP,NT],8), dspace_id, error)
   !Creating UXN dataspace and dataset
   call h5dcreate_f(file_id, 'UXN', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)
   !Creating UZN dataspace and dataset
   call h5dcreate_f(file_id, 'UZN', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)

   call h5screate_simple_f(rnk, int([N_sd, 2, NT],8), dspace_id, error) !Creates dataspace
   !Creating X_SD dataspace and dataset
   call h5dcreate_f(file_id, 'X_SD', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)

   call h5screate_simple_f(2, int([N_sd,NT],8), dspace_id, error) !Creates dataspace
   !Creating X_SD dataspace and dataset
   call h5dcreate_f(file_id, 'R_SD', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)

   call h5screate_simple_f(1, int([NT],8), dspace_id, error) !Creates dataspace
   !Creating X_SD dataspace and dataset
   call h5dcreate_f(file_id, 'TIME', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)

   !Closing file and interface
   call h5fclose_f(file_id,error)
   call h5close_f(error)

END SUBROUTINE HDF5_CREATE_FILE

SUBROUTINE HDF5_SAVE_RESULTS (NSEC)

    !Prints the output every "NSEC" minutes
    use HDF5_VARIABLES
    use STEPPER
    use ADVECTION, ONLY: UXN,UZN
    use THERMODYNAMIC_VAR, ONLY: THETA, RV, RL, RI, DPB, SAT
    use ENVIRONMENT, ONLY: RHO
    use GRID
    use SD_VARIABLES, ONLY: N_sd, X, R
    use IO_PARAMETERS
    IMPLICIT NONE

    REAL*8  :: NSEC
    INTEGER :: IT,PERIOD
    CHARACTER(100) :: FILE_PATH
    logical, save  :: first_call = .true.
    !Debugging only
    integer(HSIZE_T) :: current_dims(3), max_dims(3)

    FILE_PATH = TRIM(output_folder)//'/'//TRIM(output_file)

    IT = NINT( TIME/DT )
    PERIOD = NINT( NSEC/DT )

    IF ( MOD(IT,PERIOD).EQ.0 ) THEN

        call h5open_f(error) !Open fortran HDF5 interface
        call h5fopen_f(file_path,H5F_ACC_RDWR_F,file_id,error) !Open file

        ! Save RHO profile only once
        if (first_call) then
            !Open RHO dataset
            call h5dopen_f(file_id, 'RHO', dset_id, error)
            call h5dget_space_f(dset_id, dspace_id, error)
            !Select subset
            call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, int([0],8), int([NZ],8), error)
            !Create memory data space
            call h5screate_simple_f(1, int([NZ],8), memspace, error)
            !Write subset to dataset
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, RHO, int([NZ],8), error, memspace, dspace_id)
            !Close dataset
            call h5dclose_f(dset_id, error)
            first_call = .false.
        end if

        !======================================= NX by NZ ==================================================
        
        count = int([NX, NZ, 1],8)
        offset = [0,0,IT/PERIOD]
        dimsm = count
        !Create memory data space
        call h5screate_simple_f(rnk, dimsm, memspace, error)

        !Open THETA dataset
        call h5dopen_f(file_id, 'THETA', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, THETA, dimsm, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        !Open RV dataset
        call h5dopen_f(file_id, 'RV', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, RV, dimsm, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        !Open RL dataset
        call h5dopen_f(file_id, 'RL', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, RL, dimsm, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        !Open RI dataset
        call h5dopen_f(file_id, 'RI', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, RI, dimsm, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        !Open DPB dataset
        call h5dopen_f(file_id, 'DPB', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, DPB, dimsm, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        !Open SAT dataset
        call h5dopen_f(file_id, 'SAT', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, SAT, dimsm, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        !====================================== NXP by NZP =================================================
        count = int([NXP, NZP, 1],8)
        dimsm = count
        !Create memory data space
        call h5screate_simple_f(rnk, dimsm, memspace, error)

        !Open UXN dataset
        call h5dopen_f(file_id, 'UXN', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        call h5sget_simple_extent_dims_f(dspace_id, current_dims, max_dims, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, UXN, dimsm, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        !Open UZN dataset
        call h5dopen_f(file_id, 'UZN', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        call h5sget_simple_extent_dims_f(dspace_id, current_dims, max_dims, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, UZN, dimsm, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)
        
        !====================================== N_sd by 2 ===================================================

        count = int([N_sd, 2, 1],8)
        dimsm = count
        !Create memory data space
        call h5screate_simple_f(rnk, dimsm, memspace, error)

        !Open X_SD dataset
        call h5dopen_f(file_id, 'X_SD', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, X, dimsm, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        !====================================== N_sd by 1 ===================================================
        count2 = int([N_sd, 1],8)
        offset2 = [0,IT/PERIOD]
        dimsm2 = count2
        !Create memory data space
        call h5screate_simple_f(2, dimsm2, memspace, error)

        !Open R_SD dataset
        call h5dopen_f(file_id, 'R_SD', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset2, count2, error, stride2, block2)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, R(2,:), dimsm2, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        !====================================== NT by 1 ===================================================

        !Create memory data space
        call h5screate_simple_f(1, int([1],8), memspace, error)

        !Open TIME dataset
        call h5dopen_f(file_id, 'TIME', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, int([IT/PERIOD],8), int([1],8), error, int([1],8), int([1],8))
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time, int([1],8), error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)
        !Close everything opened
        call h5sclose_f(dspace_id, error)
        call h5sclose_f(memspace , error)
        call h5fclose_f(file_id  , error)
        call h5close_f(error)
    END IF
END SUBROUTINE HDF5_SAVE_RESULTS

SUBROUTINE HDF5_READ_VELOCITIES(filename)
    use HDF5_VARIABLES
    use STEPPER
    use GRID
    use ADVECTION, ONLY: UXT, UZT

    implicit NONE

    character(100), intent(in) :: filename

    offset = [0,0,NINT(TIME/DT_ADV)]

    call h5open_f(error) !Open fortran HDF5 interface
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error) !Open file
    if (error /= 0) then
        call system('clear')
        stop 'HDF5_READ_VELOCITIES: Velocities file not found.'
    end if
    !Open UXT dataset
    call h5dopen_f(file_id, 'UXT', dset_id, error)
    !Get dataspace info
    call h5dget_space_f(dset_id, dspace_id, error)
    !Select subset
    count = int([NXP, NZ, 1],8)
    dimsm = count
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
    !Create memory data space
    call h5screate_simple_f(rnk, dimsm, memspace, error)
    !Read subset from dataset
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, UXT, dimsm, error, memspace, dspace_id)
    !Close dataset
    call h5dclose_f(dset_id, error)
    !Close dataspace
    call h5sclose_f(dspace_id, error)
    !Close memory dataspace
    call h5sclose_f(memspace , error)

    !Open UZT dataset
    call h5dopen_f(file_id, 'UZT', dset_id, error)
    call h5dget_space_f(dset_id, dspace_id, error)
    !Select subset
    count = int([NXP, NZ, 1],8)
    dimsm = count
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
    !Create memory data space
    call h5screate_simple_f(rnk, dimsm, memspace, error)
    !Write subset to dataset
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, UZT, dimsm, error, memspace, dspace_id)
    !Close dataset
    call h5dclose_f(dset_id, error)
    !Close dataspace
    call h5sclose_f(dspace_id, error)
    !Close memory dataspace
    call h5sclose_f(memspace , error)

    !Close file and interface
    call h5fclose_f(file_id, error)
    call h5close_f(error)
END SUBROUTINE HDF5_READ_VELOCITIES

SUBROUTINE HDF5_READ_PSI(filename)
   use HDF5_VARIABLES
   use STEPPER
   use GRID
   use ADVECTION, ONLY: PSI

   implicit NONE

   character(100), intent(in) :: filename
   integer, save :: offset_max 
   integer, save :: call_count = 0

   if (call_count==0) offset_max = nint(T_MAX/DT_ADV) + 1
   offset     = [0, 0, min(call_count, offset_max)]
   call_count = call_count + 1

   call h5open_f(error) !Open fortran HDF5 interface
   call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error) !Open file
   if (error /= 0) then
       call system('clear')
       stop 'HDF5_READ_PSI: file not found.'
   end if
   !Open PSI dataset
   call h5dopen_f(file_id, 'PSI', dset_id, error)
   !Get dataspace info
   call h5dget_space_f(dset_id, dspace_id, error)
   !Select subset
   count = int([NXP, NZP, 1],8)
   dimsm = count
   call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
   !Create memory data space
   call h5screate_simple_f(rnk, dimsm, memspace, error)
   !Read subset from dataset
   call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, PSI, dimsm, error, memspace, dspace_id)
   !Close dataset
   call h5dclose_f(dset_id, error)
   !Close dataspace
   call h5sclose_f(dspace_id, error)
   !Close memory dataspace
   call h5sclose_f(memspace , error)

   !Close file and interface
   call h5fclose_f(file_id, error)
   call h5close_f(error)
END SUBROUTINE HDF5_READ_PSI

SUBROUTINE GET_ENVIRONMENT_PROFILES_DYCOMS_II

   !Follows Stevens et al. (2005) with an inversion at z_inv

   use ENVIRONMENT
   use GRID
   use CONSTANTS
   use FUNCTIONS, only: linspace, interpol, SVP, EXNER
   use ADVECTION, only: eps, eps_input
   IMPLICIT NONE

   REAL*8  :: EE
   REAL*8  :: S, P
   real*8, allocatable :: z_table(:),theta_table(:),rv_table(:)
   real*8, allocatable :: Z(:)
   integer :: K
    
   !Mixed-phase profiles
   allocate(z_table(6), theta_table(6), rv_table(6))
   z_table     = [  0.0D0, 500.00D0, 900.0D0, 1100.0D0, 1150.00D0, 1800.D0]
   theta_table = [264.0D0, 265.50D0, 266.0D0,  266.5D0,  271.00D0,  276.D0]
   rv_table    = [  1.8D0,   1.45D0,   1.3D0,    1.2D0,    1.65D0,    1.D0]/1000

   Z = linspace(0.5*DZ, LZ-0.5*DZ, NZ)

   ALLOCATE ( TH_E(NZ), RV_E(NZ) )
   TH_E = interpol(z_table,theta_table,Z)
   RV_E = interpol(z_table,rv_table,Z)

   CALL GET_ENV_PRESSURE_TEMP

   !Base state density profile
   ALLOCATE (RHO(NZ))

   OPEN (UNIT = 99, FILE = "environ_profile.out")
   DO K = 1,NZ
      EE = PR_E(K)*( RV_E(K)/( RV_E(K) + R_D/R_V ) ) !Partial vapor pressure
      S  = EE/SVP(T_E(K),'L') - 1.D0 !Supersaturation

      RHO(K) = (PR_E(K)-EE)/(R_D*T_E(K))

      WRITE(99,*) Z(K), TH_E(K), RV_E(K)*1.D3, PR_E(K), T_E(K), RHO(K), S, eps(k), 0.D0
   END DO
   CLOSE (99, STATUS = "KEEP")

   !Gathering value of Pressure and Exner function at grid box center
   !from environmental pressure profile
   if (.not.allocated(Exn)) allocate(Exn(NZ))
   do k = 1,NZ
      P      = interpol(Z_E,PR_E,(k - 0.5D0)*DZ)
      Exn(k) = EXNER(P)
   end do

   ! Building the eddy diffusivity profile
   allocate(EDDY_DIFF(NZ))
   TT        = DZ !Transition thickness across inversion level
   EDDY_MAX   = (1.D2**(4.D0/3))*eps_input**(1.D0/3) ! Melhorar: Importar DELTA e EPS do gerador de velocidades.
   EDDY_DIFF = 0.5D0*EDDY_MAX*ERFC((Z - Z_INV)/(TT*sqrt(2.D0)))
   
END SUBROUTINE GET_ENVIRONMENT_PROFILES_DYCOMS_II

SUBROUTINE GET_ENV_PRESSURE_TEMP

   use ENVIRONMENT
   use GRID
   use CONSTANTS
   use FUNCTIONS, only: linspace

   implicit none
   integer :: K

   REAL*8 :: KAPPA
   REAL*8 :: THV,THV_P,THV_M

   KAPPA = R_D/CP_D

   ALLOCATE ( PR_E(NZ) )
   ALLOCATE (  T_E(NZ) )

   Z_E = linspace(0.5*DZ, LZ-0.5*DZ, NZ)

   PR_E(1) = 1.015D5

   DO K = 2,NZ
      THV_M = TH_E(K-1)*( 1 + 0.608D0*RV_E(K-1) )
      THV_P = TH_E( K )*( 1 + 0.608D0*RV_E(K)   )
      THV   = 0.5D0*( THV_M + THV_P )

      PR_E(K) = ( PR_E(K-1)**KAPPA - ( (G*P00**KAPPA)/(CP_D*THV) )*(Z_E(K)-Z_E(K-1)) )**(1.D0/KAPPA)
   END DO

   T_E = TH_E*(PR_E/P00)**KAPPA

END SUBROUTINE GET_ENV_PRESSURE_TEMP