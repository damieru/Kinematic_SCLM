PROGRAM MAIN
    use STEPPER
    use THERMODYNAMIC_VAR
    !use ADVECTION
    use ADV_FUNCTIONS
    use HDF5_VARIABLES
    use SD_FUNCTIONS
    use BOX
    use FUNCTIONS, only: progress_bar, exner
    use OMP_LIB

    IMPLICIT NONE

    INTEGER*8 :: I

    !Initialization -----------------------------------------------------
    CALL PRINT_REMARKS
    CALL READIN                              !Reads the input parameters
    CALL GET_ENVIRONMENT_PROFILES_DYCOMS_II  !Sets up the environmental base state with a strong temperature inversion at z_inv
    CALL INIT_THERMODYNAMIC_VAR_SD           !Initializes the thermodynamic fields using environmental profiles
    CALL INIT_RANDOM_SEED                    !For the (random) prescribed fluid flow
    CALL INIT_SD
    CALL UPDATE_BOXES
    CALL UPDATE_PARTICLE_BOX_MAP
    CALL SAT_FIELD

    TIME = 0.D0
    T_MAX = 60 !Maximum simulation time in seconds
    OUT_PER = 1.D0/60 !Output period in minutes
    N_STEPS = NINT(T_MAX/DT) + 1

    CALL GET_VELOCITIES
    CALL HDF5_CREATE_FILE
    CALL HDF5_SAVE_RESULTS(OUT_PER)
    TIME = TIME + DT
    DO I = 2,N_STEPS
        CALL CHECK_ADV
        IF (ADVECT) THEN
            !Updates thermodynamic fields using Euler scheme ------------------
            CALL ADVECTION_MPDATA
            IF (FLAG_TRANSIENT_FLOW.EQ.1) THEN
                CALL GET_VELOCITIES
            END IF
        END IF
        CALL ADVECTION_SD
        CALL UPDATE_PARTICLE_BOX_MAP
        CALL GROW_DROPLETS      ! Solve growth equations
        CALL UPDATE_BOXES       ! Update grid boxes mean properties after condensation
        CALL SAT_FIELD          ! Diagnostic only
        
        call HDF5_SAVE_RESULTS(OUT_PER)

        if (mod(int(TIME/DT),10) == 0) then
            call progress_bar(TIME/T_MAX)
        end if
        TIME = TIME + DT
    END DO
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
   use PERIOD3D
   use VELOCITY_MODES, ONLY: MK,NK,LLX,LLZ,DZ_INV_W
   use SD_VARIABLES,   ONLY: im_freezing, init_ice_frac
   use ADVECTION
   use STEPPER
   use ENVIRONMENT, ONLY: Z_INV
   use IO_PARAMETERS
   use SD_VARIABLES, only: DPB_input

   IMPLICIT NONE   

   character(len=100) :: file_line

   ! GRID ------------------------------------------------------------------------
   file_line = get_line(); READ(file_line,*) NX, NZ
   file_line = get_line(); READ(file_line,*) DX, DZ

   NXP = NX + 1
   NZP = NZ + 1

   !Period 3-D
   LX = DX*NX
   LZ = DZ*NZ

   ! VELOCITY MODES --------------------------------------------------------------
   file_line = get_line(); READ(file_line,*) MK, NK
   file_line = get_line(); READ(file_line,*) FLAG_TRANSIENT_FLOW

   LLX = LX
   LLZ = 2.D0*LZ

   WRITE(*,'(/A,F6.2)') "Turbulence resolution length scale: ", SQRT( LLX*LLZ/REAL(MK*NK) )

   ! STEPPER ---------------------------------------------------------------------
   file_line = get_line(); READ(file_line,*) DT
   file_line = get_line(); READ(file_line,*) DT_ADV

   ! INVERSION -------------------------------------------------------------------
   file_line = get_line(); READ(file_line,*) Z_INV    !Thermodynamic initial inversion
   file_line = get_line(); READ(file_line,*) DZ_INV_W  !Velocity profile overshoot

   ! OUTPUT PARAMETERS -----------------------------------------------------------
   file_line = get_line(); READ(file_line,*) output_folder
   file_line = get_line(); READ(file_line,*) output_file

   ! Number of super droplets per box
   file_line = get_line(); READ(file_line,*) DPB_input

   !Eddy size and TKE dissipation rate--------------------------------------------
   file_line = get_line(); READ(file_line,*) L
   file_line = get_line(); READ(file_line,*) eps

   ! Uses the deterministic model if either L or eps = 0
   if (eps*L == 0.D0) then
      L   = 0.D0
      eps = 0.D0
   end if

   ! Immersion freezing ON / OFF
   file_line = get_line(); READ(file_line,*) im_freezing

   ! Initial ice fraction
   file_line = get_line(); READ(file_line,*) init_ice_frac

   ! Path to velocities file
   file_line = get_line(); READ(file_line,*) field
   file_line = get_line(); READ(file_line,*) input_file

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

SUBROUTINE GET_STREAM_FUNCTION

  !Gets the streamfunction "PSI"

  !Time dependence is incorporated into the coefficients
  !{a_mn,b_mn}. Their time evolution is controlled by the subroutine
  !UPDATE_AMN_BMN

  use ADVECTION
  use GRID

  IMPLICIT NONE

  REAL*8  :: X,Z
  INTEGER*8 :: I,K

  INTEGER :: INIT = 0
  SAVE INIT

  IF ( INIT.EQ.0 ) THEN
     ALLOCATE ( PSI(NXP,NZP) )
     CALL INIT_VELOCITY_MODES !Initializes the modes stuff
     INIT = -1
  END IF

  DO I = 1,NXP
     DO K = 1,NZP
        X = DX*(I-1)
        Z = DZ*(K-1)
        CALL STREAM_FUNCTION_PINSKY ( PSI(I,K), X, Z )
     END DO
  END DO

END SUBROUTINE GET_STREAM_FUNCTION

SUBROUTINE INIT_THERMODYNAMIC_VAR_SD

  use GRID
  use ENVIRONMENT
  use THERMODYNAMIC_VAR

  IMPLICIT NONE

  REAL*8  :: RHO_BAR
  INTEGER*8 :: IZ

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

  !-----------------------------------------------------------------


  !******************** A T T E N T I O N **************************
  !*****************************************************************
  !For Boussinesq approximation as in Pinsky et al.
  !Sets a constant value for dry air density distribution

  RHO_BAR = 0.D0
  DO IZ = 1,NZ
     RHO_BAR = RHO_BAR + RHO(IZ)
  END DO
  RHO_BAR = RHO_BAR/REAL(NZ)

  !Or simply: RHO_BAR = RHO(1)

  RHO = RHO_BAR

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

  SEED = CLOCK + 37 * (/ (I - 1, I = 1, N) /)

  CALL RANDOM_SEED(PUT = SEED)

  DEALLOCATE(SEED)

END SUBROUTINE INIT_RANDOM_SEED

SUBROUTINE HDF5_CREATE_FILE

   use HDF5_VARIABLES
   use GRID
   use STEPPER, ONLY: T_MAX, OUT_PER
   use SD_VARIABLES, only: N_sd
   use IO_PARAMETERS
   
   integer*8 :: NT
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

   NT = int(NINT(T_MAX)/(60*OUT_PER) + 1,8) ! Number of saved instants
   call h5open_f(error) !Opens HDF5 Fortran interface
   call h5fcreate_f(file_path, H5F_ACC_TRUNC_F, file_id, error) !Creates file

   call h5screate_simple_f(1, (/NZ/), dspace_id, error) !Creates dataspace
   !Creating RHO dataset
   call h5dcreate_f(file_id, 'RHO', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)

   call h5screate_simple_f(rnk, (/NX,NZ,NT/), dspace_id, error) !Creates dataspace
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
   !Creating UX dataspace and dataset
   call h5dcreate_f(file_id, 'UX', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)
   !Creating UZ dataspace and dataset
   call h5dcreate_f(file_id, 'UZ', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)
   !Creating DPB dataspace and dataset
   call h5dcreate_f(file_id, 'DPB', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)
   !Creating SAT dataspace and dataset
   call h5dcreate_f(file_id, 'SAT', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)

   call h5screate_simple_f(rnk, (/int(N_sd,8),int(2,8),NT/), dspace_id, error) !Creates dataspace
   !Creating X_SD dataspace and dataset
   call h5dcreate_f(file_id, 'X_SD', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)

   call h5screate_simple_f(2, (/int(N_sd,8),NT/), dspace_id, error) !Creates dataspace
   !Creating X_SD dataspace and dataset
   call h5dcreate_f(file_id, 'R_SD', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)

   call h5screate_simple_f(1, (/NT/), dspace_id, error) !Creates dataspace
   !Creating X_SD dataspace and dataset
   call h5dcreate_f(file_id, 'TIME', H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
   call h5dclose_f(dset_id,error)

   !Closing file and interface
   call h5fclose_f(file_id,error)
   call h5close_f(error)

END SUBROUTINE HDF5_CREATE_FILE

SUBROUTINE HDF5_SAVE_RESULTS (NMIN)

    !Prints the output every "NMIN" minutes
    use HDF5_VARIABLES
    use STEPPER
    use ADVECTION, ONLY: UX,UZ
    use THERMODYNAMIC_VAR, ONLY: THETA, RV, RL, RI, DPB, SAT
    use ENVIRONMENT, ONLY: RHO
    use GRID
    use SD_VARIABLES, ONLY: N_sd, X, R
    use IO_PARAMETERS
    IMPLICIT NONE

    REAL*8  :: NMIN
    INTEGER :: IT,FREQ
    CHARACTER(100) :: FILE_PATH

    FILE_PATH = TRIM(output_folder)//'/'//TRIM(output_file)

    IT = NINT( TIME/DT )
    FREQ = NINT( NMIN * 60.D0/DT )

    IF ( MOD(IT,FREQ).EQ.0 ) THEN

        count = (/NX,NZ,int(1,8)/)
        offset = (/0,0,IT/FREQ/)
        dimsm = count

        call h5open_f(error) !Open fortran HDF5 interface
        call h5fopen_f(file_path,H5F_ACC_RDWR_F,file_id,error) !Open file
        if (NINT(TIME/60.D0) == 0) then
            !Open RHO dataset
            call h5dopen_f(file_id, 'RHO', dset_id, error)
            call h5dget_space_f(dset_id, dspace_id, error)
            !Select subset
            call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, (/int(0,8)/), (/NZ/), error)
            !Create memory data space
            call h5screate_simple_f(1, (/NZ/), memspace, error)
            !Write subset to dataset
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, RHO, dimsm, error, memspace, dspace_id)
            !Close dataset
            call h5dclose_f(dset_id, error)
        end if

        !Open THETA dataset
        call h5dopen_f(file_id, 'THETA', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Create memory data space
        call h5screate_simple_f(rnk, dimsm, memspace, error)
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

        !Open UX dataset
        call h5dopen_f(file_id, 'UX', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, UX, dimsm, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        !Open UZ dataset
        call h5dopen_f(file_id, 'UZ', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, UZ, dimsm, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        !Open DPB dataset
        call h5dopen_f(file_id, 'DPB', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, DPB, dimsm, error, memspace, dspace_id)
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
        
        count = (/int(N_sd,8),int(2,8),int(1,8)/)
        offset = (/0,0,IT/FREQ/)
        dimsm = count
        !Open X_SD dataset
        call h5dopen_f(file_id, 'X_SD', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, count, error, stride, block)
        !Create memory data space
        call h5screate_simple_f(rnk, dimsm, memspace, error)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, X, dimsm, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        count2 = (/int(N_sd,8),int(1,8)/)
        offset2 = (/0,IT/FREQ/)
        dimsm2 = count2
        !Open R_SD dataset
        call h5dopen_f(file_id, 'R_SD', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset2, count2, error, stride2, block2)
        !Create memory data space
        call h5screate_simple_f(2, dimsm2, memspace, error)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, R(2,:), dimsm2, error, memspace, dspace_id)
        !Close dataset
        call h5dclose_f(dset_id, error)

        !Open TIME dataset
        call h5dopen_f(file_id, 'TIME', dset_id, error)
        call h5dget_space_f(dset_id, dspace_id, error)
        !Select subset
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, int((/IT/FREQ/),8), int((/1/),8), error, int((/1/),8), int((/1/),8))
        !Create memory data space
        call h5screate_simple_f(1, int((/1/),8), memspace, error)
        !Write subset to dataset
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time, int((/1/),8), error, memspace, dspace_id)
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

    offset = (/0,0,NINT(TIME/DT_ADV)/)

    call h5open_f(error) !Open fortran HDF5 interface
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error) !Open file
    if (error /= 0) then
        call system('clear')
        print*, 'Velocities file not found.'
        stop
    end if
    !Open UXT dataset
    call h5dopen_f(file_id, 'UXT', dset_id, error)
    !Get dataspace info
    call h5dget_space_f(dset_id, dspace_id, error)
    !Select subset
    count = (/NXP,NZ,int(1,8)/)
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
    count = (/NX,NZP,int(1,8)/)
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

   offset = (/0,0,NINT(TIME/DT_ADV)/)

   call h5open_f(error) !Open fortran HDF5 interface
   call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error) !Open file
   if (error /= 0) then
       call system('clear')
       print*, 'PSI file not found.'
       stop
   end if
   !Open PSI dataset
   call h5dopen_f(file_id, 'PSI', dset_id, error)
   !Get dataspace info
   call h5dget_space_f(dset_id, dspace_id, error)
   !Select subset
   count = (/NXP,NZP,int(1,8)/)
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

SUBROUTINE INIT_VELOCITY_MODES

  use VELOCITY_MODES
  use FUNCTIONS, ONLY: BOX_MULLER
  use PERIOD3D
  use GRID      !In case the grid resolution dictates the number of Fourier modes

  IMPLICIT NONE

  REAL*8  :: PI
  INTEGER :: N,M

  REAL*8  :: SUM
  REAL*8  :: Z,ZMIN,ZMAX,DDZ
  REAL*8  :: X,XMIN,XMAX,DDX
  INTEGER :: NNZ,IZ
  INTEGER :: NNX,IX

  REAL*8, ALLOCATABLE :: G(:),H(:)
  REAL*8  :: G1,G2,H1,H2

  REAL*8  :: F!,FMAX

  REAL*8  :: FZ,DFZ

  !*****************************************************************************
  !*****************************************************************************

  ! * * * THESE ARE NOW SET AT "READIN" FROM PARAMETERS IN THE INPUT FILE * * *

  !Periods ---------------------------------------------------------------------

  !LLX = LX
  !LLZ = 2.D0*LZ

  !Parameters ------------------------------------------------------------------

  !MK = 25 !NINT( LLX/DX ) !50
  !NK = 17 !NINT( LLZ/DZ ) !50

  !******************************************************************************
  !******************************************************************************

  WRITE(*,*) "MK = ", MK
  WRITE(*,*) "NK = ", NK
  WRITE(*,*)
  WRITE(*,*) "T_MIN = ", (1.D-3)**(-1.D0/3.D0)*(LLZ/REAL(NK))**(2.D0/3.D0)
  WRITE(*,*)

  !*****************************************************************************

  ALLOCATE ( KM(MK), KN(NK) )
  ALLOCATE ( DM(MK), CN(NK) )

  ALLOCATE ( AMN(MK,NK) )
  ALLOCATE ( BMN(MK,NK) )

  !Wavenumbers -----------------------------------------------------------------

  PI = ACOS(-1.D0)

  DO M = 1,MK
     KM(M) = (2.D0*PI/LLX)*M !x-direction
  END DO

  DO N = 1,NK
     KN(N) = (2.D0*PI/LLZ)*N !z-direction
  END DO

  !Random "unsteadiness" coefficients - Initialization with
  !independent unit Gaussians

  DO M = 1,MK
     DO N = 1,NK
        CALL BOX_MULLER (G1,G2)
        AMN(M,N) = G1
        BMN(M,N) = G2
     END DO
  END DO

  ! * * *

  !Fourier coefficients --------------------------------------------------------

  !{C_n} from the vertical profile of vertical velocity variance

  ZMIN = - 0.5D0*LLZ
  ZMAX =   0.5D0*LLZ
  NNZ  = 2**10
  DDZ  = (ZMAX-ZMIN)/REAL(NNZ)

  !We assume a simple vertical profile of vertical velocity variance

  ALLOCATE ( G(0:NNZ) )

  open (unit = 59, file = "f.dat")
  DO IZ = NNZ/2,NNZ
     Z = ZMIN + DDZ*IZ
     CALL VELOCITY_VARIANCE_PROFILE ( F, Z )
     write(59,*) z, f
     G(IZ) = SQRT( F )
  END DO
  close (59,status = "keep")

  !Antisymmetrical reflection of "G"

  DO IZ = 0,NNZ/2-1
     G(IZ) = - G(NNZ-IZ)
  END DO

  open (unit = 59, file = "g.dat")
  do iz = 0,nnz
     z = zmin + ddz*iz
     write(59,*) z, g(iz)
  end do
  close (59,status = "keep")

  CN = 0.D0
  DO N = 1,NK
     SUM = 0.D0
     DO IZ = 0,NNZ-1
        Z = ZMIN + DDZ*IZ
        G1 = SIN(KN(N)* Z)     *G(IZ  )
        G2 = SIN(KN(N)*(Z+DDZ))*G(IZ+1)
        SUM = SUM + 0.5D0*(G1+G2)*DDZ
     END DO
     CN(N) = 2.D0*SUM/LLZ
  END DO

  !-----------------------------------------------------------------------------

  !{D_m} from the (vertical-velocity) transverse correlation function

  XMIN = - 0.5D0*LLX !0.D0
  XMAX =   0.5D0*LLX !LLX
  NNX  = 2**10
  DDX  = (XMAX-XMIN)/REAL(NNX)

  !We assume a simple (normlized) transverse correlation function for
  !the vertical velocity in the for of an exponential decay with
  !integral length scale L = 200m

  ALLOCATE ( H(0:NNX) )

  DO IX = NNX/2,NNX
     X = XMIN + DDX*IX
     IF (X.LE.250.D0) THEN
        H(IX) = ( COS(0.5D0*PI*X/250.D0) )**2
     ELSE
        H(IX) = 0.D0
     END IF
  END DO
  !"Symmetrizing" --------------------
  DO IX = 0,NNX/2-1
     H(IX) = H(NNX-IX)
  END DO
  !-----------------------------------

  open (unit = 59, file = "h.dat")
  do ix = 0,nnx
     x = xmin + ddx*ix
     write(59,*) x, h(ix)
  end do
  close (59,status = "keep")

  DM = 0.D0
  DO M = 1,MK
     SUM = 0.D0
     DO IX = 0,NNX-1
        X = XMIN + DDX*IX
        H1 = COS(KM(M)* X)     *H(IX  )
        H2 = COS(KM(M)*(X+DDX))*H(IX+1)
        SUM = SUM + 0.5D0*(H1+H2)*DDX
     END DO
    !DM(M) = 2.D0*SUM/LLX
     DM(M) = SUM/LLX !Like in Pinsky
  END DO

  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !SOME TESTS...
  !-----------------------------------------------------------------------------

  !Check normalization of the coefficients {D_m}

  SUM = 0.D0
  DO M = 1,MK
     SUM = SUM + DM(M)**2
  END DO

  DM = DM/SQRT(SUM)

  SUM = 0.D0
  DO M = 1,MK
     SUM = SUM + DM(M)**2
  END DO

  WRITE(*,*) "NORMALIZATION D_M: ", SUM

  !-----------------------------------------------------------------------------
  !Prints the vertical profile of the vertical velocity variance

  OPEN (UNIT = 39, FILE = "w_variance.dat")

  DO IZ = NNZ/2+1,NNZ-1

     Z = ZMIN + DDZ*IZ

     CALL VERTICAL_FUNCTION (FZ,DFZ,Z)

     SUM = 0.D0
     DO N = 1,NK
        SUM = SUM + ( CN(N)*SIN(KN(N)*Z) )**2
     END DO

     WRITE(39,*) Z, FZ*FZ*SUM

  END DO

  CLOSE (39, STATUS = "KEEP")

  !-----------------------------------------------------------------------------
  !Prints the transverse vertical velocity correlation function

  OPEN (UNIT = 39, FILE = "cw.dat")

  DO IX = 0,NNX
     X = XMIN + DDX*IX
     SUM = 0.D0
     DO M = 1,MK
        SUM = SUM + DM(M)**2*COS(KM(M)*X)
     END DO
     WRITE(39,*) X, SUM, H(IX)
  END DO

  CLOSE (39, STATUS = "KEEP")

END SUBROUTINE INIT_VELOCITY_MODES

SUBROUTINE UPDATE_AMN_BMN ( DT )

  !Updates the random "unsteadiness" coefficients

  use VELOCITY_MODES
  use FUNCTIONS, ONLY: BOX_MULLER

  IMPLICIT NONE

  REAL*8  :: DT
  REAL*8  :: G1,G2
  REAL*8  :: TAU
  REAL*8  :: PI
  REAL*8  :: EPS = 1.D-3 !Local dissipation rate

  INTEGER :: M,N

  PI = ACOS(-1.D0)

  DO N = 1,NK

    !Model based on dimensional analysis----------------------
     TAU = EPS**(-1.D0/3.D0) * ( 2.D0*PI/KN(N) )**(2.D0/3.D0)
    !---------------------------------------------------------

     DO M = 1,MK
        CALL BOX_MULLER (G1,G2)
        AMN(M,N) = AMN(M,N)*EXP(-DT/TAU) + G1*SQRT(1.D0-EXP(-2.D0/TAU))
        BMN(M,N) = BMN(M,N)*EXP(-DT/TAU) + G2*SQRT(1.D0-EXP(-2.D0/TAU))
     END DO
  END DO

END SUBROUTINE UPDATE_AMN_BMN

SUBROUTINE VELOCITY_FIELD_PINSKI_ET_AL(U,W,X,Z)

  !Incompressible synthetic flow following Pinsky et al. (JAS, 2008).

  use VELOCITY_MODES

  IMPLICIT NONE

  REAL*8 :: U,W,X,Z
  REAL*8 :: UX,WX
  REAL*8 :: F,DF
  REAL*8 :: SHEAR_RATE
  INTEGER :: N,M

  CALL VERTICAL_FUNCTION (F,DF,Z)

  !Vertical component -----------------------------------------------------------

  W = 0.D0
  DO N = 1,NK

     WX = 0.D0
     DO M = 1,MK
        WX = WX + DM(M)*( AMN(M,N)*SIN( KM(M)*X ) + BMN(M,N)*COS( KM(M)*X ) )
     END DO

     W = W + CN(N)*WX*SIN( KN(N)*Z )

  END DO

  W = F*W

  !Horizontal component ---------------------------------------------------------

  SHEAR_RATE = 3.D0/850.D0

  U = 0.D0

  DO N = 1,NK

     UX = 0.D0
     DO M = 1,MK
        UX = UX + ( DM(M)/KM(M) )*( AMN(M,N)*COS(KM(M)*X) - BMN(M,N)*SIN(KM(M)*X) )
     END DO

     U = U + CN(N)*UX*( DF*SIN(KN(N)*Z) + F*KN(N)*COS(KN(N)*Z) )

  END DO

  U = U + SHEAR_RATE*Z


END SUBROUTINE VELOCITY_FIELD_PINSKI_ET_AL

SUBROUTINE STREAM_FUNCTION_PINSKY (PSI,X,Z)

  !Generates the streamfunction that produces the synthetic flow
  !proposed in Pinsky et al. (JAS, 2008).

  use VELOCITY_MODES

  IMPLICIT NONE

  REAL*8 :: PSI,X,Z
  REAL*8 :: F,DF,SX
  REAL*8 :: SHEAR_RATE
  INTEGER :: N,M

  SHEAR_RATE = 0.D0 !1.D0/850.D0

  CALL VERTICAL_FUNCTION (F,DF,Z) !Only "F" will be used here

  PSI = 0.D0
  DO N = 1,NK
     SX = 0.D0
     DO M = 1,MK
        SX = SX + DM(M)*( BMN(M,N)*SIN(KM(M)*X) - AMN(M,N)*COS(KM(M)*X) )/KM(M)
     END DO
     PSI = PSI + CN(N)*SX*SIN( KN(N)*Z )
  END DO

  PSI = F*PSI - 0.5D0*SHEAR_RATE*Z**2

END SUBROUTINE STREAM_FUNCTION_PINSKY

SUBROUTINE VERTICAL_FUNCTION (F,DF,Z)

  use VELOCITY_MODES

  IMPLICIT NONE

  REAL*8 :: F,DF,Z
  REAL*8 :: F1,F2,F22,F3
  REAL*8 :: DF1,DF2

  INTEGER :: N

  IF ( (Z.EQ.0.D0).OR.(Z.EQ.(0.5D0*LLZ)) ) THEN

     F  = 0.D0
     DF = 0.D0

  ELSE

     F1 = 0.D0
     DO N = 1,NK
        F1 = F1 + CN(N)*SIN( KN(N)*Z )
     END DO

     F22 = 0.D0
     DO N = 1,NK
        F22 = F22 + ( CN(N)*SIN( KN(N)*Z ) )**2
     END DO

     F2 = SQRT(F22)

     !--------------------------------------------------------------------

     F = F1/F2

     !--------------------------------------------------------------------

     F3 = 0.D0
     DO N = 1,NK
        F3 = F3 + 2.D0*KN(N)*CN(N)*CN(N)*SIN( KN(N)*Z )*COS( KN(N)*Z )
     END DO

     DF1 = 0.D0
     DO N = 1,NK
        DF1 = DF1 + KN(N)*CN(N)*COS( KN(N)*Z )
     END DO

     DF2 = 0.5D0*F3/F2

     !--------------------------------------------------------------------

     DF = DF1/F2 - F1*DF2/F22

     !--------------------------------------------------------------------

  END IF

END SUBROUTINE VERTICAL_FUNCTION

SUBROUTINE VELOCITY_VARIANCE_PROFILE (W2,Z)

  !Approximates experimental or simulation data for "<w^2>(z)" using a
  !Lagrange interpolating polynomial

  use ENVIRONMENT, ONLY: Z_INV
  use VELOCITY_MODES, ONLY: DZ_INV_W
  
  IMPLICIT NONE

  REAL*8  :: W2,Z
  REAL*8, ALLOCATABLE :: X(:),Y(:)
  REAL*8 :: P,PJ
  INTEGER :: N
  INTEGER :: J,K

  N = 5

  ALLOCATE (X(N),Y(N))

  X(1) = 0.D0
  Y(1) = 0.D0

  X(2) = 200.D0
  Y(2) = 0.34D0

  X(3) = 400.D0
  Y(3) = 0.48D0

  X(4) = 600.D0
  Y(4) = 0.5D0

  X(5) = Z_INV + DZ_INV_W !850.D0
  Y(5) = 0.D0

  ! * * *

  P = 0.D0
  DO J = 1,N
     PJ = 1.D0
     DO K = 1,N
        IF (K.NE.J) THEN
           PJ = PJ*( Z - X(K) )/( X(J) - X(K) )
        END IF
     END DO
     PJ = Y(J)*PJ
     P = P + PJ
  END DO

  ! * * *

  W2 = (0.5D0/0.514D0)*P !"Gambiarra"

  !Inversion layer *********************

  IF ( Z.GE.(Z_INV + DZ_INV_W) ) THEN
     W2 = 0.D0
  END IF

  !*************************************

END SUBROUTINE VELOCITY_VARIANCE_PROFILE

FUNCTION E_VAP(RV,P)

  !Partial vapor pressure

  use CONSTANTS

  IMPLICIT NONE

  REAL*8 :: E_VAP
  REAL*8 :: RV,P
  REAL*8 :: EPS

  EPS = R_D/R_V

  E_VAP = P*( RV/(RV + EPS) )

END FUNCTION E_VAP

FUNCTION R_MOIST (Q)

  !Gas constant of moist air

  use CONSTANTS

  IMPLICIT NONE

  REAL*8 :: R_MOIST
  REAL*8 :: Q !VAPOR MIXING RATIO

  R_MOIST = (R_D + Q*R_V)/(1.D0 + Q)

END FUNCTION R_MOIST

FUNCTION CP_MOIST(Q)

  !Specific heat of moist air

  use CONSTANTS

  IMPLICIT NONE

  REAL*8 :: CP_MOIST
  REAL*8 :: Q !VAPOR MIXING RATIO

  CP_MOIST = (CP_D + Q*CP_V)/(1.D0 + Q)

END FUNCTION CP_MOIST

FUNCTION ES_LIQ(T)

  !SATURATION VAPOR PRESSURE OVER LIQUID WATER
  !ASSUME CONSTANT LATENT HEAT "LV = LV0"

  use CONSTANTS

  IMPLICIT NONE

  REAL*8 :: ES_LIQ
  REAL*8 :: T

  ES_LIQ = ES0*EXP( (LV0/R_V)*(1.D0/T0 - 1.D0/T) )

END FUNCTION ES_LIQ

FUNCTION ES_ICE(T)

  !SATURATION VAPOR PRESSURE OVER ICE
  !ASSUME CONSTANT LATENT HEAT "LS = LS0"

  use CONSTANTS

  IMPLICIT NONE

  REAL*8 :: ES_ICE
  REAL*8 :: T

  ES_ICE = ES0*EXP( (LS0/R_V)*(1.D0/T0 - 1.D0/T) )

END FUNCTION ES_ICE

FUNCTION S_LIQ(T,P,Q)

  !SUPERSATURATION OVER LIQUID

  IMPLICIT NONE

  REAL*8 :: S_LIQ
  REAL*8 :: T,P,Q
  REAL*8 :: E_VAP,ES_LIQ

  S_LIQ = E_VAP(Q,P)/ES_LIQ(T) - 1.D0

END FUNCTION S_LIQ

FUNCTION RV_S(T,P)

  !Vapor mixing ratio at saturation

  use CONSTANTS

  IMPLICIT NONE

  REAL*8 :: RV_S
  REAL*8 :: T,P
  REAL*8 :: ES_LIQ
  REAL*8 :: EPS

  EPS = R_D/R_V

  RV_S = EPS*ES_LIQ(T)/( P - ES_LIQ(T) )

END FUNCTION RV_S

SUBROUTINE GET_ENVIRONMENT_PROFILES_DYCOMS_II

    !Follows Stevens et al. (2005) with an inversion at z_inv

    use ENVIRONMENT
    use GRID
    use CONSTANTS
    use FUNCTIONS
    use ADVECTION, only: eps
    IMPLICIT NONE

    !REAL*8  :: R_MOIST
    REAL*8  :: EE
    REAL*8  :: ES_LIQ,S
    REAL*8  :: EPSILON
    real*8, allocatable :: z_table(:),theta_table(:),rv_table(:)
    real*8, allocatable :: Z(:)
    INTEGER*8 :: K
    
    !Mixed-phase profiles
    allocate(z_table(6), theta_table(6), rv_table(6))
    z_table     = (/0.D0  , 500.D0 , 900.D0, 1100.D0, 1150.D0, 1800.D0 /)
    theta_table = (/264.D0, 265.5D0, 266.D0, 266.5D0, 271.D0 , 276.D0  /)
    rv_table    = (/1.8D0 , 1.45D0 , 1.3D0 , 1.2D0  , 1.65D0 , 1.D0    /)/1000
   !rv_table    = (/1.8D0 , 1.45D0 , 1.3D0 , 1.2D0  , 1.15D0 , 0.7D0   /)


    !Warm Profile
    !z_table     = (/0.D0, 840.D0, 841.D0, 850.D0, 900.D0, 1250.D0 /)
    !theta_table = (/289.7D0, 298.7D0, 302.4D0, 308.2D0, 311.85D0, 370.D0/)
    !rv_table    = (/17.0D0, 16.3D0, 10.7D0, 4.2D0, 3.0D0, 0.1D0/)/1000.D0

    !Paper do Abade
    ! z_table = linspace(0.D0,LZ,1251)
    ! allocate( theta_table(size(z_table)) , rv_table(size(z_table)) )
    ! do k = 1,size(z_table)
    !     if (z_table(k) < z_inv) then
    !         theta_table(k) = 289.D0
    !         rv_table(k)    = 8.9D0/1000
    !     else
    !         theta_table(k) = 289.D0 + (z_table(k) - z_inv)**(1.D0/3)
    !         rv_table(k)    = 1.5D0/1000
    !     end if
    ! end do

    Z = linspace(DZ/2.0,(NZ-0.5D0)*DZ,int(NZ,4))

    ALLOCATE ( TH_E(NZ), RV_E(NZ) )
    TH_E = interpol(z_table,theta_table,Z)
    RV_E = interpol(z_table,rv_table,Z)



    ! OLD PROFILES--------------------------------
    !TH_E = 289.D0
    !RV_E = 9.D-3  !7.5D-3

    !--------------------------------------------
    !Inversion layer
    !--------------------------------------------

    !DO K = 1,NZ
        ! Z = DZ*K - 0.5D0*DZ
        ! IF (Z.GT.Z_INV) THEN
        !    TH_E(K) = 289.D0 + (Z-Z_INV)**(1.D0/3.D0)
        !    RV_E(K) = 1.5D-3
        ! END IF
    !END DO

    ! OLD PROFILES-------------------------------

    CALL GET_ENV_PRESSURE_TEMP

    !Base state density profile

    ALLOCATE ( RHO(NZ) )

    DO K = 1,NZ
        !Total density (dry air + water vapor)
        !RHO(K) = PR_E(K) / ( R_MOIST(RV_E(k))*T_E(K) )

        !Partial density of dry air
        EE = PR_E(K)*( RV_E(K)/( RV_E(K) + R_D/R_V ) ) !Partial vapor pressure
        RHO(K) = (PR_E(K)-EE)/(R_D*T_E(K))
    END DO

    ! * * *

    OPEN (UNIT = 99, FILE = "environ_profile.out")
    DO K = 1,NZ
        EE = PR_E(K)*( RV_E(K)/( RV_E(K) + R_D/R_V ) ) !Partial vapor pressure
        S  = EE/ES_LIQ(T_E(K)) - 1.D0 !Supersaturation

        EPSILON = 1.D-3
        IF (Z(K).GT.Z_INV) THEN
            EPSILON = 0.D0
        END IF

        WRITE(99,*) Z(K), TH_E(K), RV_E(K)*1.D3, PR_E(K), T_E(K), RHO(K), S, EPSILON, 0.D0

    END DO
    CLOSE (99, STATUS = "KEEP")

    ! Building the eddy diffusivity profile
    allocate(EDDY_DIFF(NZ))
    TT        = DZ 
    EDDY_MAX  = (1.D2**(4.D0/3))*eps**(1.D0/3) ! Melhorar: Importar DELTA e EPS do gerador de velocidades.
    EDDY_DIFF = 0.5D0*EDDY_MAX*ERFC((Z - Z_INV)/(TT*sqrt(2.D0)))
  ! * * *

END SUBROUTINE GET_ENVIRONMENT_PROFILES_DYCOMS_II

SUBROUTINE GET_ENV_PRESSURE_TEMP

   use ENVIRONMENT
   use GRID
   use CONSTANTS

   implicit none
   INTEGER*8 :: K

   REAL*8 :: KAPPA
   REAL*8 :: THV,THV_P,THV_M

   KAPPA = R_D/CP_D

   ALLOCATE ( PR_E(NZ) )
   ALLOCATE (  T_E(NZ) )
   ALLOCATE (  Z_E(NZ) )

   DO K = 1,NZ
      Z_E(K) = (K - 0.5D0)*DZ
   END DO

   PR_E(1) = 1.015D5

   DO K = 2,NZ

      THV_M = TH_E(K-1)*( 1 + 0.608D0*RV_E(K-1) )
      THV_P = TH_E(K)  *( 1 + 0.608D0*RV_E(K)   )
      THV = 0.5D0*( THV_M + THV_P )

      PR_E(K) = ( PR_E(K-1)**KAPPA - ( (G*P00**KAPPA)/(CP_D*THV) )*(Z_E(K)-Z_E(K-1)) )**(1.D0/KAPPA)
   END DO

  ! * * *

   DO K = 1,NZ
      T_E(K) = TH_E(K)*(PR_E(K)/P00)**KAPPA
   END DO

END SUBROUTINE GET_ENV_PRESSURE_TEMP