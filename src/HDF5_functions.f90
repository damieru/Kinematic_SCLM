module HDF5_functions

   use HDF5_VARIABLES

   implicit none

contains

   SUBROUTINE HDF5_CREATE_DATA_SET(dlabel, dsize, int_data)

      character(len=*), intent(in)  :: dlabel
      integer, intent(in)           :: dsize(:)
      logical, optional, intent(in) :: int_data
      integer   :: drank, err
      integer*8 :: dtype
      
      drank = size(dsize)

      if (present(int_data).and.int_data) then
         dtype = H5T_NATIVE_INTEGER
      else
         dtype = H5T_NATIVE_DOUBLE
      end if

      call h5screate_simple_f(drank, int(dsize,8), dspace_id, err)      ! Create dataspace
      call h5dcreate_f(file_id, dlabel, dtype, dspace_id, dset_id, err) ! Create dataset
      call h5dclose_f(dset_id,err)                                      ! Close dataset
      call h5sclose_f(dspace_id,err)                                    ! Close dataspace
      
   END SUBROUTINE HDF5_CREATE_DATA_SET

   SUBROUTINE HDF5_CREATE_FILE

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

      NT = NINT(T_MAX/OUT_PER) + 1 ! Number of saved instants
      call h5open_f(error) !Opens HDF5 Fortran interface
      call h5fcreate_f(file_path, H5F_ACC_TRUNC_F, file_id, error) !Creates file

      call HDF5_CREATE_DATA_SET('TIME' ,[NT]      )
      call HDF5_CREATE_DATA_SET('RHO'  ,[NZ]      )
      call HDF5_CREATE_DATA_SET('THETA',[NX,NZ,NT])
      call HDF5_CREATE_DATA_SET('RV'   ,[NX,NZ,NT])
      call HDF5_CREATE_DATA_SET('RL'   ,[NX,NZ,NT])
      call HDF5_CREATE_DATA_SET('RI'   ,[NX,NZ,NT])
      call HDF5_CREATE_DATA_SET('SAT'  ,[NX,NZ,NT])
      call HDF5_CREATE_DATA_SET('DPB'  ,[NX,NZ,NT],int_data=.true.)

      call HDF5_CREATE_DATA_SET('UXN'  ,[NXP,NZP,NT])
      call HDF5_CREATE_DATA_SET('UZN'  ,[NXP,NZP,NT])

      call HDF5_CREATE_DATA_SET('X_SD' ,[N_sd,2 ,NT])
      call HDF5_CREATE_DATA_SET('R_SD' ,[N_sd,NT])
     

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
               !Close memspace
               call h5sclose_f(memspace, error)
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

end module HDF5_functions