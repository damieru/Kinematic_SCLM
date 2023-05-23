module HDF5_functions


   implicit none

contains

   SUBROUTINE HDF5_CREATE_DATA_SET(dlabel, dsize, dtype_in)

      use HDF5_VARIABLES

      character(len=*)   , intent(in) :: dlabel
      integer            , intent(in) :: dsize(:)
      integer*8, optional, intent(in) :: dtype_in
      integer   :: drank, err
      integer*8 :: dtype
      
      drank = size(dsize)

      if (present(dtype_in)) then
         dtype = dtype_in
      else
         dtype = H5T_NATIVE_DOUBLE
      end if

      call h5screate_simple_f(drank, int(dsize,8), dspace_id, err)      ! Create dataspace
      call h5dcreate_f(file_id, dlabel, dtype, dspace_id, dset_id, err) ! Create dataset
      call h5dclose_f(dset_id,err)                                      ! Close dataset
      call h5sclose_f(dspace_id,err)                                    ! Close dataspace
      
   END SUBROUTINE HDF5_CREATE_DATA_SET

   SUBROUTINE HDF5_WRITE_DOUBLE(data, dset_name, time_offset)

      use HDF5_VARIABLES, only: memspace, error, file_id, dset_id, dspace_id
      use HDF5

      real*8           , intent(in) :: data(..)
      character(len=*) , intent(in) :: dset_name
      integer, optional, intent(in) :: time_offset

      integer                       :: dspace_rank
      integer*8       , allocatable :: dcount(:), offset(:)
      integer(HSIZE_T), allocatable :: stride(:), block_size(:)
      

      if (present(time_offset)) then
         dspace_rank = rank(data) + 1
         dcount  = [shape(data)  , 1]
         offset = [shape(data)*0, time_offset]
         stride = dcount*0 + 1
         block_size = stride;
      else 
         dspace_rank = rank(data)
         dcount = shape(data)
         offset = shape(data)*0
         stride = dcount*0 + 1
         block_size = stride
      end if

      !Create memory data space
      call h5screate_simple_f(dspace_rank, dcount, memspace, error)
      !Open THETA dataset
      call h5dopen_f(file_id, dset_name, dset_id, error)
      call h5dget_space_f(dset_id, dspace_id, error)
      !Select subset
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, dcount, error, stride, block_size)
      !Write subset to dataset (can't do without the "select rank" statement because of the shitty h5dwrite interface)
      select rank(data)
         rank(0)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dcount, error, memspace, dspace_id)
         rank(1)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dcount, error, memspace, dspace_id)
         rank(2)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dcount, error, memspace, dspace_id)
         rank default
            stop 'Data rank must be either 0, 1 or 2'
      end select
      !Close dataset
      call h5dclose_f(dset_id, error)
      !Close memory dataspace
      call h5sclose_f(memspace , error)

   END SUBROUTINE HDF5_WRITE_DOUBLE

   SUBROUTINE HDF5_WRITE_INTEGER(data, dset_name, time_offset)

      use HDF5_VARIABLES, only: memspace, error, file_id, dset_id, dspace_id
      use HDF5

      integer          , intent(in) :: data(..)
      character(len=*) , intent(in) :: dset_name
      integer, optional, intent(in) :: time_offset

      integer                       :: dspace_rank
      integer*8       , allocatable :: dcount(:), offset(:)
      integer(HSIZE_T), allocatable :: stride(:), block_size(:)
      

      if (present(time_offset)) then
         dspace_rank = rank(data) + 1
         dcount  = [shape(data)  , 1]
         offset = [shape(data)*0, time_offset]
         stride = dcount*0 + 1
         block_size = stride;
      else 
         dspace_rank = rank(data)
         dcount = shape(data)
         offset = shape(data)*0
         stride = dcount*0 + 1
         block_size = stride
      end if

      !Create memory data space
      call h5screate_simple_f(dspace_rank, dcount, memspace, error)
      !Open THETA dataset
      call h5dopen_f(file_id, dset_name, dset_id, error)
      call h5dget_space_f(dset_id, dspace_id, error)
      !Select subset
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, dcount, error, stride, block_size)
      !Write subset to dataset (can't do without the "select rank" statement because of the shitty h5dwrite interface)
      select rank(data)
         rank(0)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dcount, error, memspace, dspace_id)
         rank(1)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dcount, error, memspace, dspace_id)
         rank(2)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dcount, error, memspace, dspace_id)
         rank default
            stop 'Data rank must be either 0, 1 or 2'
      end select
      !Close dataset
      call h5dclose_f(dset_id, error)
      !Close memory dataspace
      call h5sclose_f(memspace , error)

   END SUBROUTINE HDF5_WRITE_INTEGER

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
      call HDF5_CREATE_DATA_SET('DPB'  ,[NX,NZ,NT], H5T_NATIVE_INTEGER)

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
      use STEPPER          , only: TIME, DT
      use ADVECTION        , only: UXN,UZN
      use THERMODYNAMIC_VAR, only: THETA, RV, RL, RI, DPB, SAT
      use ENVIRONMENT      , only: RHO
      use SD_VARIABLES     , only: X, R
      use IO_PARAMETERS
      IMPLICIT NONE

      REAL*8, intent(in) :: NSEC
      INTEGER            :: iteration, period, time_offset
      CHARACTER(100)     :: file_path
      logical, save      :: first_call = .true.

      iteration = NINT( TIME / DT ) 
      period    = NINT( NSEC / DT )

      if (MOD(iteration, period) /= 0) return

      time_offset = iteration / period

      file_path = TRIM(output_folder)//'/'//TRIM(output_file)

      call h5open_f(error) !Open fortran HDF5 interface
      call h5fopen_f(file_path,H5F_ACC_RDWR_F,file_id,error) !Open file

      ! Save RHO profile only once
      if (first_call) then
            call HDF5_WRITE_DOUBLE(RHO ,'RHO')
            first_call = .false.
      end if

      call HDF5_WRITE_DOUBLE(TIME  ,'TIME' ,time_offset)

      call HDF5_WRITE_DOUBLE (THETA ,'THETA',time_offset)
      call HDF5_WRITE_DOUBLE (RV    ,'RV'   ,time_offset)
      call HDF5_WRITE_DOUBLE (RL    ,'RL'   ,time_offset)
      call HDF5_WRITE_DOUBLE (RI    ,'RI'   ,time_offset)
      call HDF5_WRITE_DOUBLE (SAT   ,'SAT'  ,time_offset)
      call HDF5_WRITE_INTEGER(DPB   ,'DPB'  ,time_offset)

      call HDF5_WRITE_DOUBLE (UXN   ,'UXN'  ,time_offset)
      call HDF5_WRITE_DOUBLE (UZN   ,'UZN'  ,time_offset)

      call HDF5_WRITE_DOUBLE (  X   ,'X_SD' ,time_offset)
      call HDF5_WRITE_DOUBLE (R(2,:),'R_SD' ,time_offset)

      call h5fclose_f(file_id  , error)
      call h5close_f(error)
   
   END SUBROUTINE HDF5_SAVE_RESULTS

   SUBROUTINE HDF5_READ_DOUBLE(data, dset_name, time_offset)

      use HDF5_VARIABLES, only: file_id, dset_id, dspace_id, memspace, error
      use HDF5

      real*8           , intent(out) :: data(..)
      character(len=*) , intent(in ) :: dset_name 
      integer, optional, intent(in ) :: time_offset

      integer                        :: dspace_rank
      integer(HSIZE_T) , allocatable :: dcount(:), offset(:)
      integer(HSIZE_T) , allocatable :: stride(:), block_size(:)

      !dcount = int([NXP, NZ, 1],8)
      ! rnk = 3

      if (present(time_offset)) then
         dspace_rank = rank(data) + 1
         dcount  = [shape(data)  , 1]
         offset = [shape(data)*0, time_offset]
         stride = dcount*0 + 1
         block_size = stride;
      else 
         dspace_rank = rank(data)
         dcount = shape(data)
         offset = shape(data)*0
         stride = dcount*0 + 1
         block_size = stride
      end if
      
      call h5dopen_f(file_id, dset_name, dset_id, error)             !Open dataset
      call h5dget_space_f(dset_id, dspace_id, error)                 !Get dataspace id
      call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, dcount, error, stride, block_size) !Select subset      
      call h5screate_simple_f(dspace_rank, dcount, memspace, error)  !Create memory data space
      select rank(data)                                              !Read subset from dataset
         rank(0)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dcount, error, memspace, dspace_id)
         rank(1)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dcount, error, memspace, dspace_id)
         rank(2)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, dcount, error, memspace, dspace_id)
         rank default
            stop 'Data rank must be either 0, 1 or 2'
      end select
      call h5dclose_f(dset_id, error)     !Close dataset
      call h5sclose_f(dspace_id, error)   !Close dataspace
      call h5sclose_f(memspace , error)   !Close memory dataspace

   END SUBROUTINE HDF5_READ_DOUBLE

   SUBROUTINE HDF5_READ_VELOCITIES(filename)

      use HDF5_VARIABLES
      use STEPPER
      use GRID
      use ADVECTION, ONLY: UXT, UZT

      implicit NONE
      

      character(100), intent(in) :: filename
      integer, save :: call_count = 0, offset_max
      integer :: time_offset
      
      if (call_count==0) offset_max = nint(T_MAX/DT_ADV) + 1
      time_offset = min(call_count, offset_max)
      call_count  = call_count + 1

      call h5open_f(error) !Open fortran HDF5 interface
      call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error) !Open file
      if (error /= 0) then
         call system('clear')
         stop 'HDF5_READ_VELOCITIES: Velocities file not found.'
      end if

      call HDF5_READ_DOUBLE(UXT,'UXT',time_offset)
      call HDF5_READ_DOUBLE(UZT,'UZT',time_offset)

      !Close file and interface
      call h5fclose_f(file_id, error)
      call h5close_f(error)
   END SUBROUTINE HDF5_READ_VELOCITIES

   SUBROUTINE HDF5_READ_PSI(filename)
      use HDF5_VARIABLES 
      use STEPPER
      use GRID
      use ADVECTION, only: PSI

      implicit NONE
     
      character(100), intent(in) :: filename
      integer, save :: call_count = 0, offset_max
      integer       :: time_offset 
      

      if (call_count==0) offset_max = nint(T_MAX/DT_ADV) + 1
      time_offset = min(call_count, offset_max)
      call_count  = call_count + 1

      call h5open_f(error) !Open fortran HDF5 interface
      call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,error) !Open file
      if (error /= 0) then
         call system('clear')
         stop 'HDF5_READ_PSI: file not found.'
      end if
      
      call HDF5_READ_DOUBLE(PSI,'PSI',time_offset)

      !Close file and interface
      call h5fclose_f(file_id, error)
      call h5close_f(error)
      
   END SUBROUTINE HDF5_READ_PSI

end module HDF5_functions