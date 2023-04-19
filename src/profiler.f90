module profiler
   
   implicit none

contains

   subroutine tictoc(section_name, finish,time_unit)

      logical         , optional, intent(in) :: finish
      character(len=*), optional, intent(in) :: section_name
      character(len=*), optional, intent(in) :: time_unit

      character(len=100), save, allocatable :: names(:)
      character(len= 7 ), save              :: valid_units(3)
      real*8            , save              :: conv_factors(3)
      real*8            , save, allocatable :: elapsed(:)
      integer           , save, allocatable :: counter(:)
      logical           , save :: first_call = .true.
      integer           , save :: rate
      integer           , save :: time_now, time_old = 0
      logical                  :: match
      real*8 , allocatable :: avg_times(:)
      integer:: i
      integer, save :: j

      ! Some checks
      if (.not.present(section_name).and.present(finish)) then
         stop "Cant' finish without a section name."
      elseif ( (.not.present(section_name)) .and. (.not.first_call) ) then
         stop "This is not the first call. You have to provide a section name."
      elseif (.not.first_call .and. present(time_unit)) then
         stop "Time unit should be set in the first call only."
      end if
      
      call system_clock(time_now)
      
      if (first_call) then
         call system_clock(count_rate=rate)
         allocate(names(0),elapsed(0),counter(0))
         valid_units(1)  = 's'
         valid_units(2)  = 'ms'
         valid_units(3)  = 'micro_s'
         conv_factors = [1, 1000, 1000000]
         !Determine conversion factor based on chosen time unit
         if (present(time_unit)) then
            match = .false.
            do j = 1,size(valid_units)
               if (time_unit == valid_units(j)) then
                  match = .true.
                  exit
               end if
            end do
            if (.not.match) stop "Invalid time unit. Please choose: s, ms or micro_s"
         else
            j = 1;
         end if
         first_call = .false.
      end if
      
      if (present(section_name)) then
         !Check if section already exists
         match = .false.
         if (size(names)/=0) then
            do i = 1,size(names)
               if (section_name == names(i)) then
                  match = .true.
                  exit 
               end if
            end do
         end if
         
         if (.not.match) then
            names   = [names, section_name]
            elapsed = [elapsed, real(time_now-time_old,8)/rate]  
            counter = [counter, 1]  
         else
            elapsed(i) = elapsed(i) + real(time_now-time_old,8)/rate
            counter(i) = counter(i) + 1
         end if
      end if
      
      if (present(finish)) then
         if (finish) then
            allocate( avg_times(size(names)) )
            
            !Compute average times
            do i = 1,size(names)
               avg_times(i) = elapsed(i)/counter(i)
            end do

            !Save average times into a file
            open(unit=42, file='tictoc.txt')
            write(42,*) 'Average time per section in '//trim(valid_units(j))//':'
            do i = 1,size(names)
               write(42,'(A,F0.2)') trim(names(i))//char(9), avg_times(i)*conv_factors(j)
            end do
            write(42,'(A,F0.2,A)') 'Average time per cycle: ', sum(avg_times)*conv_factors(j), trim(valid_units(j))//' on average.'
            write(42,'(A,F0.2,A)') 'Total elapsed time: ', sum(elapsed)*conv_factors(j), trim(valid_units(j))//'.'
            close(42)
         end if
      end if

      time_old = time_now
   end subroutine tictoc

end module profiler