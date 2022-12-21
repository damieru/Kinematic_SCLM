

program write_input
    use String_Functions
    
    implicit none
    real*8 :: Delta(2), I_Frac(4)
    character(len= 50)  :: name, out_name
    character(len=100)  :: lish
    character(len= 50)  :: VAE(3)
    integer :: i, j, n
    integer :: lines(3), linha_atual, IO_STATUS = 0
    

    !Used filter size
    Delta = (/0.D0, 50.D0/)
    !Used Ice Fractions
    I_Frac = (/0.001D0, 0.01D0, 0.05D0, 0.1D0/)
    ! Lines to be editted
    lines = (/23, 32, 41/)
    
    do i = 1,size(Delta)
        do j = 1,size(I_Frac)
            name     = GET_FILE_NAME(Delta(i),I_Frac(j))
            out_name = Replace_Text(name,'.in','.h5')
            write(VAE(1),   '(A)'   ) out_name
            if (Delta(i) < 10) then
                write(VAE(2),   '(F3.1)') Delta(i)
            else
                write(VAE(2),   '(F4.1)') Delta(i)
            end if
            write(VAE(3),'(ES6.1E1)') I_Frac(j)            
            open(unit=2,file='warm_no_var_50_5050.in')
            open(unit=1,file=name)
            linha_atual = 1
            n = 1
            IO_STATUS = 0
            do while (IO_STATUS == 0)
                read(2,'(A)',IOSTAT=IO_STATUS) lish
                if (IO_STATUS==0) then
                    if (linha_atual == lines(n)) then
                        if (n == 1) then
                            write(1,'(A)') '/media/daniel/Bkp Mint/Results' 
                            write(1,'(A)') VAE(n)
                            linha_atual = linha_atual + 1
                            read(2,'(A)',IOSTAT=IO_STATUS) lish
                        else
                            write(1,'(A)') VAE(n)(1:12)//lish(13:len(lish))
                        end if
                        n = n + 1
                    else
                        write(1,'(A)') lish
                    end if
                end if
                linha_atual = linha_atual + 1
            end do
            close(1)
            close(2)
        end do
    end do

    print*,'Files created successfully.'
    call system('code ice*in')


    contains
        function GET_FILE_NAME(Delta,I) result(file_name)
            real*8, intent(in) :: Delta, I
            character(len=50) :: file_name
            if (Delta < 10) then
                write(file_name,'(A,F3.1,A,F5.3,A)') 'ice_D',Delta,'_i',I,'.in'
            else
                write(file_name,'(A,F4.1,A,F5.3,A)') 'ice_D',Delta,'_i',I,'.in'
            end if
        end function  
end program 