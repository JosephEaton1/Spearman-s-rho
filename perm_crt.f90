program perm_crt
	!This code converts frequency tables modulo m1 and m2 into full tables
    implicit none
    integer, parameter :: i64 = selected_int_kind(18)
    integer, parameter :: qp  = selected_real_kind(33)

    integer(i64), parameter :: M1 = 144115188075855871_i64   ! 2^57 - 1
    integer(i64), parameter :: M2 = 144115188075855869_i64   ! 2^57 - 3

    integer :: n, deg, max_deg, ios
    integer(i64) :: val, t
    integer(i64), allocatable :: r1(:), r2(:)
    real(qp) :: x
    character(len=256) :: path, fname
    character(len=64)  :: numstr

    path = '/exports/eddie/scratch/s2642520/perm_poly/' !you will have to change path to where
							!you stored data

    do n = 18, 26

        max_deg = n * n * (n + 1) / 2
        allocate(r1(0:max_deg), r2(0:max_deg))
        r1 = 0_i64
        r2 = 0_i64

        !this reads residue modulo m1 treating missing frequencies as 0, again here you will
	!have to change path accordingly
        write(fname, '(A,A,I0,A,I0,A)') trim(path), 'mod/permanent_poly_n', n, '_m', M1, '.csv'
        open(10, file=trim(fname), status='old', action='read')
        read(10, *)                                  ! skip header
        do
            read(10, *, iostat=ios) deg, val
            if (ios /= 0) exit
            r1(deg) = val
        end do
        close(10)

        !same for m2
        write(fname, '(A,A,I0,A,I0,A)') trim(path), 'mod/permanent_poly_n', n, '_m', M2, '.csv'
        open(11, file=trim(fname), status='old', action='read')
        read(11, *)
        do
            read(11, *, iostat=ios) deg, val
            if (ios /= 0) exit
            r2(deg) = val
        end do
        close(11)
	!again here change to where you want to store on your computer

        write(fname, '(A,A,I0,A)') trim(path), 'final/permanent_poly_n', n, '.csv'
        open(20, file=trim(fname), status='replace', action='write')
        write(20, '(A)') 'degree,coefficient'

        do deg = 0, max_deg
            if (r1(deg) == 0_i64 .and. r2(deg) == 0_i64) cycle

            
	    ! t = (r2 - r1) * M1^(-1) mod M2  =  (r2 - r1) / 2 mod M2 (simplified garner)
            
	    t = modulo(r2(deg) - r1(deg), M2)
            if (mod(t, 2_i64) /= 0_i64) t = t + M2
            t = t / 2_i64

            ! x = r1 + M1 * t 
            x = real(r1(deg), qp) + real(M1, qp) * real(t, qp)

            ! this prints as an integer instead of floating point in qp
            write(numstr, '(F45.0)') x
            numstr = adjustl(numstr)
            numstr = numstr(1:len_trim(numstr) - 1)
            write(20, '(I0,A,A)') deg, ',', trim(numstr)
        end do

        close(20)
        deallocate(r1, r2)
        write(*, '(A,I0)') "Done n = ", n

    end do

end program perm_crt
