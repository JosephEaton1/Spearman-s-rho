program permanent_polynomial_mod
	!This code computes frequency table of S stat modulo some m (you can change m) for
	!given range of n
    implicit none
    integer, parameter :: i64 = selected_int_kind(18)

    integer :: n, n_min, n_max
    integer :: i, j, l, deg, max_prev_deg, max_new_deg, max_deg, s
    integer(i64), allocatable :: row_prod(:), temp_poly(:), result_poly(:), delta(:)
    integer(i64) :: sign_val, m, k, total, fact
    integer :: flip_row, csv_unit
    character(len=256) :: filename, arg
    character(len=4096) :: line

    n_min = 18
    n_max = 26
    m = 2_i64**57 - 1_i64       !must be odd to allow division later
				!and m*n< n^63 to stay in 64bit

    do n = n_min, n_max

        write(*, '(A, I3, A)') "  Starting n =", n, " ..."
        max_deg = n * n * (n + 1) / 2

        allocate(row_prod(0:max_deg), temp_poly(0:max_deg), &
                 result_poly(0:max_deg), delta(n))

        result_poly = 0_i64
        delta(1)    = 1_i64  !d_1=1, then use trailz function to generate
        do i = 2, n          !gray code (1 flip per turn)
            delta(i) = -1_i64
        end do
        sign_val = 1_i64
        if (mod(n - 1, 2) == 1) sign_val = -1_i64 !signval iterates each flip

        do k = 0_i64, 2_i64**(n - 1) - 1_i64

            if (k > 0_i64) then
                flip_row        = trailz(k) + 2
                delta(flip_row) = -delta(flip_row)
                sign_val        = -sign_val
            end if

            row_prod    = 0_i64
            row_prod(0) = 1_i64  !otherwise will collapse to 0

            do j = 1, n
                max_prev_deg = n * (j - 1) * j / 2
                max_new_deg  = n * j * (j + 1) / 2
                temp_poly    = 0_i64
                do i = 0, max_prev_deg 
                    if (row_prod(i) == 0_i64) cycle
                    do l = 1, n
                        deg            = i + l * j
                        temp_poly(deg) = temp_poly(deg) + row_prod(i) * delta(l)
                    end do
                end do
				!each degree recieved at most n of size m so < 2^63
                do deg = 0, max_new_deg
                    temp_poly(deg) = modulo(temp_poly(deg), m)
                end do
                row_prod = temp_poly
            end do

            do deg = 0, max_deg  !here can do addition inside modulo as both <m so sum < 2^63
                result_poly(deg) = modulo(result_poly(deg) + sign_val * row_prod(deg), m)
            end do

        end do

        !to divide we half n-1 times, to half find h st 2h=rmod m
	!so for even r, h = r/2, for odd r, h = r+m/2 
        do s = 1, n - 1
            do deg = 0, max_deg
                if (mod(result_poly(deg), 2_i64) == 0_i64) then
                    result_poly(deg) = result_poly(deg) / 2_i64
                else
                    result_poly(deg) = (result_poly(deg) + m) / 2_i64
                end if
            end do
        end do
       
	write(*, '(A, I3, A)') "  Completed n =", n, ". Writing CSV..."

	!Following code is only relevant to where to store data, edit ** code accordingly

        write(filename, '(A, I0, A, I0, A)') &
            '**The file branch you want it stored in**', n, '_m', m, '.csv'

        csv_unit = 20
        open(unit=csv_unit, file=trim(filename), status='replace', action='write')
        write(csv_unit, '(A)') 'degree,residue'
        do deg = 0, max_deg
            if (result_poly(deg) /= 0_i64) then
                write(line, '(I0,A,I0)') deg, ',', result_poly(deg)
                write(csv_unit, '(A)') trim(line)
            end if
        end do
        close(csv_unit)

        write(*, '(A, A)') "  Saved: ", trim(filename)
        deallocate(row_prod, temp_poly, result_poly, delta)

    end do

    print *, "All done."

end program permanent_polynomial_mod

