Module utils
    use, intrinsic :: ieee_arithmetic
    implicit none

Contains

    subroutine digital_filter(sequence, b_yin, b_zin, filtered, s_shape)
      implicit none
      integer, parameter :: dp = kind(1.0d0)
      real(dp), dimension(:,:,:), intent(in) :: sequence
      real(dp), dimension(:), intent(in) :: b_yin, b_zin
      real(dp), dimension(:,:,:), intent(out) :: filtered
      integer, dimension(3), intent(in) :: s_shape

      integer :: N_y, N_z
      integer :: i, j, k, ii, jj
      integer :: seq_y, seq_z, seq_k
      real(dp) :: sum

      ! Get filter lengths
      N_y = size(b_yin) / 2
      N_z = size(b_zin) / 2

      ! Get sequence dimensions
      seq_y = s_shape(1)
      seq_z = s_shape(2)
      seq_k = s_shape(3)

      ! Ensure the filter coefficients are valid
      if (size(b_yin) <= 0 .or. size(b_zin) <= 0) then
          print *, "Error: Filter coefficients size must be positive."
          error stop
      end if

      ! Initialize the filtered array to zero
      filtered = 0.0_dp

      do k = 1, seq_k
          do j = 1, seq_z
              do i = 1, seq_y
                  sum = 0.0_dp

                  ! Apply the filter, accounting for edge conditions
                  do jj = -N_z, N_z
                      do ii = -N_y, N_y
                          ! Check bounds for sequence array
                          if ((i + ii >= 1 .and. i + ii <= seq_y) .and. &
                              (j + jj >= 1 .and. j + jj <= seq_z)) then                      
                                sum = sum + b_yin(N_y + 1 + ii) * b_zin(N_z + 1 + jj) * sequence(i + ii, j + jj, k)
                          end if
                      end do
                  end do

                  ! Assign the result to the filtered array
                  filtered(i, j, k) = sum

                  
              end do
          end do
      end do

    end subroutine digital_filter 


    subroutine generate_normal_random_numbers(mean, stddev, n, random_array)
      implicit none
      integer, parameter :: dp = kind(1.0d0)
      real(dp), intent(in) :: mean
      real(dp), intent(in) :: stddev
      integer, intent(in) :: n
      real(dp), dimension(:), intent(out) :: random_array

      integer :: i
      real(dp) :: u1, u2
      real(dp), dimension(:), allocatable :: uniform_random_numbers

      ! Define pi
      real(dp) :: pi

      allocate(uniform_random_numbers(2 * n))
      
      ! Initialize the random seed
      call random_seed()

      ! Generate uniform random numbers
      call random_number(uniform_random_numbers)
      
      ! Define pi with the correct precision
      pi = acos(-1.0_dp)

      do i = 1, n, 2
          u1 = uniform_random_numbers(i)
          u2 = uniform_random_numbers(i + 1)
          
          ! Avoid log(0) and invalid operations
          if (u1 <= 0.0_dp) then
              print *, "Warning: u1 is non-positive. Adjusting value to avoid log(0)."
              u1 = 1.0_dp / 1.0d-30
          end if
          
          ! Apply the Box-Muller transform
          random_array(i) = mean + stddev * sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
          if (i + 1 <= n) then
              random_array(i + 1) = mean + stddev * sqrt(-2.0_dp * log(u1)) * sin(2.0_dp * pi * u2)
          endif
      end do

      deallocate(uniform_random_numbers)
    end subroutine generate_normal_random_numbers


      subroutine cholesky(a,r,n) 
        implicit none
        integer, parameter :: dp = kind(1.0d0)
        integer(dp) :: n
        real(dp),intent(out) :: a(n,n)
        real(dp),intent(in) :: r(n,n)

        integer :: i,j,k
        a(1:n,1:n) = 0.0

        do i = 1,n
          do j = 1,i
            a(i,j) = r(i,j)

            if (i==j) then
              do k = 1,j-1
                a(i,j) = a(i,j) - a(j,k)**2
              end do
              a(i,j) = sqrt(abs(a(i,j)))

            else
              do k = 1,j-1
                a(i,j) = a(i,j) - a(i,k)*a(j,k)
              end do
              a(i,j) = a(i,j)/(a(j,j) + 1.0e-16)

            end if

          end do
        end do

      end subroutine cholesky

      subroutine matrix_vector_dot_product(matrix, vector, result)
        implicit none
        integer, parameter :: dp = kind(1.0d0)
        real(dp), dimension(3,3), intent(in) :: matrix
        real(dp), dimension(3), intent(in) :: vector
        real(dp), dimension(3), intent(out) :: result
        integer(dp) :: i, j
    
        ! Initialize the result vector to zero
        result = 0.0d0
    
        ! Compute the dot product
        do i = 1, 3
          do j = 1, 3
            result(i) = result(i) + matrix(i,j) * vector(j)
          end do
        end do

      end subroutine matrix_vector_dot_product

      subroutine read_file_skip_first(filename, data, num_columns)
        implicit none
        integer, parameter :: dp = kind(1.0d0)
        character(len=*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: data(:,:)
        integer(dp), intent(out) :: num_columns
    
        integer(dp) :: i, io_status, num_rows, unit
        character(len=256) :: line
        integer(dp) :: num_items
        real(dp), dimension(:), allocatable :: temp_row
    
        ! Open the file
        unit = 10
        open(unit=unit, file=filename, status='old', action='read')
    
        ! Skip the first row
        read(unit, '(A)', iostat=io_status) line
        if (io_status /= 0) then
            print *, 'Error reading file: ', io_status
            stop
        endif
    
        ! Read the second row to determine the number of columns
        read(unit, '(A)', iostat=io_status) line
        if (io_status /= 0) then
            print *, 'Error reading file: ', io_status
            stop
        endif
    
        ! Count number of columns (whitespace-separated values)
        num_columns = 0
        num_items = 0
        do i = 1, len_trim(line)
          if (line(i:i) /= ' ' .and. (i == 1 .or. line(i-1:i-1) == ' ')) then
            num_columns = num_columns + 1
          endif
        enddo
    
        ! Allocate temporary row storage
        allocate(temp_row(num_columns))
    
        ! Read the entire file to determine the number of rows
        num_rows = 0
        do
          read(unit, *, iostat=io_status) temp_row
          if (io_status /= 0) exit
          num_rows = num_rows + 1
        enddo
        ! Supplement the total number of lines with an additional line
        num_rows = num_rows + 1

        ! Allocate the data array
        allocate(data(num_rows, num_columns))
    
        ! Rewind the file and skip the first row again
        rewind(unit)
        read(unit, '(A)') line
    
        ! Read the actual data
        i = 0
        do
          read(unit, *, iostat=io_status) temp_row
          if (io_status /= 0) exit
          i = i + 1
          data(i, :) = temp_row
        enddo
    
        ! Close the file
        close(unit)
    
        ! Deallocate temporary storage
        deallocate(temp_row)
      end subroutine read_file_skip_first

      subroutine write_2d_array_to_file(filename, array)
        implicit none
        integer, parameter :: dp = kind(1.0d0)
        character(len=*), intent(in) :: filename
        real(dp), dimension(:,:), intent(in) :: array
        integer(dp) :: i, j, rows, cols
        integer(dp) :: unit
    
        ! Determine the dimensions of the array
        rows = size(array, 1)
        cols = size(array, 2)
    
        ! Open the file for writing
        unit = 10
        open(unit=unit, file=filename, status='unknown', action='write')
    
        ! Write the array to the file
        do i = 1, rows
          do j = 1, cols
            write(unit, '(F8.4)', advance='no') array(i, j)
            if (j < cols) then
              write(unit, '(A)', advance='no') ' '  ! Add a space between columns
            end if
          end do
          write(unit, *)  ! New line at the end of each row
        end do
    
        ! Close the file
        close(unit)
      end subroutine write_2d_array_to_file


      subroutine read_data_from_file(filename, data_array, num_elements, status)
        implicit none
        
        ! Input parameters
        character(len=*), intent(in) :: filename  ! File name
        integer, parameter :: dp = kind(1.0d0)
        integer(dp), intent(out) :: num_elements      ! Number of elements read
        real(dp), dimension(:), allocatable, intent(out) :: data_array ! Array to store data
        
        ! Output parameter
        integer(dp), intent(out) :: status             ! Status of the read operation
        
        ! Local variables
        integer(dp) :: unit_number
        integer(dp) :: i
        real(dp) :: temp_value
        character(len=100) :: line
        
        ! Initialize variables
        status = 0
        num_elements = 0
        
        ! Open the file
        open(unit=unit_number, file=filename, status='old', action='read', iostat=status)
        
        if (status /= 0) then
            print *, 'Error opening file ', trim(filename)
            return
        end if
        
        ! Count number of lines in the file
        do
            read(unit_number, '(A)', iostat=status) line
            if (status /= 0) exit ! Exit loop if end of file or error
            num_elements = num_elements + 1
        end do
        
        ! Allocate memory for data_array
        allocate(data_array(num_elements))
        
        ! Rewind the file to read data
        rewind(unit_number)
        
        ! Read the data into the array
        do i = 1, num_elements
            read(unit_number, *) data_array(i)
        end do
        
        ! Close the file
        close(unit_number)
        
    end subroutine read_data_from_file


    subroutine tick(t)
      integer, intent(OUT) :: t
      call system_clock(t)
    end subroutine tick

    ! returns time in seconds from now to time described by t 
    ! Source: https://stackoverflow.com/questions/24395686/best-way-to-write-a-large-array-to-file-in-fortran-text-vs-other
    real function tock(t)
      integer, intent(in) :: t
      integer :: now, clock_rate
      call system_clock(now,clock_rate)
      tock = real(now - t)/real(clock_rate)
    end function tock

End Module utils
