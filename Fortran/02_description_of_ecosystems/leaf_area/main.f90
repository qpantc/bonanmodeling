program main

    use MathFunctions
    implicit none
    
    ! =============================================================================================
    ! Parameters for beta distribution
    ! =============================================================================================
    real(8),dimension(3) :: p = [2.5, 3.5, 11.5]
    real(8),dimension(3) :: q = [2.5, 2.0, 3.5]
    
    ! =============================================================================================
    ! Canopy parameters
    ! =============================================================================================
    real(8) :: LAI = 5.0
    real(8) :: hc = 10.0

    real(8),dimension(:), allocatable :: z
    real(8) :: dz = 0.1
        
    ! =============================================================================================
    ! Loop variables
    ! =============================================================================================
    integer :: i,j ! i: loop for beta parameters; j: loop for scaled height
    real(8) :: z_max, z_min
    real(8) :: x, lad, sum ! x: scaled height from 0 to 1; lad: leaf area density; sum: sum lai
    real(8), dimension(:), allocatable :: y1, y2, y3
    character(80) :: filename

    ! Define the filename for the CSV file
    filename = "leaf_area_density.csv"

    ! =============================================================================================
    ! Define the minimum and maximum heights
    ! =============================================================================================
    z_min = 0
    z_max = hc

    ! =============================================================================================
    ! Calculate the number of heights and allocate memory
    ! =============================================================================================
    allocate(z(INT((z_max-z_min)/dz)))
    allocate(y1(size(z)), y2(size(z)), y3(size(z)))

    ! =============================================================================================
    ! Create a vector of heights
    ! =============================================================================================
    do j=0,size(z)
        z(j) = z_min + dz * j
    end do

    ! Open the CSV file for writing
    open(unit=10, file=filename, status="replace")

    ! Write header to the CSV file
    write(10, *) "Height, LAD1, LAD2, LAD3"

    ! =============================================================================================
    ! Calculate the leaf area density profile for each (p, q) with a loop
    ! =============================================================================================
    do i = 1, size(p)
        sum = 0.0
        do j = 0, size(z)
            x = z(j)/hc ! x: scaled height from 0 to 1;
            lad = (LAI / hc) * (x**(p(i) - 1.) * (1. - x)**(q(i) - 1.)) / beta(p(i), q(i))
            ! write(*,*) j,x,lad,size(z)
            ! Numerically sum leaf area for each height
            sum = sum + lad * dz

            ! Save output for graphing for each beta distribution
            if (i == 1) then
            y1(j) = lad
            elseif (i == 2) then
            y2(j) = lad
            elseif (i == 3) then
            y3(j) = lad
            end if
        end do

        write(*, '(A,F6.2,F6.2)') "p, q = ", p(i), q(i)
        write(*, '(A,F8.4)') "Leaf area index (numerical) = ", sum
    end do
    
    ! Write the results to the CSV file
    do j = 0, SIZE(z)
        write(10, '(3(F8.4,","),(F8.4))') z(j), y1(j), y2(j), y3(j)
    end do
    
    ! Close the CSV file
    close(10)

    deallocate(z, y1, y2, y3)

end program main
