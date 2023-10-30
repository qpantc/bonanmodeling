program main

    use MathFunctions
    implicit none
    
    ! =============================================================================================
    ! Parameters for beta distribution
    ! =============================================================================================
    real(4),dimension(3) :: p = [2.5, 3.5, 11.5]
    real(4),dimension(3) :: q = [2.5, 2.0, 3.5]
    
    ! =============================================================================================
    ! Canopy parameters
    ! =============================================================================================
    real(4) :: LAI = 5.0
    real(4) :: hc = 10.0

    INTEGER(8) :: layer !层数
    real(8),dimension(:), allocatable :: z
    real(4) :: dz = 0.05
        
    ! =============================================================================================
    ! Loop variables
    ! =============================================================================================
    integer :: i,j ! i: loop for beta parameters; j: loop for scaled height
    real(4) :: z_max, z_min
    real(4) :: x, lad, sum ! x: scaled height from 0 to 1; lad: leaf area density; sum: sum lai
    real(4), dimension(:), allocatable :: y1, y2, y3
    character(80) :: filename

    ! Define the filename for the CSV file
    filename = "leaf_area_density.csv"

    ! =============================================================================================
    ! Define the minimum and maximum heights
    ! =============================================================================================
    z_min = 0
    z_max = hc
    layer = INT((z_max-z_min)/dz)

    ! =============================================================================================
    ! Calculate the number of heights and allocate memory
    ! 指定数组的索引范围，避免歧义，但目前来看，即使明确了数组z范围从1到200layer，但还是可以额对z(0)的位置进行赋值，
    ! 并且不会改变size(z)的结果
    ! =============================================================================================
    allocate(z(1:layer), y1(1:layer), y2(1:layer), y3(1:layer)) 

    ! =============================================================================================
    ! Create a vector of heights
    ! =============================================================================================
    do j=0,layer
        z(j) = z_min + dz * j
    end do
    write(*,*) layer,size(z),z(1),z(200) !fortran 只记录1开始的元素个数
    ! Open the CSV file for writing
    open(unit=10, file=filename, status="replace")

    ! Write header to the CSV file
    write(10, *) "Height, LAD1, LAD2, LAD3"

    ! =============================================================================================
    ! Calculate the leaf area density profile for each (p, q) with a loop
    ! =============================================================================================
    do i = 1, size(p)
        sum = 0.0
        do j = 0, layer
            x = z(j)/hc ! x: scaled height from 0 to 1;
            lad = (LAI / hc) * (x**(p(i) - 1.) * (1. - x)**(q(i) - 1.)) / beta(p(i), q(i))
            ! write(*,*) j,z(j),x,lad,size(z)
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
