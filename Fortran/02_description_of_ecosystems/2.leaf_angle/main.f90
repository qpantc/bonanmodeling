
! =============================================================================================
! 使用数组保存结果的方式
! =============================================================================================
program main

    ! =============================================================================================
    ! Evaluate the leaf angle probability density function (PDF) from the beta
    ! distribution using the mean and standard deviation of the leaf inclination 
    ! angle. The leaf angle PDF is calculated for 9 angle classes between 5 and 85 
    ! degrees in increments of 10 degrees.
    ! 
    ! PDF <- mean and standard deviation of the leaf inclination angle. 倾斜角度
    ! lad_ave - mean leaf angle (radians)
    ! lad_std - standard deviation of leaf angle (radians)
    ! =============================================================================================

    use MathFunctions
    implicit none

    ! =============================================================================================
    ! the constant will be used in this case
    ! =============================================================================================
    real(8), parameter :: pi = 3.141592653
    real(8), parameter :: dangle = 5.0 * (pi/180.0)! 只有设置为参数之后才可以在变量计算中随意使用
    ! real(8) :: dangle = 5.0 ! 这样定义的话，无法在定义变量时使用
    integer, parameter :: N = int(90.* (pi/180.0)/dangle)

    ! =============================================================================================
    ! output file
    ! =============================================================================================
    character(80) :: filename

    ! =============================================================================================
    ! lean name leaf angle mean and standard deviation
    ! =============================================================================================
    character(20) :: lad_type
    real(8), dimension(5) :: lad_ave = [26.76,63.5,45.0,45.00,57.3]*(pi/180)
    real(8), dimension(5) :: lad_std = [18.51,18.5,16.3,25.98,21.55]*(pi/180)

    ! =============================================================================================
    ! arrar for angle probability density function
    ! =============================================================================================
    real(8), dimension(N) :: angle,beta_pdf,beta_lad,exact_pdf,exact_lad

    ! =============================================================================================
    ! p,q parameters for the beta distribution
    ! =============================================================================================
    real(8) :: num,den
    real(8) :: p,q

    ! =============================================================================================
    ! parameters for the loop
    ! i loop for angle; j loop for input leaf value; x: scaled angle value
    ! =============================================================================================
    integer :: i,j
    real(8) :: x
    ! =============================================================================================
    ! 18 5-degree bins for calculating the leaf inclination angle probability density function (PDF) 
    ! and fractional abundance (lad) 
    ! =============================================================================================
    do i = 1, size(angle)
        angle(i) = dangle * i
        write(*,*) angle(i)
    end do

    ! =============================================================================================
    ! prepare file to save the outputs
    ! =============================================================================================
    filename = "leaf_angle_distribution.csv"
    ! Open the CSV file for writing
    open(unit=10, file=filename, status="replace")
    ! Write header to the CSV file
    write(10, *) "Name, lad_ave, lad_std, angle, beta_pdf, beta_lad, exact_pdf, exact_lad"

    ! =============================================================================================
    ! convert leaf angle data to the p,q parameters
    ! =============================================================================================
    do j = 1, size(lad_ave)
        num = 1 - (lad_std(j)*lad_std(j) + lad_ave(j)*lad_ave(j)) / (lad_ave(j) * pi / 2)
        den = (lad_std(j)*lad_std(j) + lad_ave(j)*lad_ave(j)) / (lad_ave(j)*lad_ave(j)) - 1
        p = num / den
        q = ((pi/2) / lad_ave(j) - 1) * p

        write(*,"('p,q value for led_ave =',F8.4,' led_std =',F8.4,' is:',2F8.4, ' num and den:',2F8.4)")  &
            lad_ave(j),lad_std(j),p,q,num,den

        ! =============================================================================================
        ! loop through each angle
        ! =============================================================================================
        do i = 1, size(angle)
            x = angle(i)/(pi/2.)

            beta_pdf(i) = 2./pi * x**(p - 1.) * (1. - x)**(q - 1.) / beta(p,q)   ! Leaf angle probability density function 
            beta_lad(i) = beta_pdf(i) * dangle          ! Fraction of leaves in this angle bin
            ! write(*,*) "The radians value of angle is: ", angle(i),x, beta_pdf(i),dangle

            ! using different way to calculate ['Planophile','Erectophile','Plagiophile','Uniform','Spherical']
            if ( j == 1) then
                lad_type = 'Planophile'
                exact_pdf(i) = 2 / pi * (1 + cos(2 * angle(i)))
            end if 
            if ( j == 2) then
                lad_type = 'Erectophile'
                exact_pdf(i) = 2 / pi * (1 - cos(2 * angle(i)))
            end if 

            if ( j == 3) then
                lad_type = 'Plagiophile'
                exact_pdf(i) = 2 / pi * (1 - cos(4 * angle(i)))
            end if 

            if ( j == 4) then
                lad_type = 'Uniform'
                exact_pdf(i) = 2 / pi
            end if 

            if ( j == 5) then
                lad_type = 'Spherical'
                exact_pdf(i) = sin(angle(i))
            end if

            exact_lad(i) = exact_pdf(i) * dangle
            write(10, "(A12,',',6(F8.4,','),(F8.4))") &
                lad_type,lad_std(j),lad_ave(j),angle(i), &
                beta_pdf(i), beta_lad(i), exact_pdf(i), exact_lad(i)

        end do

    end do


end program main

! program LeafAnglePDF

!     use MathFunctions
!   implicit none

!   character(12) :: leaf
!   real(8) :: lad_ave, lad_std
!   real(8) :: num, den, p, q
!   real(8) :: dangle
!   real(8), dimension(9) :: angle, beta_pdf, beta_lad, exact_pdf, exact_lad
!   real(8) :: beta_sum, beta_ave, exact_sum, exact_ave
!   real(8) :: analytical_ave
!   real(8) :: F1, F2, F3, beta_xl, exact_xl
!   integer :: i
!   real(8) :: fp, fq ! Declare fp and fq as real variables


!   ! --- The variable "leaf" defines the leaf angle distribution type with parameters:
!   ! lad_ave - mean leaf angle (radians)
!   ! lad_std - standard deviation of leaf angle (radians)

!   leaf = 'Planophile'
!   ! leaf = 'Erectophile'
!   ! leaf = 'Plagiophile'
!   ! leaf = 'Uniform'
!   ! leaf = 'Spherical'

!   select case (trim(leaf))
!   case ('Planophile')
!     lad_ave = 26.76 * (pi/180)
!     lad_std = 18.5068 * (pi/180)
!   case ('Erectophile')
!     lad_ave = 63.24 * (pi/180)
!     lad_std = 18.4960 * (pi/180)
!   case ('Plagiophile')
!     lad_ave = 45.00 * (pi/180)
!     lad_std = 16.2681 * (pi/180)
!   case ('Uniform')
!     lad_ave = 45.00 * (pi/180)
!     lad_std = 25.9808 * (pi/180)
!   case ('Spherical')
!     lad_ave = 57.30 * (pi/180)
!     lad_std = 21.5485 * (pi/180)
!   end select

!   ! --- Convert these to the p,q parameters for the beta distribution

!   num = 1 - (lad_std*lad_std + lad_ave*lad_ave) / (lad_ave * pi / 2)
!   den = (lad_std*lad_std + lad_ave*lad_ave) / (lad_ave*lad_ave) - 1
!   p = num / den
!   q = ((pi/2) / lad_ave - 1) * p

!   ! --- Calculate leaf inclination angle probability density function (PDF) and fractional abundance (lad) for 9 10-degree bins

!   dangle = 10 * (pi/180) ! Leaf inclination angle increment (radians)
!   angle = [5, 15, 25, 35, 45, 55, 65, 75, 85] * (pi/180) ! Leaf inclination angle (degrees -> radians)

!   ! --- Loop through each angle

!   do i = 1, 9
!      x = angle(i) / (pi/2)
!      fp = x ** (p - 1)
!      fq = (1 - x) ** (q - 1)
!      beta_pdf(i) = 2 / pi * fp * fq / beta(p, q) ! Leaf angle probability density function
!      beta_lad(i) = beta_pdf(i) * dangle ! Fraction of leaves in this angle bin
!   end do

!   ! --- Calculate the known solution

!   do i = 1, 9
!      select case (trim(leaf))
!      case ('Planophile')
!         exact_pdf(i) = 2 / pi * (1 + cos(2 * angle(i)))
!      case ('Erectophile')
!         exact_pdf(i) = 2 / pi * (1 - cos(2 * angle(i)))
!      case ('Plagiophile')
!         exact_pdf(i) = 2 / pi * (1 - cos(4 * angle(i)))
!      case ('Uniform')
!         exact_pdf(i) = 2 / pi
!      case ('Spherical')
!         exact_pdf(i) = sin(angle(i))
!      end select

!      ! Exact relative leaf angle distribution (fraction)
!      exact_lad(i) = exact_pdf(i) * dangle
!   end do

!   ! --- Print out fractional abundance and compare with known solution

!   beta_sum = 0
!   beta_ave = 0
!   exact_sum = 0
!   exact_ave = 0

!   write(*, *)
!   write(*, *) 'Leaf type =', leaf
!   write(*, *) '     Angle         beta            exact'
!   do i = 1, 9
!      beta_sum = beta_sum + beta_lad(i)
!      beta_ave = beta_ave + angle(i) * beta_lad(i)
!      exact_sum = exact_sum + exact_lad(i)
!      exact_ave = exact_ave + angle(i) * exact_lad(i)
!      write(*, '(3(F10.2, 1X))') angle(i)*180/pi, beta_lad(i), exact_lad(i)
!   end do

!   write(*, *)
!   write(*, *) 'beta distribution'
!   write(*, *) 'Sum of leaf angle distribution =', beta_sum
!   write(*, *) 'Mean leaf angle =', beta_ave*180/pi
!   write(*, *)
!   write(*, *) 'Exact solution'
!   write(*, *) 'Sum of leaf angle distribution =', exact_sum
!   write(*, *) 'Mean leaf angle =', exact_ave*180/pi

!   ! --- Analytical mean leaf angle

!   select case (trim(leaf))
!   case ('Planophile')
!      analytical_ave = integral(0, pi/2, fx)
!   case ('Erectophile')
!      analytical_ave = integral(0, pi/2, fx)
!   case ('Plagiophile')
!      analytical_ave = integral(0, pi/2, fx)
!   case ('Uniform')
!      analytical_ave = integral(0, pi/2, fx)
!   case ('Spherical')
!      analytical_ave = integral(0, pi/2, fx)
!   end select

!   write(*, *)
!   write(*, *) 'Analytical solution'
!   write(*, *) 'Mean leaf angle =', analytical_ave*180/pi

!   ! --- Calculate Ross index

!   F1 = beta_lad(1) + beta_lad(2) + beta_lad(3)
!   F2 = beta_lad(4) + beta_lad(5) + beta_lad(6)
!   F3 = beta_lad(7) + beta_lad(8) + beta_lad(9)
!   beta_xl = 0.5 * (abs(0.134-F1) + abs(0.366-F2) + abs(0.5-F3))
!   if ((0.5-F3) < 0) beta_xl = -beta_xl

!   F1 = exact_lad(1) + exact_lad(2) + exact_lad(3)
!   F2 = exact_lad(4) + exact_lad(5) + exact_lad(6)
!   F3 = exact_lad(7) + exact_lad(8) + exact_lad(9)
!   exact_xl = 0.5 * (abs(0.134-F1) + abs(0.366-F2) + abs(0.5-F3))
!   if ((0.5-F3) < 0) exact_xl = -exact_xl

!   write(*, *)
!   write(*, *) 'Ross index'
!   write(*, *) 'beta distribution =', beta_xl
!   write(*, *) '   Exact solution =', exact_xl

!   ! --- Graph PDFs

!   do i = 0, 100
!      x = i * (pi/2) / 100
!      select case (trim(leaf))
!      case ('Planophile')
!         y = 2 / pi * (1 + cos(2 * x))
!      case ('Erectophile')
!         y = 2 / pi * (1 - cos(2 * x))
!      case ('Plagiophile')
!         y = 2 / pi * (1 - cos(4 * x))
!      case ('Uniform')
!         y = 2 / pi
!      case ('Spherical')
!         y = sin(x)
!      end select
!      x = x * 180/pi ! radians -> degrees
!      write(*, *) x, beta_pdf(i), exact_pdf(i), y
!   end do

! end program LeafAnglePDF


