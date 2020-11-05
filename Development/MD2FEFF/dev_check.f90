subroutine dev_check(feff_k,feff_chi,points_num,exp_k,exp_chi,exp_points_num,deviation)
!----------------------------------------------------------------------------------------!
!
! dev_check.f90
! Calculate the quantitative consistence between the experimental k-chi data
! and the configuration averaged data from FEFF simulation.
! Taking parameters:    Array for k points from FEFF                    -       Input
!                       Array for 'chi' from FEFF                       -       Input
!                       Number of k points in FEFF output               -       Input
!                       Array for experimental k points                 -       Input
!                       Array for experimental 'chi' data               -       Input
!                       Number of k points in experimental data         -       Input
!                       Output the deviation                            -       Output
!
!----------------------------------------------------------------------------------------!
!
!----------------------------------------------------------------------------------------!
! 
! Author: Yuanpeng Zhang
! School of Physics and Astronomy
! Queen Mary, University of London
! Nov-12-2015
!
!----------------------------------------------------------------------------------------!

implicit none
include "defs.h"

integer, intent(in)                                             ::      points_num, exp_points_num
real (KIND=8), intent(in), dimension(POINTS_NUM_INI)            ::      feff_k, feff_chi
real (KIND=8), intent(in), dimension(EXP_POINTS_NUM_INI)        ::      exp_k, exp_chi
real (KIND=8), intent(out)                                      ::      deviation

integer                                         ::      i
real (KIND=8), dimension(points_num)            ::      feff_k_ini
real (KIND=8), dimension(1,points_num)          ::      feff_chi_ini
real (KIND=8), dimension(exp_points_num)        ::      exp_k_ini
real (KIND=8), dimension(1,exp_points_num)      ::      exp_chi_ini
real (KIND=8), dimension(INTERPO_INI)           ::      k_interpo
real (KIND=8), dimension(1,INTERPO_INI)         ::      feff_chi_interpo, exp_chi_interpo

do i = 1, points_num
        feff_k_ini(i) = feff_k(i)
        feff_chi_ini(1,i) = feff_chi(i)
end do

do i =1, exp_points_num
        exp_k_ini(i) = exp_k(i)
        exp_chi_ini(1,i) = exp_chi(i)
end do

deviation = 0

! Make the standard k-grid, on which the interpolation will be targeted.
do i = 1, INTERPO_INI
        k_interpo(i) = kmin + (i - 1)*((kmax - kmin)/(INTERPO_INI - 1))
end do


! Interpolate the experimental and simulated results to the same k-grid.
! Previously, the linear and Lagrange interpolation method was used, which
! turned out to not work properly in current situation.
call interp_nearest(1, points_num, feff_k_ini, feff_chi_ini, &
                   INTERPO_INI, k_interpo, feff_chi_interpo)
call interp_nearest(1, exp_points_num, exp_k_ini, exp_chi_ini, & 
                   INTERPO_INI, k_interpo, exp_chi_interpo)

! Calculate the deviation.
do i = 1, INTERPO_INI
        deviation = deviation + (feff_chi_interpo(1,i) - exp_chi_interpo(1,i))**2
end do

deviation = sqrt(deviation)

end subroutine dev_check