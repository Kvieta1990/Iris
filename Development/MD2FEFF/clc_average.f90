subroutine clc_average(NATOMS, &
                       feff_omega_temp, feff_e_temp, feff_k_temp, &
                       feff_mu_temp, feff_mu0_temp, feff_chi_temp, &
                       feff_omega, feff_e, feff_k, feff_mu, feff_mu0, &
                       feff_chi, points_num,iteration)
!---------------------------------------------------------------------------------------!
!
! clc_average.f90
! Calculate the configuration average for all FEFF results where each 
! single atom in the
! system was taken as the absorber, respectively.
! Taking parameters:    Number of atoms                                 -       Input
!                       Input array for 'omega' for all clusters        -       Input
!                       Input array for 'e' for all clusters            -       Input
!                       Input array for 'k' for all clusters            -       Input
!                       Input array for 'mu' for all clusters           -       Input
!                       Input array for 'mu0' for all clusters          -       Input
!                       Input array for 'chi' for all clusters          -       Input
!                       Output array for 'omega' of feff output         -       Output
!                       Output array for 'e' of feff output             -       Output
!                       Output array for 'k' of feff output             -       Output
!                       Output array for 'mu' of feff output            -       Output
!                       Output array for 'mu0' of feff output           -       Output
!                       Output array for 'chi' of feff output           -       Output
!                       Number of calculation points from FEFF          -       Input
!                       Total iteration number                          -       Input
!
!---------------------------------------------------------------------------------------!
!
!---------------------------------------------------------------------------------------!
!
! Author: Yuanpeng Zhang
! School of Physics and Astronomy
! Queen Mary, University of London
! Nov-05-2015
!
!---------------------------------------------------------------------------------------!

implicit none
include "defs.h"

integer, intent(in)                                                     ::      NATOMS, points_num, &
                                                                                iteration
real (KIND=8), intent(in), dimension(NATOMS_INI,POINTS_NUM_INI)         ::      feff_omega_temp, &
                                                                                feff_e_temp, &
                                                                                feff_k_temp, &
                                                                                feff_mu_temp, &
                                                                                feff_mu0_temp, &
                                                                                feff_chi_temp
real (KIND=8), intent(out), dimension(POINTS_NUM_INI)                   ::      feff_omega, &
                                                                                feff_e, feff_k, &
                                                                                feff_mu, feff_mu0, &
                                                                                feff_chi

integer                                         ::      i, j, NATOMS_EFF, io, p_index
integer, dimension(NATOMS_INI)                  ::      flag
character(len=10)                               ::      str_iteration,stri
character(len=20)                               ::      error_file
real (KIND=8)                                   ::      feff_omega_ptemp, feff_e_ptemp, feff_k_ptemp, &
                                                        feff_mu_ptemp, feff_mu0_ptemp, feff_chi_ptemp
logical                                         ::      error_exist

NATOMS_EFF = NATOMS
do i = 1, NATOMS
        flag(i) = 0
        do j = 1, points_num
                if (feff_chi_temp(i,j) > 1) then
                        flag(i) = 1
                        exit
                end if
        end do
        if (flag(i) == 1) then
                NATOMS_EFF = NATOMS_EFF - 1
        end if
end do

do i = 1, NATOMS
        write(stri,'(I0)') i
        error_file = "Cluster_" // trim(stri) // ".error"
        error_file = trim(error_file)
        inquire(file=error_file,exist=error_exist)
        if (error_exist) NATOMS_EFF = NATOMS_EFF - 1
end do

do i = 1, POINTS_NUM_INI
        feff_omega(i) = 0
        feff_e(i) = 0
        feff_k(i) = 0
        feff_mu(i) = 0
        feff_mu0(i) = 0
        feff_chi(i) = 0
end do

! Calculate the configuration average for each quantity.
do i = 1, NATOMS
        do j = 1, points_num
                feff_omega(j) = feff_omega(j) + feff_omega_temp(i,j)
                feff_e(j) = feff_e(j) + feff_e_temp(i,j)
                feff_k(j) = feff_k(j) + feff_k_temp(i,j)
                if (flag(i) == 0) then
                        feff_mu(j) = feff_mu(j) + feff_mu_temp(i,j)
                        feff_mu0(j) = feff_mu0(j) + feff_mu0_temp(i,j)
                        feff_chi(j) = feff_chi(j) + feff_chi_temp(i,j)
                end if
        end do
end do

do j = 1, points_num
        feff_omega(j) = feff_omega(j)/NATOMS
        feff_e(j) = feff_e(j)/NATOMS
        feff_k(j) = feff_k(j)/NATOMS
        feff_mu(j) = feff_mu(j)/NATOMS_EFF
        feff_mu0(j) = feff_mu0(j)/NATOMS_EFF
        feff_chi(j) = feff_chi(j)/NATOMS_EFF
end do

if (iteration > 0) then
        write(str_iteration,'(I0)') iteration-1
        open(unit=41,file="Config_Archive/xmu_conf_aver_" // trim(str_iteration) // ".dat", &
        status="old", action="read")
        do i = 1, 4
                read(41,*)
        end do
        p_index = 1
        do
                read(41,*,IOSTAT=io) feff_omega_ptemp,feff_e_ptemp, &
                feff_k_ptemp, feff_mu_ptemp, feff_mu0_ptemp, &
                feff_chi_ptemp
                if (io > 0) then
                        write(*,*) "Error encountered while reading 'xmu_conf_aver_" // &
                                    trim(str_iteration) // ".dat' file."
                        stop
                else if (io < 0) then
                        exit
                else
                        feff_omega(p_index) = (feff_omega_ptemp * iteration + &
                                               feff_omega(p_index))/(iteration + 1)
                        feff_e(p_index) = (feff_e_ptemp * iteration + &
                                           feff_e(p_index))/(iteration + 1)
                        feff_k(p_index) = (feff_k_ptemp * iteration + &
                                           feff_k(p_index))/(iteration + 1)
                        feff_mu(p_index) = (feff_mu_ptemp * iteration + &
                                            feff_mu(p_index))/(iteration + 1)
                        feff_mu0(p_index) = (feff_mu0_ptemp * iteration + &
                                             feff_mu0(p_index))/(iteration + 1)
                        feff_chi(p_index) = (feff_chi_ptemp * iteration + &
                                             feff_chi(p_index))/(iteration + 1)
                        p_index = p_index + 1
                end if
        end do
        p_index = p_index - 1
        close(41)
        if (p_index /= points_num) then
                write(*,*) "Points number inconsistent!"
                write(*,*) "Check the 'xmu_conf_aver_" // trim(str_iteration) // ".dat' file."
                stop
        end if
end if

end subroutine clc_average