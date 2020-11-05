subroutine read_xmu(i,stri,unit_num_rxmu,unit_num_w_erlog_xmu,NATOMS,&
           feff_omega_temp, feff_e_temp, feff_k_temp, &
           feff_mu_temp,feff_mu0_temp,feff_chi_temp,points_num_temp)
use omp_lib
!--------------------------------------------------------------------!
!
! xmu_read.f90
! The subroutine for reading in the 'xmu.dat' file corresponding to 
! specified cluster.
! Taking parameters:    Index for cluster       -       Input
!                       Cluster index in
!                       in character form       -       Input
!                       Unit number for 
!                       FEFF output file        -       Input
!                       Unit number for
!                       error log file          -       Input
!                       Number of atoms         -       Input
!                       Array containing        
!                       all the columns of 
!                       'xmu.dat' file for
!                       all clusters            -       Output
!                       Array containing
!                       number of points
!                       for each cluster        -       Output
!
!--------------------------------------------------------------------!
!
!--------------------------------------------------------------------!
!
! Author: Yuanpeng Zhang
! School of Physics and Astronomy
! Queen Mary, University of London
! Nov-20-2015
!
!--------------------------------------------------------------------!

implicit none
include "defs.h"

integer, intent(in)                                                     ::      i, NATOMS, &
                                                                                unit_num_rxmu, &
                                                                                unit_num_w_erlog_xmu
integer, intent(out), dimension(NATOMS_INI)                             ::      points_num_temp
real (KIND=8), intent(out), dimension(NATOMS_INI,POINTS_NUM_INI)        ::      feff_omega_temp, &
                                                                                feff_e_temp, &
                                                                                feff_k_temp, &
                                                                                feff_mu_temp, &
                                                                                feff_mu0_temp, &
                                                                                feff_chi_temp
character(len=10), intent(in)                                           ::      stri

integer                 ::      j, io
character(len=25)       ::      xmu_file
character(len=150)      ::      comm_line
logical                 ::      file_exist

j = 1
xmu_file = "Cluster_" // trim(stri) // "/xmu.dat"

inquire(file=xmu_file,exist=file_exist)
if (.NOT. file_exist) then
        go to 100
end if

open(unit=unit_num_rxmu, file=trim(xmu_file), status="old", action="read")

do
        read(unit_num_rxmu,'(A)',IOSTAT=io) comm_line
        if (io > 0) then
                call system("touch " // error_log // "." // trim(stri))
                open(unit=unit_num_w_erlog_xmu,file=error_log // "." // trim(stri),status="old",action="write")
                write(unit_num_w_erlog_xmu,'(A)') "!----------------------------------------------------------!"
                write(unit_num_w_erlog_xmu,'(A)') "! Error encountered while reading " // trim(xmu_file) // "!"
                write(unit_num_w_erlog_xmu,'(A)') "!----------------------------------------------------------!"
                write(unit_num_w_erlog_xmu,'(A)') ""
                close(unit_num_w_erlog_xmu)
                write(*,*) "Error encountered! See '" // error_log // "." // trim(stri) // "' for details."
                stop 
        end if
        if (comm_line(4:8) == "omega") exit
end do

do
        read(unit_num_rxmu,*,IOSTAT=io) &
                feff_omega_temp(i,j), feff_e_temp(i,j), feff_k_temp(i,j), &
                feff_mu_temp(i,j), feff_mu0_temp(i,j), feff_chi_temp(i,j)
        if (io > 0) then
                call system("touch " // error_log // "." // trim(stri))
                open(unit=unit_num_w_erlog_xmu,file=error_log // "." // trim(stri),status="old",action="write")
                write(unit_num_w_erlog_xmu,'(A)') "!----------------------------------------------------------!"
                write(unit_num_w_erlog_xmu,'(A)') "! Error encountered while reading " // trim(xmu_file) // "!"
                write(unit_num_w_erlog_xmu,'(A)') "!----------------------------------------------------------!"
                write(unit_num_w_erlog_xmu,'(A)') ""
                close(unit_num_w_erlog_xmu)
                write(*,*) "Error encountered! See '" // error_log // "." // trim(stri) // "' for details."
                stop 
        else if (io < 0) then
                exit
        else
                if (feff_k_temp(i,j) /= 0.000) j = j + 1
        end if
end do
points_num_temp(i) = j - 1

close(unit_num_rxmu)
go to 200

100     call system("touch Cluster_" // trim(stri) // ".error")
        do j = 1, POINTS_NUM_INI
                feff_omega_temp(i,j) = 0
                feff_e_temp(i,j) = 0
                feff_k_temp(i,j) = 0
                feff_mu_temp(i,j) = 0
                feff_mu0_temp(i,j) = 0
                feff_chi_temp(i,j) = 0
        end do
        points_num_temp(i) = POINTS_NUM_INI

200     end subroutine read_xmu