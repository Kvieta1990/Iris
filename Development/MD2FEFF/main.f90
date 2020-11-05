program main
use omp_lib
!--------------------------------------------------------------------------!
!
! main.f90
! The main program carrying out the MD simulation and iterating with FEFF
! calculation until the deviation from experimental data smaller than the
! pre-defined criterion.
! 
!--------------------------------------------------------------------------!
!
!--------------------------------------------------------------------------!
!
! Author: Yuanpeng Zhang
! School of Physics and Astronomy
! Queen Mary, University of London
! Nov-13-2015
!
!--------------------------------------------------------------------------!

implicit none
include "defs.h"

integer                                                 ::      i, id, NATOMS, cluster_id, &
                                                                points_num, exp_points_num, &
                                                                CHUNK, NATOMS_FINAL
integer                                                 ::      iteration=0, temp_step=0
integer                                                 ::      unit_num_rlog = 40, &
                                                                unit_num_cf = 45, &
                                                                unit_num_w_erlog = 50
integer                                                 ::      unit_num_rxmu = 0, &
                                                                unit_num_w_erlog_xmu = 0
integer, dimension(NATOMS_INI)                          ::      points_num_temp
real (KIND=8)                                           ::      deviation, deviation_temp
real, dimension(NATOMS_INI,3)                           ::      Atom_Coor
real (KIND=8), dimension(POINTS_NUM_INI)                ::      feff_omega, feff_e, feff_k, &
                                                                feff_mu, feff_mu0, feff_chi
real (KIND=8), dimension(EXP_POINTS_NUM_INI)            ::      exp_k, exp_chi
real (KIND=8), dimension(NATOMS_INI,POINTS_NUM_INI)     ::      feff_omega_temp, feff_e_temp, &
                                                                feff_k_temp, feff_mu_temp, &
                                                                feff_mu0_temp, feff_chi_temp
character(len=2)                                        ::      Abs_Ele
character(len=10)                                       ::      input_file, str_iteration, stri
character(len=100)                                      ::      title_line
character(len=10), dimension(NATOMS_INI)                ::      Atom_Sym
character(len=20)                                       ::      cluster
character(len=10)                                       ::      str_INIT_TEMP
logical                                                 ::      bak_dir_e

! Read in the element type of the absorption atom.
if (iargc() == 0) then
        write(*,*) "Error! Absorption atom type should be given as the argument!"
        stop
end if
call getarg(1,Abs_Ele)

! Initial run to check the deviation between the theoretical and experimental spectra.
input_file = "CONFIG"
call config_filein(input_file, NATOMS, Atom_Coor, Atom_Sym, title_line)
call system("touch " // running_log)
open(unit=unit_num_rlog, file=running_log, status="old", action="write")
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "! MD configuration file - 'CONFIG' - successfully read in."
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""

call pre_feff(NATOMS, Atom_Coor, Atom_Sym, Abs_Ele, title_line,NATOMS_FINAL)
NATOMS = NATOMS_FINAL
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "! FEFF input files preparation finishes successfully."
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""

call system("touch clc_feff.sh")
open(unit=unit_num_cf, file="clc_feff.sh", status="old", action="write")
write(unit_num_cf,'(A)') "#!/bin/bash"
write(unit_num_cf,'(A)') ""
write(unit_num_cf,'(A)') "temp_dir=$1"
write(unit_num_cf,'(A)') ""
write(unit_num_cf,'(A)') "cd $temp_dir"
write(unit_num_cf,'(A)') "feff feff.inp"
write(unit_num_cf,'(A)') "find . -type f ! -name 'feff.inp' ! -name 'xmu.dat' -exec rm -rf {} \;"
write(unit_num_cf,'(A)') "cd .."
close(unit_num_cf)

CHUNK = NATOMS/NTHREADS/5

call omp_set_dynamic(.FALSE.)
! Do the FEFF calculation in parallel by calculating several different configurations
! at the same time.
!$omp parallel private(i, stri, cluster) num_threads(NTHREADS)
!$omp do schedule(dynamic, CHUNK)
do i = 1, NATOMS
        write(stri,'(I0)') i
        cluster = "Cluster_" // trim(stri)
        write(*,*) trim(cluster)
        call system("sh clc_feff.sh " // trim(cluster))
end do
!$omp end do nowait
!$omp end parallel
write(unit_num_rlog,'(A)') "!----------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "! FEFF simulation for each single configuration successfully finishes."
write(unit_num_rlog,'(A)') "!----------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""

                            !------------------------------------!
                            !  The code block for testing only!  !
                            !------------------------------------!

                            !NATOMS = 6346
                            !CHUNK = NATOMS/NTHREADS/5
                            !title_line="Test line."
                            !call system("touch " // running_log)
                            !open(unit=unit_num_rlog, file=running_log, &
                            !     status="old", action="write")

                            !------------------------------------!
                            !  The code block for testing only!  !
                            !------------------------------------!

! Read in all FEFF output files in parallel.
!$omp parallel &
!$omp private(i,stri,unit_num_rxmu,unit_num_w_erlog_xmu) &
!$omp shared(feff_omega_temp,feff_e_temp,feff_k_temp,feff_mu_temp,feff_mu0_temp,feff_chi_temp) &
!$omp num_threads(NTHREADS)
!$omp do schedule(dynamic, CHUNK)
do i = 1, NATOMS
        stri = ""
        write(stri,'(I0)') i
        unit_num_rxmu = NATOMS + i
        unit_num_w_erlog_xmu = 2 * NATOMS + i
        call read_xmu(i,stri,unit_num_rxmu,unit_num_w_erlog_xmu,NATOMS,&
                      feff_omega_temp,feff_e_temp,feff_k_temp,&
                      feff_mu_temp,feff_mu0_temp,feff_chi_temp,points_num_temp)
end do
!$omp end do
!$omp end parallel

! Check the consistence of the the points number for each cluster.
do i = 2, NATOMS
        if(points_num_temp(i) /= points_num_temp(i-1) .AND. &
                                 points_num_temp(i) /= POINTS_NUM_INI .AND. &
                                 points_num_temp(i-1) /= POINTS_NUM_INI) then
                write(stri,'(I0)') i
                call system("touch " // error_log)
                open(unit=unit_num_w_erlog, file=error_log, status="old", action="write")
                write(unit_num_w_erlog,'(A)') &
                "!------------------------------------------------------------------------------!"
                write(unit_num_w_erlog,'(A)') &
                "! The points number of 'Cluster_" // trim(stri) // "' is not consistent with others."
                write(unit_num_w_erlog,'(A)') &
                "!------------------------------------------------------------------------------!"
                write(unit_num_w_erlog,'(A)') ""
                close(unit_num_w_erlog)
                stop "Error encountered! See the error log file for details."
        end if
end do

points_num = points_num_temp(1)

write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "! All FEFF output successfully read in."
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""

call clc_average(NATOMS, feff_omega_temp, feff_e_temp, feff_k_temp, &
                 feff_mu_temp, feff_mu0_temp, feff_chi_temp, &
                 feff_omega, feff_e, feff_k, feff_mu, feff_mu0, feff_chi, points_num,iteration)
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "! Configuration average calculation successfully finishes."
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""

call ca_fileout(feff_omega, feff_e, feff_k, feff_mu, &
                feff_mu0, feff_chi, points_num, NATOMS, title_line)
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "! Configuration average output to 'xmu_conf_aver.dat'."
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""

call exp_filein(exp_k, exp_chi, exp_points_num)
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "! Experimental spectra successfully read in."
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""

call dev_check(feff_k,feff_chi,points_num,exp_k,exp_chi,exp_points_num,deviation)
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "! Initial deviation checking finishes successfully."
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""

! Make the back up directory, where all the medium-step files will be copied to.
inquire(file="./" // "Config_Archive" // "/.", exist=bak_dir_e)
if (.NOT. bak_dir_e) then
        call system("mkdir Config_Archive")
end if
write(str_iteration,'(I0)') iteration
call system("cp xmu_conf_aver.dat ./Config_Archive/xmu_conf_aver_" // trim(str_iteration) // ".dat")
call system("cp CONFIG ./Config_Archive/CONFIG_" // trim(str_iteration))
write(unit_num_rlog,'(A)') "                        -------------------                       "
write(unit_num_rlog,'(A)') "                        | Initialization  |                       "
write(unit_num_rlog,'(A)') "                        |    finishes     |                       "
write(unit_num_rlog,'(A)') "                        -------------------                       "
write(unit_num_rlog,'(A)') ""
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "! The following is the running log for MD+FEFF simulation cycles."
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A,I0,A)') "! There are ", NATOMS, " atoms in the system."
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""
write(unit_num_rlog,'(A,I0)') "     Iteration - ", iteration
write(unit_num_rlog,'(A)') ""
write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
write(unit_num_rlog,'(A,I0,A)') "          ! The current iteration number: ", iteration, "."
write(unit_num_rlog,'(A,F15.5)') "          ! The initial deviation is: ", deviation
write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""
write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "          ! The configuration averaged k-chi data has been moved to:"
write(unit_num_rlog,'(3A)') "          ! /Config_Archive/xmu_conf_aver_" // trim(str_iteration) // ".dat"
write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "          ! The MD configuration file was copied to:"
write(unit_num_rlog,'(2A)') "          ! /Config_Archive/CONFIG_" // trim(str_iteration)
write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""

write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "          ! Step to the cycle involving MD evolution and FEFF simulation."
write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""

call system("/data/home/apw247/Packages/DL_POLY4/execute/DLPOLY.Z")
call system("cp REVCON ./Config_Archive/REVCON_" // trim(str_iteration))
call system("cp REVIVE ./Config_Archive/REVIVE_" // trim(str_iteration))
call system("cp STATIS ./Config_Archive/STATIS_" // trim(str_iteration))
call system("cp OUTPUT ./Config_Archive/OUTPUT_" // trim(str_iteration))
write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "          ! All MD simulation output files (OUTPUT, REVCON, REVIVE, STATIS)"
write(unit_num_rlog,'(A)') "          ! were copied to:"
write(unit_num_rlog,'(A)') "          ! /Config_Archive/*_" // trim(str_iteration)
write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""

do while(deviation > criterion)
        iteration = iteration + 1
        write(str_iteration,'(I0)') iteration

        input_file = "REVCON"
        call config_filein(input_file, NATOMS, Atom_Coor, Atom_Sym, title_line)

        call pre_feff(NATOMS, Atom_Coor, Atom_Sym, Abs_Ele, title_line, NATOMS_FINAL)
        NATOMS = NATOMS_FINAL

        CHUNK = NATOMS/NTHREADS/5
        call omp_set_dynamic(.FALSE.)
        !$omp parallel private(i, stri) num_threads(NTHREADS)
        !$omp do schedule(dynamic, CHUNK)
        do i = 1, NATOMS
                write(stri,'(I0)') i
                cluster = "Cluster_" // trim(stri)
                write(*,*) trim(cluster)
                call system("sh clc_feff.sh " // trim(cluster))
        end do
        !$omp end do nowait
        !$omp end parallel

        !$omp parallel &
        !$omp private(i,stri,unit_num_rxmu,unit_num_w_erlog_xmu) &
        !$omp shared(feff_omega_temp,feff_e_temp,feff_k_temp,feff_mu_temp,feff_mu0_temp,feff_chi_temp) &
        !$omp num_threads(NTHREADS)
        !$omp do schedule(dynamic, CHUNK)
        do i = 1, NATOMS
                stri = ""
                write(stri,'(I0)') i
                unit_num_rxmu = NATOMS + i
                unit_num_w_erlog_xmu = 2 * NATOMS + i
                call read_xmu(i,stri,unit_num_rxmu,unit_num_w_erlog_xmu,NATOMS,&
                              feff_omega_temp,feff_e_temp,feff_k_temp,&
                              feff_mu_temp,feff_mu0_temp,feff_chi_temp,points_num_temp)
        end do
        !$omp end do
        !$omp end parallel
        
        do i = 2, NATOMS
                if (points_num_temp(i) /= points_num_temp(i-1) .AND. &
                                          points_num_temp(i) /= POINTS_NUM_INI .AND. &
                                          points_num_temp(i-1) /= POINTS_NUM_INI) then
                        write(stri,'(I0)') i
                        call system("touch " // error_log)
                        open(unit=unit_num_w_erlog, file=error_log, status="old", action="write")
                        write(unit_num_w_erlog,'(A)') &
                        "!------------------------------------------------------------------------------!"
                        write(unit_num_w_erlog,'(A)') &
                        "! The points number of 'Cluster_" // trim(stri) // "' is not consistent with others."
                        write(unit_num_w_erlog,'(A)') &
                        "!------------------------------------------------------------------------------!"
                        write(unit_num_w_erlog,'(A)') ""
                        close(unit_num_w_erlog)
                        stop "Error encountered! See the error log file for details."
                end if
        end do
        
        points_num = points_num_temp(1)

        call clc_average(NATOMS, &
                  feff_omega_temp, feff_e_temp, feff_k_temp, &
                  feff_mu_temp, feff_mu0_temp, feff_chi_temp, &
                  feff_omega, feff_e, feff_k, feff_mu, &
                  feff_mu0, feff_chi, points_num,iteration)
        write(unit_num_rlog,'(A,I0)') "     Iteration - ", iteration
        write(unit_num_rlog,'(A)') ""
        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
        write(unit_num_rlog,'(A)') "          ! Configuration average calculation successfully finishes."
        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
        write(unit_num_rlog,'(A)') ""

        call ca_fileout(feff_omega, feff_e, feff_k, feff_mu, &
                        feff_mu0, feff_chi, points_num, NATOMS, title_line)
        
        call dev_check(feff_k,feff_chi,points_num,exp_k,exp_chi,exp_points_num,deviation_temp)

        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
        write(unit_num_rlog,'(A,I0,A)') "          ! The current iteration number: ", iteration, "."
        write(unit_num_rlog,'(A,F15.5)') "          ! The current deviation is: ", deviation_temp
        write(unit_num_rlog,'(A,F15.5)') "          ! The pre-defined criterion is: ", criterion
        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
        write(unit_num_rlog,'(A)') ""
        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
        write(unit_num_rlog,'(A)') "          ! The configuration averaged k-chi data has been moved to:"
        write(unit_num_rlog,'(3A)') "          ! /Config_Archive/xmu_conf_aver_" // trim(str_iteration) // ".dat"
        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
        write(unit_num_rlog,'(A)') "          ! The MD configuration file was copied to:"
        write(unit_num_rlog,'(2A)') "          ! /Config_Archive/CONFIG_" // trim(str_iteration)
        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
        write(unit_num_rlog,'(A)') ""
        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
        write(unit_num_rlog,'(A)') "          ! All MD simulation output files (OUTPUT, REVCON, REVIVE, STATIS)"
        write(unit_num_rlog,'(A)') "          ! were copied to:"
        write(unit_num_rlog,'(A)') "          ! /Config_Archive/*_" // trim(str_iteration)
        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
        write(unit_num_rlog,'(A)') ""

        if (TEMP_EVOLVE == 0) then
                if (iteration > 100) then
                        write(unit_num_rlog,'(A)') &
                        "          !------------------------------------------------------------------!"
                        write(unit_num_rlog,'(A)') "          ! Iteration number uplimite reached! Program existing!"
                        write(unit_num_rlog,'(A)') &
                        "          !------------------------------------------------------------------!"
                        write(unit_num_rlog,'(A)') ""
                        exit
                end if
                deviation = deviation_temp
                call system("mv REVCON CONFIG")
                write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
                write(unit_num_rlog,'(A)') "          ! 'REVCON' moved to 'CONFIG' -> Move to next round of MD evolution." 
                write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
                write(unit_num_rlog,'(A)') ""

                call system("/data/home/apw247/Packages/DL_POLY4/execute/DLPOLY.Z")

                call system("cp xmu_conf_aver.dat ./Config_Archive/xmu_conf_aver_" // trim(str_iteration) // ".dat")
                call system("cp CONFIG ./Config_Archive/CONFIG_" // trim(str_iteration))
                call system("cp REVCON ./Config_Archive/REVCON_" // trim(str_iteration))
                call system("cp REVIVE ./Config_Archive/REVIVE_" // trim(str_iteration))
                call system("cp STATIS ./Config_Archive/STATIS_" // trim(str_iteration))
                call system("cp OUTPUT ./Config_Archive/OUTPUT_" // trim(str_iteration))

                write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
                write(unit_num_rlog,'(A)') "          ! Current MD evolution finished! -> Back to check the deviation."
                write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
                write(unit_num_rlog,'(A)') ""
        else if (TEMP_EVOLVE == 1) then
                if (deviation_temp <= deviation) then
                        deviation = deviation_temp
                        call system("mv REVCON CONFIG")
                        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
                        write(unit_num_rlog,'(A)') "          ! Current deviation smaller than that from last step!"
                        write(unit_num_rlog,'(A)') "          ! MD evolution accepted!"
                        write(unit_num_rlog,'(A)') "          ! 'REVCON' moved to 'CONFIG' -> Move to next round of MD evolution." 
                        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
                        write(unit_num_rlog,'(A)') ""

                        call system("/data/home/apw247/Packages/DL_POLY4/execute/DLPOLY.Z")

                        call system("cp xmu_conf_aver.dat ./Config_Archive/xmu_conf_aver_" // trim(str_iteration) // ".dat")
                        call system("cp CONFIG ./Config_Archive/CONFIG_" // trim(str_iteration))
                        call system("cp REVCON ./Config_Archive/REVCON_" // trim(str_iteration))
                        call system("cp REVIVE ./Config_Archive/REVIVE_" // trim(str_iteration))
                        call system("cp STATIS ./Config_Archive/STATIS_" // trim(str_iteration))
                        call system("cp OUTPUT ./Config_Archive/OUTPUT_" // trim(str_iteration))
                else
                        if (temp_step >= 20) then
                                write(unit_num_rlog,'(A)') &
                                "------------------------------------------------------------------"
                                write(unit_num_rlog,'(A)') "          ! Temperature uplimit (Initial temperature + 100) reached!"
                                write(unit_num_rlog,'(A)') "          ! Program existing!"
                                write(unit_num_rlog,'(A)') &
                                "------------------------------------------------------------------"
                                exit
                        else
                                temp_step = temp_step + 1
                        end if
                        INIT_TEMP = INIT_TEMP + 5
                        write(str_INIT_TEMP,'(F0.1)') INIT_TEMP
                        call system("sed -i -e '/temperature/c temperature          " // trim(str_INIT_TEMP) // "' CONTROL")

                        call system("/data/home/apw247/Packages/DL_POLY4/execute/DLPOLY.Z")

                        call system("cp xmu_conf_aver.dat ./Config_Archive/xmu_conf_aver_" // trim(str_iteration) // ".dat")
                        call system("cp CONFIG ./Config_Archive/CONFIG_" // trim(str_iteration))
                        call system("cp REVCON ./Config_Archive/REVCON_" // trim(str_iteration))
                        call system("cp REVIVE ./Config_Archive/REVIVE_" // trim(str_iteration))
                        call system("cp STATIS ./Config_Archive/STATIS_" // trim(str_iteration))
                        call system("cp OUTPUT ./Config_Archive/OUTPUT_" // trim(str_iteration))

                        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
                        write(unit_num_rlog,'(A)') "          ! Current deviation is larger than thar from last step."
                        write(unit_num_rlog,'(A)') "          ! MD evolution rejected!"
                        write(unit_num_rlog,'(A)') "          ! The MD simulation temperature was increased by 5 K."
                        write(unit_num_rlog,'(A)') "          ! The 'CONFIG' file from last step not changed!"
                        write(unit_num_rlog,'(A)') "          ! Move to next round of MD simulation."
                        write(unit_num_rlog,'(A)') "          !------------------------------------------------------------------!"
                        write(unit_num_rlog,'(A)') ""
                end if
        else
                call system("touch " // error_log)
                open(unit=unit_num_w_erlog, file=error_log, status="old", action="write")
                write(unit_num_w_erlog,'(A)') &
                "!------------------------------------------------------------------------------!"
                write(unit_num_w_erlog,'(A)') "! Wrong 'TEMP_EVOLVE' flag given in 'defs.h' file!"
                write(unit_num_w_erlog,'(A)') "! The 'TEMP_EVOLE' value should be either '0' or '1'!"
                write(unit_num_w_erlog,'(A)') "! '0' -> No temperature increasing."
                write(unit_num_w_erlog,'(A)') "! '1' -> Increase temperature if deviation not decreasing."
                write(unit_num_w_erlog,'(A)') &
                "!------------------------------------------------------------------------------!"
                write(unit_num_w_erlog,'(A)') ""
                close(unit_num_w_erlog)
                stop "Error encountered! See the error log file for details."
        end if 
end do

call system("rm -rf clc_feff.sh")
call system("rm -rf fort*")
write(unit_num_rlog,'(A)') "!---------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "! The deviation smaller than the pre-defined criterion."
write(unit_num_rlog,'(A)') "! Exiting the MD+FEFF cycle..."
write(unit_num_rlog,'(A)') "!---------------------------------------------------------------!"
write(unit_num_rlog,'(A)') ""
write(unit_num_rlog,'(A)') "                        -------------------                       "
write(unit_num_rlog,'(A)') "                        |    Converged!   |                       "
write(unit_num_rlog,'(A)') "                        -------------------                       "
write(unit_num_rlog,'(A)') ""
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------------------------!"
write(unit_num_rlog,'(A,F15.5)') "! The final deviation (=SQRT((exp1-feff1)^2+(exp2-feff2)^2+...)) is : ", deviation
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------------------------!"
write(unit_num_rlog,'(A)') "! The final configuration and the FEFF configuration average results are stored in"
write(unit_num_rlog,'(A)') "! 'xmu_conf_aver.dat' and 'CONFIG' file, respectively."
write(unit_num_rlog,'(A)') "!------------------------------------------------------------------------------------!"
close(unit_num_rlog)

end program main
