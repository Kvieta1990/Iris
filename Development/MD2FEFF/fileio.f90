subroutine config_filein(input_file,NATOMS,Atom_Coor,Atom_Sym,title_line)
!-----------------------------------------------------------------------------------------------!
!
! fileio.f90
! Read in the MD CONFIG/REVCON file and extract the atomic coordinates and
! symbols.
! Taking parameters:    File name                                               -       Input
!                       Number of atoms                                         -       Output
!                       Output array for atomic coordinates                     -       Output
!                       Output array for atomic symbol                          -       Output
!                       Array to accept the DL_POLY CONFIG/REVCON title line    -       Output
!
!-----------------------------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------------------------------!
!
! Author: Yuanpeng Zhang
! School of Physics and Astronomy
! Queen Mary, University of London
! Nov-04-2015
!
!-----------------------------------------------------------------------------------------------!

implicit none
include "defs.h"

integer, intent(out)                                    ::      NATOMS
real, intent(out), dimension(NATOMS_INI,3)              ::      Atom_Coor
character(len=10), intent(in)                           ::      input_file
character(len=100), intent(out)                         ::      title_line
character(len=10), intent(out), dimension(NATOMS_INI)   ::      Atom_Sym

integer        ::      unit_num_r = 10, unit_num_revr = 11, unit_num_erlog = 15
integer        ::      i, io, atom_index
logical        ::      file_exists

! Read in the MD initial 'CONFIG' file.
if (input_file == "CONFIG") then
        open(unit=unit_num_r, file=input_file, status='old', action='read')
        read(unit_num_r,'(A)') title_line
        do i = 1,4
                read(unit_num_r,*)
        end do
        atom_index = 1
        do
                read(unit_num_r,*,IOSTAT=io) Atom_Sym(atom_index),NATOMS
                read(unit_num_r,*,IOSTAT=io) Atom_Coor(atom_index,1),Atom_Coor(atom_index,2),Atom_Coor(atom_index,3)
                if (CONFIG_TYPE == 1) then
                        read(unit_num_r,*,IOSTAT=io)
                        read(unit_num_r,*,IOSTAT=io)
                else if (CONFIG_TYPE /= 0) then
                        call system("touch " // error_log)
                        open(unit=unit_num_erlog, file=error_log, status='old', action='write')
                        write(unit_num_erlog,'(A)') "!---------------------------------------------------------!"
                        write(unit_num_erlog,'(A)') "!   The 'CONFIG' file type should be either '0' or '1'!   !"
                        write(unit_num_erlog,'(A)') "!    '0'      ->      Velocity and force not specified.   !"
                        write(unit_num_erlog,'(A)') "!     '1'      ->      Velocity and force specified.      !"
                        write(unit_num_erlog,'(A)') "!---------------------------------------------------------!"
                        write(unit_num_erlog,'(A)') ""
                        close(unit_num_erlog)
                        stop 'Error encountered! Check the error log file for possible reason!'
                end if
                if (io > 0) then
                        call system("touch " // error_log)
                        open(unit=unit_num_erlog, file=error_log, status='old', action='write')
                        write(unit_num_erlog,'(A)') "!---------------------------------------------------------!"
                        write(unit_num_erlog,'(A)') "!   Error encountered while reading the 'CONFIG' file!    !"
                        write(unit_num_erlog,'(A)') "!           Please check the 'CONFIG' file!               !"
                        write(unit_num_erlog,'(A)') "!---------------------------------------------------------!"
                        write(unit_num_erlog,'(A)') ""
                        close(unit_num_erlog)
                        stop 'Error encountered! Check the error log file for possible reason!'
                else if (io < 0) then
                        exit
                else
                        atom_index = atom_index + 1
                end if
        end do
        atom_index = atom_index - 1
        close(unit_num_r)
        print *, Atom_Coor(NATOMS,1),Atom_Coor(NATOMS,2),Atom_Coor(NATOMS,3)
! Read in the MD "REVCON" file - the medium configuration file during the simulation running.
else if (input_file == "REVCON") then
        open(unit=unit_num_revr, file=input_file, status='old', action='read')
        read(unit_num_revr,'(A)') title_line
        do i = 1,4
                read(unit_num_revr,*)
        end do
        atom_index = 1
        do
                read(unit_num_revr,*,IOSTAT=io) Atom_Sym(atom_index),NATOMS
                read(unit_num_revr,*,IOSTAT=io) Atom_Coor(atom_index,1),Atom_Coor(atom_index,2),Atom_Coor(atom_index,3)
                read(unit_num_revr,*,IOSTAT=io)
                read(unit_num_revr,*,IOSTAT=io)
                if (io > 0) then
                        call system("touch " // error_log)
                        open(unit=unit_num_erlog, file=error_log, status='old', action='write')
                        write(unit_num_erlog,'(A)') "!---------------------------------------------------------!"
                        write(unit_num_erlog,'(A)') "!   Error encountered while reading the 'REVCON' file!    !"
                        write(unit_num_erlog,'(A)') "!           Please check the 'REVCON' file!               !"
                        write(unit_num_erlog,'(A)') "!---------------------------------------------------------!"
                        write(unit_num_erlog,'(A)') ""
                        close(unit_num_erlog)
                        stop 'Error encountered! Check the error log file for possible reason!'
                else if (io < 0) then
                        exit
                else
                        atom_index = atom_index + 1
                end if
        end do
        close(unit_num_revr)
else
       call system("touch " // error_log)
       open(unit=unit_num_erlog, file=error_log, status='old', action='write')
       write(unit_num_erlog,'(A)') "!----------------------------------------------------------------!"
       write(unit_num_erlog,'(A)') "!   Error encountered while reading the 'REVCON/CONFIG' file!    !"
       write(unit_num_erlog,'(A)') "!          Please check the 'REVCON' or 'CONFIG' file!           !"
       write(unit_num_erlog,'(A)') "!----------------------------------------------------------------!"
       write(unit_num_erlog,'(A)') ""
       close(unit_num_erlog)
       stop 'Error encountered! Check the error log file for possible reason!'
end if

end subroutine config_filein

subroutine ca_fileout(feff_omega,feff_e,feff_k,feff_mu,feff_mu0,feff_chi,points_num, NATOMS, title_line)
!-----------------------------------------------------------------------------!
!
! fileio.f90
! Output the configuration average FEFF simulation result to file.
! Taking parameters:    Array for 'omega'                       -       Input
!                       Array for 'e'                           -       Input
!                       Array for 'k'                           -       Input
!                       Array for 'mu'                          -       Input
!                       Array for 'mu0'                         -       Input
!                       Array for 'chi'                         -       Input
!                       Number of points in FEFF output         -       Input
!                       Number of atoms                         -       Input
!                       Title line in MD CONFIG/REVCON file     -       Input
!
!-----------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------------!
!
! Author: Yuanpeng Zhang
! School of Physics and Astronomy
! Queen Mary, University of London
! Nov-11-2015
!
!-----------------------------------------------------------------------------!

implicit none
include "defs.h"

integer, intent(in)                                     ::      points_num, NATOMS
real (KIND=8), intent(in), dimension(POINTS_NUM_INI)    ::      feff_omega, feff_e, &
                                                                feff_k, feff_mu, feff_mu0, feff_chi
character(len=100), intent(in)                          ::      title_line

integer        ::       i
integer        ::       unit_num_aw = 25

call system("touch xmu_conf_aver.dat")
open(unit=unit_num_aw, file="xmu_conf_aver.dat", status="old", action="write")
write(unit_num_aw,'(A)') "# Configuration averaged FEFF simulation for: ", title_line
write(unit_num_aw,'(A)') "#  -----------------------------------------------------------------------"
write(unit_num_aw,'(A)') "#  omega    e    k    mu    mu0     chi     @#"

do i = 1, points_num
        write(unit_num_aw,'(F15.3,F15.3,F15.3,E15.5,E15.5,E15.5)') & 
        feff_omega(i), feff_e(i), feff_k(i), feff_mu(i), feff_mu0(i), feff_chi(i)
end do

close(unit_num_aw)

end subroutine ca_fileout

subroutine exp_filein(exp_k,exp_chi,exp_points_num)
!-----------------------------------------------------------------------------!
!
! fileio.f90
! Read in the experimental EXAFS spectra in k-space
! Taking parameters:    Array for exp 'k' values        -       Output
!                       Array for exp 'chi' values      -       Output
!                       Number of exp k points          -       Output
!
!-----------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------------!
!
! Author: Yuanpeng Zhang
! School of Physics and Astronomy
! Queen Mary, University of London
! Nov-12-2015
!
!-----------------------------------------------------------------------------!

implicit none
include "defs.h"

integer, intent(out)                                            ::      exp_points_num
real (KIND=8), intent(out), dimension(EXP_POINTS_NUM_INI)       ::      exp_k
real (KIND=8), intent(out), dimension(EXP_POINTS_NUM_INI)       ::      exp_chi

integer        ::       unit_num_r = 30, unit_num_erlog = 35
integer        ::       io

open(unit=unit_num_r, file=EXP_CHI_FILE, status="old", action="read")

exp_points_num = 1
do
         read(unit_num_r,*,IOSTAT=io) exp_k(exp_points_num), exp_chi(exp_points_num)
         if (io > 0) then
                 call system("touch " // error_log)
                 open(unit=unit_num_erlog, file=error_log, status="old", action="write")
                 write(unit_num_erlog,'(A)') "!-------------------------------------------------------------!"
                 write(unit_num_erlog,'(A)') "! Error encountered while reading the experimental k-chi data !"
                 write(unit_num_erlog,'(A)') "!          Please check the experimental data file!           !"
                 write(unit_num_erlog,'(A)') "!-------------------------------------------------------------!"
                 write(unit_num_erlog,'(A)') ""
                 close(unit_num_erlog)
                 stop "Error encountered! Check the error log file for possible reasons!"
         else if (io < 0) then
                 exit
         else
                 exp_points_num = exp_points_num + 1
         end if
end do

exp_points_num = exp_points_num - 1

close(unit_num_r)

end subroutine exp_filein