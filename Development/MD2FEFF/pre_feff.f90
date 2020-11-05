subroutine pre_feff(NATOMS,Atom_Coor,Atom_Sym,Abs_Ele,title_line,NATOMS_FINAL)
!-------------------------------------------------------------------------------------!
!
! pre_feff.f90
! Taking the atomic coordinates and symbols to prepare the 'feff.inp' input
! file for feff simulation.
! Taking parameters:    Number of atoms                                 -       Input
!                       Input array containing atomic coordinates       -       Input
!                       Input array containing atomic symbols           -       Input
!                       Element type of the absorption atom             -       Input
!                       String containing MD CONFIG/REVCON title line   -       Input
!                       Final number of atoms                           -       Output
!
!-------------------------------------------------------------------------------------!
!
!-------------------------------------------------------------------------------------!
!
! Author: Yuanpeng Zhang
! School of Physics and Astronomy
! Queen Mary, University of London
! Nov-05-2015
!
!-------------------------------------------------------------------------------------!

implicit none
include "defs.h"

integer, intent(in)                                     ::      NATOMS
real, intent(in), dimension(NATOMS_INI,3)               ::      Atom_Coor
character(len=10), intent(in), dimension(NATOMS)        ::      Atom_Sym
character(len=2),intent(in)                             ::      Abs_Ele
character(len=100), intent(in)                          ::      title_line
integer, intent(out)                                    ::      NATOMS_FINAL

integer                                         ::      i, j, k, l, abs_anum, &
                                                        nei_num, anum_temp, nei_uniq, &
                                                        nei_num_final
integer                                         ::      unit_num_w = 55, unit_num_erlog = 60
integer, dimension(NATOMS)                      ::      nei_anum
real                                            ::      dist, distance
real, dimension(NATOMS,3)                       ::      nei_coor
character(len=10), dimension(73)                ::      nei_anum_RMDups
character(len=10), dimension(NATOMS_INI)        ::      str_nei_anum, nei_sym
character(len=10)                               ::      str_abs_anum
character(len=10)                               ::      stri, strj, strk, strx, stry, strz
character(len=15)                               ::      cluster
logical                                         ::      dir_e, file_e, flag, subflag

NATOMS_FINAL = 0
do i = 1, NATOMS
        if (Atom_Sym(i) == Abs_Ele) then
                NATOMS_FINAL = NATOMS_FINAL + 1
                flag = .FALSE.
                nei_num = 0
        
                write(stri,'(I0)') NATOMS_FINAL
                cluster = "Cluster" // "_" // stri
        
                inquire(file="./" // trim(cluster) // "/.", exist=dir_e)
                if (.NOT. dir_e) then
                        call system("mkdir " // trim(cluster))
                end if
        
                ! Most of the FEFF simualtion parameters are given in the 'defs.h' header file.
                ! Therefore if further changing of the simulation parameters is needed, one can
                ! then revise the 'defs.h' file.
                call system("touch " // trim(cluster) // "/feff.inp")
                open(unit=unit_num_w, file=trim(cluster) // "/feff.inp", status='old', action='write')
        
                write(unit_num_w,'(A)') "TITLE: FEFF for " // title_line
                write(unit_num_w,'(A10,A)') "EDGE      ", ABS_EDGE
                write(unit_num_w,'(A10,F3.1,A)') "S02       ", AMP, new_line('A')
                write(unit_num_w,'(A10,A)') "CONTROL   ", "1 1 1 1 1 1"
                write(unit_num_w,'(A10,A)') "PRINT     ", "0 0 0 0 0 0"
                write(unit_num_w,'(A)') ""
                write(unit_num_w,'(A10,F4.1)') "EXAFS     ", EXAFS_kmax
                write(unit_num_w,'(A10,F4.1,A)') "RPATH     ", RPATH, new_line('A')
                write(unit_num_w,*) "POTENTIALS"
                write(unit_num_w,'(4A10)') "*         ", "ipot      ", "z         ", "label     "
        
                ! Check the atomic number of the absorber for current configuration
                do j = 1, size(pe_table)
                        if (Atom_Sym(i) == pe_table(j)) then
                                abs_anum = j
                                if (abs_anum > 56) abs_anum = abs_anum + 15
                                flag = .TRUE.
                                exit
                        end if
                end do
        
                if (.NOT. flag) then
                        call system("touch " // error_log)
                        open(unit=unit_num_erlog, file=error_log, status='old', action='write')
                        write(unit_num_erlog,'(A)') "!----------------------------------------------------&
                                                 -----------------------------------!"
                        write(unit_num_erlog,'(A,A4,A)') "! The atomic symbol for atom number-", stri, &
                                                         " CANNOT be found in the elements list!          !"
                        write(unit_num_erlog,'(A)') "!            Please check the MD CONFIG/REVCON file &
                                                 for possible errors!                !"
                        write(unit_num_erlog,'(A)') "!----------------------------------------------------&
                                                 -----------------------------------!"
                        write(unit_num_erlog,'(A)') ""
                        print *,"Error encountered! Detailed info is given in the 'error.log' file."
                        close(unit_num_erlog)
                        stop 'Program stopped due to error!'
                else
                        write(str_abs_anum,'(I0)') abs_anum
                        write(unit_num_w,'(5A15)') "0         ", adjustl(str_abs_anum), &
                                                   adjustl(pe_table(abs_anum)), "-1        ", "-1        "
                end if
        
                ! Create the neighbour atoms list for current configuration, around the absorber.
                do j = 1, NATOMS
                        subflag = .FALSE.
                        dist = distance(Atom_Coor(i,1),Atom_Coor(j,1),Atom_Coor(i,2),Atom_Coor(j,2),Atom_Coor(i,3),Atom_Coor(j,3))
                        if (dist < RADIUS .AND. dist > DMIN) then
                                nei_num = nei_num + 1
                                nei_coor(nei_num,1) = Atom_Coor(j,1)
                                nei_coor(nei_num,2) = Atom_Coor(j,2)
                                nei_coor(nei_num,3) = Atom_Coor(j,3)
                                nei_sym(nei_num) = Atom_Sym(j)
        
                                do k = 1, size(pe_table)
                                        if (Atom_Sym(j) == pe_table(k)) then
                                                nei_anum(nei_num) = k
                                                if (k > 56) nei_anum(nei_num) = nei_anum(nei_num) + 15
                                                subflag = .TRUE.
                                                exit
                                        end if
                                end do
                                if (.NOT. subflag) then
                                        call system("touch " // error_log)
                                        open(unit=unit_num_erlog, file=error_log, status='old', action='write')
                                        write(unit_num_erlog,'(A)') "!--------------------------------------------&
                                                                 -------------------------------------------!"
                                        write(unit_num_erlog,'(A,A4,A)') "! The atomic symbol for atom number-", &
                                                                         stri, " CANNOT be found in the elements list!          !"
                                        write(unit_num_erlog,'(A)') "!            Please check the MD CONFIG/REVCON &
                                                                 file for possible errors!                !"
                                        write(unit_num_erlog,'(A)') "!--------------------------------------------&
                                                                 -------------------------------------------!"
                                        write(unit_num_erlog,'(A)') new_line('A')
                                        print *,"Error encountered! Detailed info is given in the 'error.log' file."
                                        close(unit_num_erlog)
                                        stop 'Program stopped due to error!'
                                end if
                        end if
                end do
        
                nei_num_final = nei_num
                ! Exclude one of the two atoms that are too close to each other.
                do j = 1, nei_num
                        if (nei_coor(j,1) == 1000 .AND. nei_coor(j,2) == 1000 .AND. nei_coor(j,3) == 1000) then
                                nei_num_final = j - 1
                                exit
                        end if
                        do k = j+1, nei_num
                                dist = distance(nei_coor(j,1),nei_coor(k,1),nei_coor(j,2),nei_coor(k,2),&
                                        nei_coor(j,3),nei_coor(k,3))
                                do while (dist < DIS_MIN)
                                        do l = k, nei_num-1
                                                nei_coor(l,1) = nei_coor(l+1,1)
                                                nei_coor(l,2) = nei_coor(l+1,2)
                                                nei_coor(l,3) = nei_coor(l+1,3)
                                                nei_sym(l) = nei_sym(l+1)
                                                nei_anum(l) = nei_anum(l+1)
                                        end do
                                        nei_coor(nei_num,1) = 1000
                                        nei_coor(nei_num,2) = 1000
                                        nei_coor(nei_num,3) = 1000
                                        dist = distance(nei_coor(j,1),nei_coor(k,1),nei_coor(j,2),nei_coor(k,2),&
                                                nei_coor(j,3),nei_coor(k,3))
                                end do
                        end do
                 end do
                
                ! Obtain the unique member in the neighbour atoms list.
                do j = 1, nei_num_final
                        write(str_nei_anum(j),'(I0)') nei_anum(j)
                end do
                call remove_dups(nei_num_final,str_nei_anum, nei_anum_RMDups, nei_uniq)
        
                ! Assign each unique neighbouring atom with independent potential index.
                do j = 1, nei_uniq
                        write(strj,'(I0)') j
                        read(nei_anum_RMDups(j),'(I10)') anum_temp
                        if (anum_temp > 56) anum_temp = anum_temp - 15
                        write(unit_num_w,'(5A15)') adjustl(strj), adjustl(nei_anum_RMDups(j)), &
                                                   adjustl(pe_table(anum_temp)), "-1        ", "-1        "
                end do
        
                ! Write the atomic coordinates to the 'feff.inp' input file for FEFF.
                write(unit_num_w,'(A)') ""
                write(unit_num_w,'(A)') "ATOMS"
                write(unit_num_w,'(A)') "*         x         y         z         ipot      atom      distance"
                write(strx,'(F10.6)') Atom_Coor(i,1)
                write(stry,'(F10.6)') Atom_Coor(i,2)
                write(strz,'(F10.6)') Atom_Coor(i,3)
                write(unit_num_w,'(5A15)') adjustl(strx), adjustl(stry), adjustl(strz), &
                                           "0         ", adjustl(pe_table(abs_anum))
                do j = 1, nei_num_final
                        do k = 1, nei_uniq
                                read(nei_anum_RMDups(k),'(I10)') anum_temp
                                if (anum_temp > 56) anum_temp = anum_temp - 15
                                if (nei_sym(j) == pe_table(anum_temp)) then
                                        write(strx,'(F10.6)') nei_coor(j,1)
                                        write(stry,'(F10.6)') nei_coor(j,2)
                                        write(strz,'(F10.6)') nei_coor(j,3)
                                        write(strk,'(I0)') k
                                        write(unit_num_w,'(5A15)') adjustl(strx), adjustl(stry), &
                                                                   adjustl(strz), adjustl(strk), adjustl(pe_table(anum_temp))
                                end if
                        end do
                end do
                write(unit_num_w,'(A)') ""
                write(unit_num_w,'(A)') "END"
                close(unit_num_w)
        end if
end do

end subroutine pre_feff

real function distance(x1,x2,y1,y2,z1,z2)
!---------------------------------------------------------------------
!
! function distance()
! The function to calculate the distance between two given points.
! Taking parameters:    'x', 'y' and 'z' coordinates for each of the
!                       two points.
!
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!
! Author: Yuanpeng Zhang
! School of Physics and Astronomy
! Queen Mary, University of London
! Nov-16-2015
!
!---------------------------------------------------------------------

        implicit none

        real, intent(in) :: x1,x2,y1,y2,z1,z2

        distance = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

end function distance