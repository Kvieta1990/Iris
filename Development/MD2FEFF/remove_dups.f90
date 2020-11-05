subroutine remove_dups(nei_num_final,input_array,output_array,uni_num)
!--------------------------------------------------------------------------------------!
!
! remove_dups.f90
! The subroutine to remove the duplicate elements in the input array. The output
! array will be assigned to the 'output_array' passed to current subroutine from
! the program that is calling the current subroutine.
! Taking parameters:    Number of effective elements                    -       Input
!                       Input array to process                          -       Input
!                       Output array to accept the processed results.   -       Output
!			Number of the unique elements.                  -       Output
!
!--------------------------------------------------------------------------------------!
!
!--------------------------------------------------------------------------------------!
!
! Author: Yuanpeng Zhang
! School of Physics and Astronomy
! Queen Mary, University of London
! Nov-05-2015
!
!--------------------------------------------------------------------------------------!
!
!--------------------------------------------------------------------------------------!
!
!         *********************************************************************
!         *    The main block of current code is from the following link:     *
!         *   http://rosettacode.org/wiki/Remove_duplicate_elements#Fortran   *
!         *                     Copyright @ GNU FDL                           *
!         *********************************************************************
!
!--------------------------------------------------------------------------------------!

implicit none
include "defs.h"

integer, intent(out)                                    ::      nei_num_final, uni_num
character(len=10), intent(in), dimension(NATOMS_INI)    ::      input_array
character(len=10), intent(out), dimension(73)           ::      output_array

integer        ::       i, j

!-------------------------------------------------------------------------!
! The following statements does not apply for current program, and it's
! just some tips about the initial value assignment for a character array.
!-------------------------------------------------------------------------!
! NOTICE that when assigning initial value to the character array all at 
! once, all the elements should be with the same length. Here is the link
! to refer to:
! http://stackoverflow.com/questions/29697314/different-character-lengths
! -3-4-in-array-constructor-how-to-trim-strings-fo
!-------------------------------------------------------------------------!
uni_num = 1
output_array(1) = input_array(1)
outer: do i=2,nei_num_final
                do j=1,uni_num
                        if (output_array(j) == input_array(i)) then
                        ! Found a match so start looking again
                                cycle outer
                        end if
                end do
                        ! No match found so add it to the output
                uni_num = uni_num + 1
                output_array(uni_num) = input_array(i)
       end do outer

end subroutine remove_dups