! defs.h
! Parameters 'header' file.

!-------------------------------------!
! Parameters for comparing in k-space
!-------------------------------------!
real :: kmin = 2.0
real :: kmax = 15.0

!------------------------------------!
! Minimal distance between two atoms
!------------------------------------!
real :: DIS_MIN = 1.0

!---------------------------------------------------!
! Initial number of interpolation points in k-space.
!---------------------------------------------------!
integer :: INTERPO_INI
parameter (INTERPO_INI = 400)

!-----------------------!
! Convergence criterion
!-----------------------!
real :: criterion = 0.01

!-------------------------!
! Initial number of atoms
!-------------------------!
integer :: NATOMS_INI
parameter ( NATOMS_INI = 7000 )

!--------------------------------------------!
! Initial number of FEFF calculation points.
!--------------------------------------------!
integer :: POINTS_NUM_INI
parameter ( POINTS_NUM_INI = 1000 )

!--------------------------------------------!
! Initial number of experimental k points.
!--------------------------------------------!
integer :: EXP_POINTS_NUM_INI
parameter ( EXP_POINTS_NUM_INI = 500 )

!---------------------------------!
! Parameters for feff calculation.
!---------------------------------!
real :: RPATH = 10.0
real :: EXAFS_kmax = 20.0
real :: RADIUS = 10.0
real :: AMP=1.0
real :: DMIN=1.0
character(len=2) :: ABS_EDGE = "K"

!-------------------------------------------!
! Initial temperature for DLPOLY simulation
!-------------------------------------------!
! Make sure the temperature given here is 
! consistent with the initial value given in
! DLPOLY 'CONTROL' file 'temperature' card.
!-------------------------------------------!
real :: INIT_TEMP = 300.0

!--------------------------------------------!
! The flag determining whether increaing the 
! temperature during the MD evovling, or not.
! '0' -> DO NOT increase temperature.
! '1' -> Increase temperature for 
!        'artificial' converging purpose.
!--------------------------------------------!
integer :: TEMP_EVOLVE = 0

!------------------------------------------------!
! The specifier for the 'CONFIG' file type.
! '0'   ->      Velocity and force not specified.
! '1'   ->      Velocity and force specified.
!------------------------------------------------!
integer :: CONFIG_TYPE = 0

!-------------------------------------------!
! Number of threads for parallel processing.
!-------------------------------------------!
integer :: NTHREADS = 12

!------------------------------!
! Experimental k-chi file name
!------------------------------!
character(len=13) :: EXP_CHI_FILE = "exp_exafs.chi"

!----------------!
! Log file names
!----------------!
character(len=9) :: error_log = 'error.log'
character(len=10) :: fileio_log = 'fileio.log'
character(len=11) :: running_log = 'running.log'

!------------------------------------------------------------------------------------!
! Table for transformation between atomic symbol and the corresponding atomic number.
!------------------------------------------------------------------------------------!
character(len=10), dimension(73) :: pe_table = (/' H', 'He', 'Li', 'Be', ' P', ' C', &
                ' N', ' O', ' F', 'Ne', 'Na', 'Mg', 'Al', 'Si', ' P', ' S', 'Cl', 'Ar', ' K', &
                'Ca', 'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', &
                'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', ' Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', &
                'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', ' I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', &
                ' W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', &
                'Fr', 'Ra'/)