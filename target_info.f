	subroutine target_info(z,a,aux)
C+______________________________________________________________________________
!
! TARGET_INFO - Calculate material thickness before and after scattering.
!
! ARGUMENTS:
!
!   THETA_TAR:  (R*8):  Target angle - 0 deg is normal to beam, increase CCW.
!			(i.e. increasing angle faces SOS).
!   THETA:	(R*8):	Scattering angle in Degrees (Input).
!   TB:		(R*8):	Thickness BEFORE scatter, in radiation lengths.
!   TA:		(R*8):	Thickness AFTER  scatter, in radiation lengths.
!   Z:		(R*8):	'Z' of target nucleus.
!   A:		(R*8):	'A' of target nucleus.
!   tgt_len:    (R*8):  target length in cm.
!   M_TGT:	(R*8):	Mass in GeV/c^2 of target nucleus.
!   AUX:	(R*8):	Array of length 10 containing:
!			aux(1) = RESOL parameter for DEEPSIG smearing (not used).
!			aux(2) = Separation energy in GeV of target nucleus.
!			aux(3) = F(0)  parameter for F(y) function
!			aux(4) = BigB  parameter for F(y)
!			aux(5) = a     parameter for F(y)
!			aux(6) = b     parameter for F(y)
!			aux(7) = alpha parameter for F(y)
!
! AUX definitions for old f(y) model:
!		!	aux(1) = Mass of recoiling nucleon in GeV/c^2.
!		!	aux(2) = Separation energy in GeV of target nucleus.
!		!	aux(3) = Y0 parameter for F(y) function
!		!	aux(4) = B  parameter for F(y)
!		!	aux(5) = C  parameter for F(y)
!		!	aux(6) = D  parameter for F(y)
!		!	aux(7) = RESOL parameter for DEEPSIG smearing.
!
C-______________________________________________________________________________

	implicit none

C Get various constants.

	include 'math_physics.inc'

	integer i,j
	logical found
C Declare arguments.

	real*8 z,a,aux(7)

C Local declarations.

C Lookup table from which we find the radiation length and AUX parameters.
C DAVE G NOTE: Stuff for He3,Be,Cu mostly placeholders

	real*8		lookup(10,10)/			!data vs. Z.

!-------------------------------------------------------------------------------
!        H       2H       3He   4He     Be      C       Al       Fe  Cu      Au
!-------------------------------------------------------------------------------
     >   1.,     1.,     2.,     2.,     4.,     6.,    13.,    26.,   29.,     79., !Z's
     >   1.,     2.,     3.,     4.,     9.,    12.,    27.,    56.,   64.,    197., !A's
     >61.28,  122.6,  65.27,  94.32,  65.19,   42.7,  24.01,  13.84, 12.86,   6.461, !R.L.s in g/cm^2.
     > 0.00,   0.14,   0.10,   0.16,   0.10,   0.25,   0.25,   0.25,  0.10,    0.25, !RESOL's
     > 0.00,0.00225,0.00549,0.02020,0.00928,0.01727, 0.0099,0.01060,0.00855, 0.00693, !E_SEP's in GeV.
     > 0.00,0.00914,0.00628,0.00445,0.00391,0.00358,0.00361,0.00367,0.00367,0.00367, !f(0)'s
     > 0.00,0.00025,0.00033,0.00065,0.00060,0.00057,0.00060,0.00062,0.00062,0.00062, !BigB's
     > 0.00,0.00610,0.00710,0.00680,0.00574,0.00510,0.00493,0.00460,0.00460,0.00460, !a's
     > 0.00, 0.0060, 0.0060, 0.0060, 0.0060, 0.0060, 0.0060, 0.0060, 0.0060, 0.0060, !b's
     > 0.00,    45.,    83.,   167.,   166.,   166.,   156.,   138.,   138.,   138./ !alpha's

!!-------------------------------------------------------------------------------
!!        H       2H       He        C       Al       Fe       Au
!!-------------------------------------------------------------------------------
!     >   1.,      1.,      2.,      6.,     13.,     26.,     79., !Z's
!     >   1.,      2.,      4.,     12.,     26.,     56.,    197., !A's
!     >61.28,   122.6,   94.32,    42.7,   24.01,   13.84,   6.461, !R.L.s in g/cm^2.
!     > 0.00,    0.14,    0.16,    0.25,    0.25,    0.25,    0.25, !RESOL's
!     > 0.00, 0.00225,    0.00, 0.01727,  0.0099, 0.01060, 0.00693, !E_SEP's in GeV.
!     > 0.00,   0.010,   0.010,  0.0040,  0.0034,  0.0034,  0.0033, !f(0)'s
!     > 0.00, 0.00101, 0.00101, 0.00079, 0.00053, 0.00053, 0.00083, !BigB's
!     > 0.00, 0.00751, 0.00751, 0.00410,  0.0040,  0.0040, 0.00412, !a's
!     > 0.00, 0.00963, 0.00963,  0.0100,  0.0100,  0.0100,  0.0105, !b's
!     > 0.00,     45.,     45.,    140.,    140.,    140.,    140./ !alpha's

!!-------------------------------------------------------------------------------
!!       H     2H      He      C       Al      Fe      Au
!!-------------------------------------------------------------------------------
!     >   1.,     1.,     2.,     6.,    13.,    26.,    79.,  !Z's
!     >   1.,     2.,     4.,    12.,    26.,    56.,   197.,  !A's
!     >61.28,  122.6,  94.32,   42.7,  24.01,  13.84,  6.461,  !R.L.s in g/cm^2.
!     >  m_p,    .88,    .88,  .9383,  .9383,  .9383,  .9383,  !Recoil masses (GeV/c^2)
!     > 0.00, 0.0022,   0.00,  0.030,  0.040,  0.041,  0.049,  !E_SEP's in GeV.
!     > 0.00,  .0322, .07727, 0.1288, 0.1541, 0.1586, 0.1587,  !Y0's
!     > 0.00,11.2302,11.2302, 12.905, 12.652, 12.652, 12.019,  !B's
!     > 0.00,  .1058,  .1058,  .1131,  .1131,  .1131,  .1131,  !C's
!     > 0.00,  6.381,  6.381,  8.247,  8.085,  8.085,  7.681,  !D's
!     > 0.00,   0.14,   0.16,    0.2,    0.2,    0.2,    0.2/  !RESOL's

C ============================ Executable Code =================================

C Get target specific stuff from lookup table.

c	rho = tgt_thick/tgt_len
	found = .false.
	do i = 1,10				!loop over known targets.
c	  if (lookup(i,1).eq.z.and.dble(int(a)).eq.lookup(i,2)) then
	  if (lookup(i,1).eq.z.and.a.eq.lookup(i,2)) then
	    found = .true.
c	    itar=i
c	    tgt_rl = tgt_thick/lookup(i,3)	!Divide by radiation length.
	    do j = 1,7
	      aux(j) = lookup(i,j+3)
	    enddo
	  endif
	enddo
	if (.not.found) then
	  write(6,*) 'cant find target in lookup table!: z,a',z,a
          return			!Quit if couldn't find info.
	endif

	return	

1001	format(a)
	end
