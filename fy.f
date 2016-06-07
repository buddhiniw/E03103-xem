	real*8	function fy(y,f0,bigB,a,b,alpha)


!!!	real*8	function fy(y,y0,b,c,d)		!DHP version
C+______________________________________________________________________________
!
! DESCRIPTION:
!
!   Calculate scaling function F(Y) from a model, which is a sum of a gaussian
!   plus 1/cosh functions. This model is originally due to R. McKeown, although
!   his version gives the wrong asymptotic behavior at large Y. The model is:
!
!   F(y) = (C/y0) * exp[-0.5(y/y0)^2] + D / (exp[B*Y] + exp[-B*Y])
!
!   At large |Y|, the asymptotic behavior is ~ D*exp[-B*|Y|].
!
!   From P. Bosted, et. al. (PRL vol. 49, p. 1380 [1982]), I read the following
!   numbers from their graph of F(y) for He4:
!
!	Y = -0.6 (GeV/c), F(Y) = .001
!	Y = -0.3        , F(Y) = .07
!
!   Solving for D and B with the above points yields D = 4.9000, B=14.162.
!   Parameter C is then determined from the normalization of F. It has the
!   value: C = 0.1821. The parameter Y0 is calculated from an ad hoc expression
!   for the fermi momentum, exactly as in McKeown's subroutine.
!
! USAGE:
!
!   FY is an R*8 function. Its arguments are:
!
!   Y:		- Scaling variable in units of GeV/c.
!   A_NUC:	- "A" of target nucleus.
!
! AUTHOR:
!
!   D. Potterveld, 10/10/86
C-______________________________________________________________________________

        implicit none

!!!	real*8 y,y0,b,c,d
	real*8 y,f0,bigB,a,b,alpha

! Parameter values as originally specified by R. McKeown:

!	parameter b = 11.513
! 	parameter c = .290
! 	parameter d = 2.0

! Parameter values derived from above reference for HE3 data:

!	parameter b = 13.040
! 	parameter c = .2788
! 	parameter d = 2.500

! Parameter values derived from above reference for HE4 data:

!	parameter b = 14.162
! 	parameter c = 0.1821
! 	parameter d = 4.9000

! Parameter values derived from NE3 iron data:

!	parameter b = 11.09679
! 	parameter c = 6.5577E-2
! 	parameter d = 5.9032

! Parameter values derived from fitting new (11/24/86) model to He data
! of 3.6 GeV and 16 deg. The fit adjusted C and Y0, holding B and D fixed.
! This results in an F(y) which is not normalized to 1. (But it fits the
! data very well!)

ccc	parameter b = 11.2303
ccc	parameter c = .1058
ccc	parameter d = 6.381

! Parameters for an obviously wrong F(y) - used to check stability of
! radiative corrections.

!	parameter b = 10.6
! 	parameter c = 0.7
! 	parameter d = 0.024

C ============================ Executable Code =================================

cC Compute Y0. - from old version that passed different quantaties
c
c	pf = (1. - exp(-a_nuc/8.))*.22 + .04
c	y0 = pf/1.638		!Formerly 1.8. The new value comes from
c				!the He fit described above.

C compute F(Y).  B,C,D are PASSED to fy.

!!	fy = C*exp(-0.5*(y/Y0)**2)/Y0 + D/(exp(B*y)+exp(-B*y))

!!! quick test - this is better for 15 degrees.
!!!	fy = (.8*c)*exp(-0.5*(y/Y0)**2)/Y0 + (.8*d)/(exp((.8*b)*y)+exp(-(.8*b)*y))

! Claudio/West's version: test for deuterium

!	f0=.0083
!	bigB=.000213
!	a=.00506
!	b=.006
!	alpha=45

!	f0=.010
!	bigB=1.01214E-03
!	A=7.51329E-03
!	B=9.63269E-03
!	alpha=45

	y=y*1000
	fy = (f0-bigB)*alpha**2*exp(-(a*y)**2)/(alpha**2+y**2)+bigB*exp(-b*abs(y))
	fy = fy*1000
	y=y/1000

	return
	end
