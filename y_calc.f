	subroutine y_calc(e1,e2,theta,matom,m_rec,es,y,y_calc_ok)
C+______________________________________________________________________________
!
! Y_CALC - Calculate Y-scaling variable from double precision passed arguments:
!
!	e1:	(R*8) - Incident energy in GeV
!	e2:	(R*8) - Scattered energy in GeV
!	theta:	(R*8) - Scattering angle in radians
!	m1:	(R*8) - Mass of target nucleus in GeV/c^2.
!	m_rec:	(R*8) - Mass of recoiling nucleon in GeV/c^2.
!	es:	(R*8) - Separation energy in GeV.
!	y:	(R*8) - Calculated value of Y-scaling variable, in GeV/c.
!
! The function returns a value of .TRUE. if it could calculate Y, or a value
! of .FALSE. if there is no real solution for Y.
!
! The equation for Y assumes that the recoiling (A-1) system has a mass (m2)
! given by:
!
!	m2 = m1 + es - m_rec
!
! where m_p is the proton mass.
!
C-______________________________________________________________________________

	implicit none

        logical y_calc_ok
	real*8 es
	real*8 e1,e2,theta,m1,m_rec,y,matom
	real*8 m2,w,th,q4_2,q2,q
	real*8 pp2
	real*8 temp,coeff_a,coeff_b,coeff_c,root2

	include 'math_physics.inc'

C ============================ Executable Code =================================

	y_calc_ok = .false.			!assume failure

	y=0.0
	m1 = matom*m_amu

c	write(6,*) 'inputs:'
c	write(6,*) 'ebeam',e1
c	write(6,*) 'eprime',e2
c	write(6,*) 'theta',theta*r_d
c	write(6,*) 'mtar',m1
c	write(6,*) 'esep',es
c	write(6,*) 'y_init',y

C Compute kinematics
	m2 = m1 + es - m_rec			!Mass of (A-1) system
	w = e1 - e2	    			!Energy loss
c	th = theta*d_r				!Theta in radians
	th=theta
	q4_2 = 4.*e1*e2*(sin(th/2.))**2.	!4-momentum transfer squared
	q2 = q4_2 + w*w				!3-momentum    "        "
	q  = sqrt(q2)

c	write(6,*) 'nu,q2,q4_2',w,q2,q4_2
C Compute approximate value for k_perp.

! K_perp now suppressed since two parameters M_rec and E_sep make it
! entirely irrelevent.

C$$$	a = m1/m_p				!Approximate value for A.
C$$$	pf = .22*(1.-exp(-a/8.)) + .04
C$$$	pp2 = pf**2.

	pp2 = 0.				!Suppress k_perp.

C Compute terms in quad. formula.

	temp = q2 + m_rec**2 - m2*m2 - (m1 + w)**2
	coeff_a = 4.*(q2 - (m1 + w)**2.)
	coeff_b = 4.*q*temp
	coeff_c = temp**2. - 4.*(m2*m2 + pp2)*(m1 + w)**2.

C If no real solution, return ERROR.

	root2 = coeff_b**2. - 4.*coeff_a*coeff_c
	if (root2.lt.0..or.coeff_a.eq.0.) return

C Otherwise, return Y, and SUCCESS.

	y = (-coeff_b - sqrt(root2))/(2.*coeff_a)

	y_calc_ok = .true.

c	write(6,*) 'y final',y

	return
	end
