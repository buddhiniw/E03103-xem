C+______________________________________________________________________________
! Math and Physics constants and conversion factors.
C-______________________________________________________________________________

	real*8 pi,d_r,r_d,alpha,m_p,m_e,m_amu,hc_2

	parameter (pi = 3.1415926535D0)
	parameter (d_r = pi/180.)		!Degrees to radians.
	parameter (r_d = 180./pi)		!Radians to degrees.
	parameter (alpha = 1./137.03599D0)	!Fine structure const.
	parameter (m_p = .93827231D0)		!Proton mass in GeV/c^2.
	parameter (m_e = .5109991D-03)		!Electron mass in GeV/c^2.
	parameter (m_amu = .9314943D0)		!Atomic mass unit in GeV/c^2.
	parameter (hc_2 = .389379292)	        !(hbar c)^2 GeV2-milibarns.

C ========================== END OF INCLUDE FILE ===============================
