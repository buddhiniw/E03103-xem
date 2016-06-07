	subroutine calc_eloss(targ_A,targ_Z,eloss,eloss_corr,ep,th_ev,ytar)
C+______________________________________________________________________________
!
! DESCRIPTION:
!
!   Subroutine to calculate energy loss in the target
!
!
!       05-12-2015 buddhini - Coppied from the original file (Aji's file) to clean up the code
!       a bit.
C-______________________________________________________________________________

	implicit none

	include 'target.inc'
	include 'math_physics.inc'

	real*8 targ_A,targ_Z,eloss,eloss_corr,ep,th_ev,ytar
	real*8 zoff,length ! hardwired stuff from the Monte Carlo

	real*8 ztar,t,atmp,btmp,ctmp
	real*8 z_can,side_path,costmp
	real*8 th_can

	real*8 energy,m_elec

	real*8 s_Al,s_Al_targ,s_air,s_kevlar,s_mylar,s_target,targ_rho
	real*8 eloss_target,eloss_al,eloss_air,eloss_mylar,eloss_kevlar

	integer typeflag 
	

	m_elec = 1000.0*m_e
	if(targ_A.le.4) then
	   zoff = 0.18
	   length = 3.926
	elseif(targ_A.eq.9) then
	   zoff = 0.25
	   length = 1.012
	elseif(targ_a.eq.12) then
	   zoff = 0.25
	   length = 0.294
	elseif(targ_A.eq.64) then
	   zoff = 0.13
	   length = 0.089
	else
	   zoff = 0.13
	   length = 0.0196
	endif

	if(targ_A.eq.1.0) then
	   targ_rho = 0.0723
	elseif(targ_A.eq.2.0) then
	   targ_rho=0.167
	elseif(targ_A.eq.3.0) then
	   targ_rho=0.0708
	elseif(targ_A.eq.4) then
	   targ_rho=0.135
	elseif(targ_A.eq.9) then
	   targ_rho=1.848
	elseif(targ_A.eq.12) then
	   targ_rho=2.265
	elseif(targ_A.eq.64) then
	   targ_rho=8.96
	elseif(targ_A.eq.27) then
	   targ_rho=2.70
	endif

	ztar= ytar/sin(th_ev)

	if (targ_A.le.4.) then   

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C DJG Here I'm just copying stuff from SIMC - hope it works
C JRA this is ugly.  Solve for z position where particle intersects can.  The
C JRA pathlength is then (z_intersect - z_scatter)/cos(theta)
C JRA Angle from center to z_intersect is acos(z_intersect/R).  Therefore the
C JRA angle between the particle and wall is pi/2 - (theta - theta_intersect)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	     t=tan(th_ev)**2
	     atmp=1+t
	     btmp=-2*(ztar-zoff)*t
	     ctmp=(ztar-zoff)**2*t-(length/2.)**2
	     z_can=(-btmp+sqrt(btmp**2-4.*atmp*ctmp))/2./atmp
	     side_path = (z_can - (ztar-zoff))/abs(cos(th_ev))
	     s_target = side_path
	     costmp=z_can/(length/2.)
	     if (abs(costmp).le.1) then
		th_can=acos(z_can/(length/2.))
	     else if (abs(costmp-1.).le.0.000001) then
		th_can=0.	!extreme_trip_thru_target can give z/R SLIGHTLY>1.0
	     else
		stop 'z_can > can radius in target.f !!!'
	     endif
	     s_Al_targ =  0.0050*2.54/abs(sin(pi/2 - (th_ev - th_can)))
	  else
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Case 2 solid target
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	     s_al_targ = 0.0
	     s_target = abs(length/2. - (ztar-zoff))/cos(th_ev)
	  endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Scattering before magnets:  Approximate all scattering as occuring AT TARGET.
C  16 mil Al scattering chamber window (X0=8.89cm)
C  15(?) cm air (X0=30420cm)
C spectrometer entrance window
C  15 mil Kevlar (X0=74.6cm)
C  5 mil Mylar  (X0=28.7cm)             Total of 0.60% rad. length.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	  s_Al = s_Al_targ + 0.016*2.54
	  s_air = 15.0
	  s_kevlar = 0.015*2.54
	  s_mylar = 0.005*2.54


	  energy=1000.0*ep
	  typeflag=1
	  eloss=0
	  call enerloss_new(s_target,targ_rho,targ_Z,targ_A,energy,m_elec,
	1    typeflag,Eloss_target)
	  call enerloss_new(s_Al,rho_Al,Z_Al,A_Al,energy,m_elec,typeflag,Eloss_Al)
	  call enerloss_new(s_air,rho_air,Z_air,A_air,energy,m_elec,typeflag,
	1    Eloss_air)
	  call enerloss_new(s_kevlar,rho_kevlar,Z_kevlar,A_kevlar,energy,
	1    m_elec,typeflag,Eloss_kevlar)
	  call enerloss_new(s_mylar,rho_mylar,Z_mylar,A_mylar,energy,m_elec,
	1    typeflag,Eloss_mylar)
	  Eloss = Eloss_target + Eloss_Al + Eloss_air + Eloss_kevlar + Eloss_mylar

	  typeflag=4
	  eloss_corr=0
	  s_target = length/2.0
	  s_Al =  s_Al - s_Al_targ
	  if(targ_A.le.4) then
	     s_Al = s_Al + 0.0050*2.54
	  endif
	  call enerloss_new(s_target,targ_rho,targ_Z,targ_A,energy,m_elec,
	1    typeflag,Eloss_target)
	  call enerloss_new(s_Al,rho_Al,Z_Al,A_Al,energy,m_elec,typeflag,Eloss_Al)
	  call enerloss_new(s_air,rho_air,Z_air,A_air,energy,m_elec,typeflag,
	1    Eloss_air)
	  call enerloss_new(s_kevlar,rho_kevlar,Z_kevlar,A_kevlar,energy,
	1    m_elec,typeflag,Eloss_kevlar)
	  call enerloss_new(s_mylar,rho_mylar,Z_mylar,A_mylar,energy,m_elec,
	1    typeflag,Eloss_mylar)
	  Eloss_corr = Eloss_target + Eloss_Al + Eloss_air + Eloss_kevlar + Eloss_mylar

	  eloss = eloss/1000.0
	  eloss_corr = eloss_corr/1000.0

	  return
	  end
