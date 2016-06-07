CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	real*8 function sig_xiscale(ebeam,eprime,the)

	implicit none
C This function will (hopefully) return the diff. xsec PER NUCLEON
C in units of microbarns/GeV/sr

	real*8 ebeam,eprime,x,q2,the  !inputs
	real*8 mp,hbarc2,alpha  ! parameters

	real*8 nuW2,sig
	real*8 xi,nu,R,w1 ! calculated along the way

	external nuW2

	parameter(mp=0.938272)
	parameter(alpha=1./137.035)
	parameter(hbarc2 = 0.389379292E03)   ! GeV2-microbarns 

	nu = ebeam-eprime
	q2 = 4.0*ebeam*eprime*sin(the/2.0)**2
	x = q2/2./mp/nu
	xi = 2.*x/(1.+sqrt(1.+4.*mp**2*x**2/q2))
	R = 0.32/q2

	W1 = nuW2(xi)/nu*(1.+nu**2/q2)/(R+1)
	sig_xiscale = 4.*hbarc2*alpha**2*eprime**2/q2**2*
	1    (nuW2(xi)/nu*cos(the/2.)**2 + 2.*w1*sin(the/2.)**2)
	sig = sig_xiscale

c	write(6,*) 'eprime,xi,sig,nuW2,the',eprime,xi,sig,nuW2(xi)
	return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	real*8 function nuW2(xi)

	implicit none
C John's parameterization of xi scaling data - gives nuW2 for the deuteron.
C I'll Divide by 2 to get nuW2 per nucleon

C DJG modified by Dave G. to be a little better behaved at lower xi

	real*8 xi,a(4)

cdg	data a/1.8854,-9.4656,12.5166,-7.6104/
cdg	if(xi.lt.0.50) then
cdg	   nuW2=-0.25403+0.65533*xi-3.0678*xi*xi
cdg	   nuW2=10**nuW2
cdg      elseif (xi.lt.1.0) then
cdg          nuw2 = a(1)+a(2)*xi+a(3)*xi*xi+a(4)*xi*xi*xi
cdg          nuw2 = 10**nuw2
cdg        else
cdg          nuw2 = 4.4257 - 7.1394*xi
cdg          nuw2 = 10**nuw2
cdg        endif

	data a/-0.0026772,-1.0912,0.54484,-2.028/
        if (xi.lt.1.0) then
          nuw2 = a(1)+a(2)*xi+a(3)*xi*xi+a(4)*xi*xi*xi
          nuw2 = 10**nuw2
        else
          nuw2 = 4.4257 - 7.1394*xi
          nuw2 = 10**nuw2
        endif

	nuW2 = nuW2/2.0

	return
	end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



