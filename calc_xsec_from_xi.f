	subroutine calc_xsec_from_xi

	implicit none

	include 'hbook.inc'
	include 'kinematics.cmn'
	include 'cuts.inc'
	include 'histos.cmn'

	integer xiid ! id for delta histograms
	integer i,id,idpos
	real*8  eprime(nxibins), xi(nxibins), xbj(nxibins),q2(nxibins),w(nxibins),y(nxbins)
        real*8  ratio(nxibins), eratio(nxibins)
	real*8  sigmodel(nxibins),sigradmodel(nxibins),sigexp(nxibins),esig(nxibins)
	real*8  cc1(nxibins),cc2(nxibins),rc1(nxibins),rc2(nxibins)
	real*8  sigpos(nxibins),esigpos(nxibins)
	real*8  xsec,y_sub
	real*8  binsize,counts
	real*8  mp,me,m_amu,mtar
	real*8  radcor1,radcor2,radcor3,coulcor1,coulcor2
	real*8  sig_xiscale,sig_vert_func,sig_vert_func2,sig_rad_func
	real*8  rad_func_2,rad_func,rad_func_4,cc_func,cc_func_2
	real*4  hi,hie
	real*8 aux(7)

	logical y_calc_ok
	parameter (xiid=4100)
	parameter (mp=0.938272)
	parameter (m_amu=0.9314943)
	parameter (me=0.0005109991)


C First, extract Data/MC ratios
	call target_info(dble(targ_Z),dble(targ_a),aux)

C       First, extract xi-binned results
C Xi bin size
	binsize = (ximax-ximin)/nxibins
	id = xiid+13
	idpos=xiid+15

	write(33,*) 'theta in calc_xsec_from_xi',th_deg,thdegcent

c	write(6,*) 'xmin etc.',ximin,ximax,nxibins
	do i=1,nxibins
	   sigexp(i) =0.0
	   esig(i) = 0
	   counts = dble(hi(id-4,i))
	   ratio(i) = dble(hi(id,i))
	   eratio(i) = dble(hie(id,i))
	   xi(i) = ximin + (i-1.0)*binsize + binsize/2.0
	   eprime(i) = (2.*ebeam_cor*xi(i)*mp+xi(i)**2*mp**2)/
	1	(4.*ebeam_cor*sin(thradcent/2)**2 + 2.*xi(i)*mp)
	   q2(i) = 4.*ebeam_cor*eprime(i)*sin(th_rad/2)**2
	   xbj(i) = q2(i)/2./mp/(ebeam_cor-eprime(i))
	   w(i) = -q2(i) + mp**2 + 2.*mp*(ebeam_cor-eprime(i))
	   mtar = targ_matom*m_amu - targ_z*me
	   call y_calc(ebeam_cor,eprime(i),thradcent,mtar,mp,aux(2),y_sub,y_calc_ok)
	   y(i) = y_sub

c	   write(6,*) 'here i am',i,xi(i),eprime(i)

	   xsec = sig_rad_func(eprime(i),thradcent)


	   radcor1 = rad_func_2(eprime(i),thradcent) !PB "externals_all"
	   radcor2 = rad_func(eprime(i),thradcent) !XEM, no CCRAD
c	   radcor3 = rad_func_4(eprime(i),thradcent) !XEM, with CCRAD
	   coulcor2 = cc_func_2(eprime(i),thradcent)
	   coulcor1 = cc_func(eprime(i),thradcent)
c	   if(radcor2.gt.0.0) then
c	      coulcor2 = radcor3/radcor2
c	   else
c	      coulcor2 = 1.0
c	   endif

	   sigmodel(i) = sig_vert_func2(eprime(i),thradcent)
	   sigradmodel(i) = xsec
	   sigexp(i) = ratio(i)*xsec
	   esig(i) = eratio(i)*xsec
	   cc1(i) = coulcor1
	   cc2(i) =coulcor2
	   rc1(i) = radcor1
	   rc2(i) = radcor2
c	   write(6,*) 'here i am', xsec,ratio(i),sigexp(i)
	enddo

C Do it again for positrons

	do i=1,nxibins
	   ratio(i) = dble(hi(idpos,i))
	   eratio(i) = dble(hie(idpos,i))
	   sigpos(i) = sigexp(i)*ratio(i)
	   esigpos(i) = sigexp(i)*eratio(i)
	enddo


	write(48,7) 'Eprime ',' xBj     ','xi     ','W2     ','Q2     ','y     ','sig_rad   ','esig      ',
	1     'sigmodel  ','sigmod_rad ','sig_pos   ','esigpos   ','cc_norad  ','cc_rad    ','rc-full   ','rc-epeak  '   
	do i=1,nxibins
	   if(sigexp(i) .gt.0.0) then
	      write(48,8) eprime(i),xbj(i),xi(i),w(i),q2(i),y(i),sigexp(i),esig(i),sigmodel(i),sigradmodel(i),
	1	   sigpos(i),esigpos(i),cc1(i),cc2(i),rc1(i),rc2(i)
	   endif
	enddo

	close(48)

 7	format(6(a7,2x),10(a10,2x))
 8	format(6(f7.4,2x),10(e10.4e2,2x))

	return 

	end

