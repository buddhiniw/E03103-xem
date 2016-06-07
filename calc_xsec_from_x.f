	subroutine calc_xsec_from_x
C+______________________________________________________________________________
!
! DESCRIPTION:
!
!   Subroutine to calculate x binned cross-sections
!
!       05-12-2015 buddhini - Coppied from the original file (Aji's file) to clean up the code
!       a bit.
C-______________________________________________________________________________


	implicit none

	include 'hbook.inc'
	include 'kinematics.cmn'
	include 'cuts.inc'
	include 'histos.cmn'

	integer xid ! id for delta histograms
	integer i,id,idpos
	real*8  eprime(nxbins), xi(nxbins), xbj(nxbins),q2(nxbins),w(nxbins),y(nxbins)
        real*8  ratio(nxbins), eratio(nxbins)
	real*8  sigmodel(nxbins),sigradmodel(nxbins),sigexp(nxbins),esig(nxbins)
	real*8  cc1(nxbins),cc2(nxbins),rc1(nxbins),rc2(nxbins)
	real*8  sigpos(nxbins),esigpos(nxbins)
	real*8  sigrq(nxbins),sigre(nxbins),sigb(nxbins)
	real*8  xsec,y_sub
	real*8  binsize,counts
	real*8  mp,me,m_amu,mtar
	real*8  radcor1,radcor2,coulcor1,coulcor2
	real*8  sigrqe,sigrel,sigvi
	real*8  sig_vert_func2,sig_rad_func
	real*8  rad_func_2,rad_func,cc_func,cc_func_2
	real*4  hi,hie
	real*8 aux(7)

	logical y_calc_ok
	parameter (xid=4000)
	parameter (mp=0.938272)
	parameter (m_amu=0.9314943)
	parameter (me=0.0005109991)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C First, extract Data/MC ratios
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	call target_info(dble(targ_Z),dble(targ_a),aux)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C X bin size
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	binsize = (xbjmax-xbjmin)/nxbins
	id = xid+13
	idpos=xid+15


	do i=1,nxbins
	   sigexp(i) =0.0
	   esig(i) = 0
	   counts = dble(hi(id-4,i))
	   ratio(i) = dble(hi(id,i))
	   eratio(i) = dble(hie(id,i))
	   xbj(i) = xbjmin + (i-1.0)*binsize + binsize/2.0
	   eprime(i) = (2.*ebeam_cor*xbj(i)*mp)/
	1	(4.*ebeam_cor*sin(thradcent/2)**2 + 2.*xbj(i)*mp)
	   q2(i) = 4.*ebeam_cor*eprime(i)*sin(thradcent/2)**2
	   xi(i) = 2.*xbj(i)/(1.0+sqrt(1.+4.0*mp**2*xbj(i)**2/q2(i)))
	   w(i) = -q2(i) + mp**2 + 2.*mp*(ebeam_cor-eprime(i))
	   mtar = targ_matom*m_amu - targ_z*me
	   call y_calc(ebeam_cor,eprime(i),thradcent,mtar,mp,aux(2),y_sub,y_calc_ok)
	   y(i) = y_sub

	   xsec = sig_rad_func(eprime(i),thradcent)


	   radcor1 = rad_func_2(eprime(i),thradcent,sigrqe,sigrel,sigvi) !PB "externals_all"
	   radcor2 = rad_func(eprime(i),thradcent) !XEM, no CCRAD

C keep things consistent - "1" goes with externals, "2" with XEM

	   coulcor1 = cc_func_2(eprime(i),thradcent)
	   coulcor2 = cc_func(eprime(i),thradcent)

	   sigmodel(i) = sig_vert_func2(eprime(i),thradcent)
	   sigradmodel(i) = sig_rad_func(eprime(i),thradcent)
	   sigexp(i) = ratio(i)*xsec
	   esig(i) = eratio(i)*xsec
	   cc1(i) = coulcor1
	   cc2(i) =coulcor2
	   rc1(i) = radcor1
	   rc2(i) = radcor2

	   sigrq(i) = sigrqe
	   sigre(i) = sigrel
	   sigb(i) = sigvi
	enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Do it again for positrons
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	do i=1,nxbins
	   ratio(i) = dble(hi(idpos,i))
	   eratio(i) = dble(hie(idpos,i))
	   sigpos(i) = sigexp(i)*ratio(i)
	   esigpos(i) = sigexp(i)*eratio(i)
	enddo


	write(46,7) 'Eprime ',' xBj     ','xi     ','W2     ','Q2     ','y     ','sig_rad   ','esig      ',
	1     'sigmodel  ','sigmod_rad ','sig_pos   ','esigpos   ','cc_extpb  ','cc_xem    ','rc-full   ','rc-epeak  ',
	2     'sigr-qe   ','sigr-el    ','sigb_in   '   
	do i=1,nxbins
	   if(sigexp(i).gt.0.0) then
	      write(46,8) eprime(i),xbj(i),xi(i),w(i),q2(i),y(i),sigexp(i),esig(i),sigmodel(i),sigradmodel(i),
	1	   sigpos(i),esigpos(i),cc1(i),cc2(i),rc1(i),rc2(i),sigrq(i),sigre(i),sigb(i)
	   endif
	enddo
	close (46)


 7	format(6(a7,2x),13(a10,2x))
 8	format(6(f7.4,2x),13(e10.4e2,2x))

	return 

	end

