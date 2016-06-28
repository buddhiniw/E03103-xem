	subroutine calc_xsec_from_delta
C+______________________________________________________________________________
!
! DESCRIPTION:
!
!       Subroutine to calculate the delta binned cross-sections.
!
!       05-18-2015 buddhini - Coppied from the original file (Aji's file) to clean up the code
!       a bit.
!       06-1602916 buddhini - Added errors on Fy
C-______________________________________________________________________________

	implicit none

	include 'hbook.inc'
	include 'kinematics.cmn'
	include 'cuts.inc'
	include 'histos.cmn'

	integer deltaid ! id for delta histograms
	integer i,id,idpos
	real*8  delta(ndeltabins), eprime(ndeltabins), xi(ndeltabins), xbj(ndeltabins)
        real*8  q2(ndeltabins),w(ndeltabins),y(ndeltabins),sig_mott,tau(ndeltabins)
        real*8  ratio(ndeltabins), eratio(ndeltabins)
	real*8  sigexp(ndeltabins),sigbornexp(ndeltabins),esig(ndeltabins),sigdis(ndeltabins)
	real*8  sigradmodel(ndeltabins),sigbornmodel(ndeltabins)
	real*8  sigpos(ndeltabins),esigpos(ndeltabins)
	real*8  rc(ndeltabins) !radiateive correction factor
	real*8  cc(ndeltabins) !coulomb correction factor
	real*8  sigrq(ndeltabins),sigre(ndeltabins),sigb(ndeltabins)
	real*8  binsize
	real*8  mp,mtar,me,m_amu
	real*8 xsec,y_sub
	real*8 rad_func,cc_func
	real*8 sig_rad_func, sig_dis_func,sig_vert_func,sig_vert_func2
	real*8 radcor
	real*8 coulcor2
	real*8 sigrqe,sigrel,sigvi
	real*4 hi,hie
	real*8 aux(7)

	logical y_calc_ok
	parameter (deltaid=3000)
	parameter (mp=0.938272)
	parameter (m_amu=0.9314943)
	parameter (me=0.0005109991)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Declare variables for the calculation on F(y) (Buddhini)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	real*8 sig_p(ndeltabins),sig_n(ndeltabins),dwdy,nu(ndeltabins),q3(ndeltabins)
	real*8 fy(ndeltabins),efy(ndeltabins),fact(ndeltabins)
	real*8 sig_p_sub,sig_n_sub,prefac

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Declare arguments for call to Hoehler(nform).
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	real*8 gep,gmp,gen,gmn,temp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C First, extract Data/MC ratios
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	call target_info(dble(targ_Z),dble(targ_a),aux)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C 1% delta bins
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	binsize = (hdeltacuthi-hdeltacutlo)/ndeltabins
	id = deltaid+13
	idpos=deltaid+15

c	write(33,*) 'theta in calc_xsec',th_deg,thdegcent

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Mott cross section.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	   sig_mott = cos(thdegcent/2.)/(274.072*ebeam_cor*sin(thdegcent/2.)**2.)
	   sig_mott = (.019733*sig_mott)**2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Loop through the histograms in the simulation file and get the kinematics. 
c hi(id,j) returns the channel contents in a given histogram bin
c
c hie(id,j) Returns the value of the error that has been stored in a given 1-dimensional histogram channel.
c This corresponds to the square root of the sum of the squares of the weights.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	do i=1,ndeltabins
	   ratio(i) = dble(hi(id,i)) 
	   eratio(i) = dble(hie(id,i))
	   delta(i) = hdeltacutlo + (i-1.0)*binsize + binsize/2.0
	   eprime(i) = pcent*(1.0+delta(i)/100.0)
	   q2(i) = 4.*ebeam_cor*eprime(i)*sin(thradcent/2)**2
	   tau(i) = q2(i)/(4*mp**2)

	   nu(i) = ebeam_cor-eprime(i)
	   q3(i)=sqrt(q2(i) + nu(i)**2)

	   xbj(i) = q2(i)/2./mp/(ebeam_cor-eprime(i))
	   xi(i) = 2.*xbj(i)/(1.+sqrt(1.+4.*mp**2*xbj(i)**2/q2(i)))
	   w(i) = -q2(i) + mp**2 + 2.*mp*(ebeam_cor-eprime(i))
	   mtar = targ_matom*m_amu - targ_z*me

c	   write(6,*) 'ebeam_cor,eprime',ebeam_cor,eprime


c calculate the y scaling variable
c y is in units of GeV/c
	   call y_calc(ebeam_cor,eprime(i),thradcent,mtar,mp,aux(2),y_sub,y_calc_ok)
	   y(i) = y_sub
c	   write(*,*) 'y_sub =  ',y_sub
	   

c Calculate the radiated cross section from the "main" 2D radcor file
C keep things consistent - "1" goes with externals, "2" with XEM

	   sigradmodel(i) = sig_rad_func(eprime(i),thradcent)!XEM model, 2d binning
	   sigexp(i) = ratio(i)*sigradmodel(i) !radiated experimental cross section
	   esig(i) = eratio(i)*sigradmodel(i) ! error on the radiated exp cross section

c	   write(6,*)'radiated experimental cross section (sigexp rad)',sigexp(i)	   

c Calculate the born experimental cross section. 
c For this assume the ratio of the radiated cross sections
c is equal to the ratio of the born cross sections.
C keep things consistent - "1" goes with externals, "2" with XEM
	   sigbornmodel(i) = sig_vert_func(eprime(i),thradcent) !XEM model, 2d binning
c	   write(6,*)'Born model cross section (sigborn model)',sigbornmodel(i)	   

	   sigbornexp(i)=ratio(i)*sigbornmodel(i)
c	   write(6,*)'Born experimental cross section (sigborn rad)',sigbornexp(i)	   

c Get the dis cross-section from XEM model
	   sigdis(i) = sig_dis_func(eprime(i),thradcent)!XEM model, 2d binning
c	   write(*,*) 'dis cross-section from XEM model (sigdis) =  ',sigdis(i)



c Get the radiative correction
	   rc(i) = rad_func(eprime(i),thradcent) !XEM, no CCRAD
	   
c Get the coulomb correction
	   cc(i) = cc_func(eprime(i),thradcent)

c Calculate the error on the model 


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Buddhini 01-08-2016
c Calculate F(Y) and Y
c These will be used to extract the scaling function F(y)
c from the data using the formula 2.27 in N. Formins thesis
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	   

c  Calculate the DeForest based off-shell cross sections averaged over the
c  angle PHI, which is the angle between the electron scattering plane and
c  the plane containing the initial and final momentum of the struck nucleon.
	   
c Get sig_p and sig_n in units of mB/ster
c The energy/momentum inputs must be in GeV
c sigexp and sigdis are in units of nB/(ster MeV)

	if(targ_a.gt.1.0) then
	   call sig_bar_df_dble(ebeam_cor, eprime(i), thdegcent, y(i), 0.0, sig_p_sub, sig_n_sub)
	   sig_p_sub=sig_p_sub
	   sig_n_sub=sig_n_sub
	else
	   call nform (12.0,q2(i),gep,gen,gmp,gmn)
	   prefac = sig_mott*eprime(i)/ebeam_cor
	   sig_p_sub = (gep**2+tau(i)*gmp**2)/(1.0+tau(i)) + 2.0*tau(i)*gmp**2*tan(thdegcent/2.)**2
	   sig_p_sub = prefac*sig_p_sub
	   sig_n_sub = 0.
	endif

	   sig_p(i) = sig_p_sub
	   sig_n(i) = sig_n_sub

c Get dwdy
c	   write(6,*) "sig_p and sig_n and targ_z, targ_a ",sig_p(i), sig_n(i),targ_Z,targ_a
c	   write(6,*) "q3, mp, y ",q3(i), mp,y(i)

	   dwdy = q3(i)/sqrt(mp**2+q3(i)**2+y(i)**2+2.0*q3(i)*y(i))

	   if(targ_a.gt.1.0) then
	      fact(i) = (targ_Z*sig_p(i)+(targ_a-targ_Z)*sig_n(i))/dwdy
	   else
	      fact(i) = targ_Z*sig_p(i)/dwdy
	   endif
c          write(*,*) 'fact =  ',fact(i)
	   
	   
c Subtract the dis cross-section from the experimental Born cross-section and divide by fact.
c from radcor_xem/sigmodel_calc.f - the sigbornexp(i) is in units of nB/(ster-MeV) or muB/(ster-GeV)
c So fy gets calculated in units of 1/(TeV/c). 
	   fy(i) = (sigbornexp(i)-sigdis(i))/fact(i)
c	   write(*,*) 'calculated fy =  ',fy(i)

c Calculate error on fy using esig
	   efy(i) = esig(i)/fact(i)
c	   write(*,*) 'calculated efy =  ',efy(i)

	enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C do it all again for positrons
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	do i=1,ndeltabins
	   ratio(i) = dble(hi(idpos,i))
	   eratio(i) = dble(hie(idpos,i))
	   sigpos(i) = sigexp(i)*ratio(i)
	   esigpos(i) = sigexp(i)*eratio(i)
	enddo
	write(44,*)'Momentums and energies and y are in GeV, cross sections are in nB/(ster MeV)'
	write(44,*)'fy is in 1/TeV, factor is in ster/mB '


	write(44,7) 'Eprime ',' xBj   ','xi     ','W2     ','Q2     ','y     ',
	1    'fy      ','efy         ',  'factor     ','sig_rad   ','sig_dis    ',
	2    'esig      ', 'sigmodel  ','sigmod_rad ',

c	write(44,7) '(GeV) ','    ','    ','(GeV^2)   ','(GeV^2)    ','(GeV/c)     ',
c	1    '(TeV^-1)     ','(TeV^-1)     ','(ster/mB)     ',' [nB/(ster MeV)]  ',' [nB/(ster MeV)]   ',
c	2    '             ','[nB/(ster MeV)] ','[nB/(ster MeV)] ',

	do i=1,nxbins
	   if(sigexp(i).gt.0.0) then
	      write(44,8) eprime(i),xbj(i),xi(i),w(i),q2(i),y(i),fy(i),efy(i), fact(i),sigexp(i),sigdis(i),
	1	   esig(i),sigbornmodel(i),sigradmodel(i),
	   endif
	enddo


	close (44)



 7	format(6(a7,2x),13(a10,2x))
 8	format(6(f7.4,2x),13(e10.4e2,2x))

	return 

	end

