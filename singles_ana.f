	program singles_ana
	
	implicit none

	include 'hbook.inc'
	include 'cuts.inc'
	include 'radcor.cmn'
	include 'extcor.cmn'
	include 'histos.cmn'
	include 'cryocor.cmn'
	include 'math_physics.inc'
	include 'kinematics.cmn'
	
	character*80 datafile,dummyfile,posifile,dumposifile
	character*80 xfile,simfile,hbookfile,radfile,extfile
	character*80 line,title
	integer i,j,datarun(200),dummyrun(200),dummyflag(200)
	integer posirun(200),dumposirun(200)
	integer numdatruns,numdumruns,numposiruns,numdumposiruns,run,filelen
	integer dataid,dumid,newid,simid
	integer ncomment
	real*8 dataweight, dummyweight, posiweight
	real*8 counts,charge,charge_elcl
	real*8 totcounts,totcharge,cleancharge
	real*8 dummycharge,dummycounts,cleandummycharge
	real*8 posicharge,posicounts
	real*8 dummyposicharge,dummyposicounts
	real*8 totalcounts
	real*8 cdummypos,posicdummypos
	real*8 cdens,ccurrent,cbeampos,cryocor,posicryocor
	real*8 hsum
	real*8 edum,q2dum,xdum,sig1dum,sig2dum,sig3dum,dum1,dum2,dum3
	real*8 q2min,q2max,epmin,epmax

c	real*8 dumchi2(6),dumpar(1)
	
	real*8 meanchi2,cutratio

	real*8 x1,x1err,x2,x2err
	real*8 xmin,xmax,ymin,ymax
	real*8 chisquare(25),nchannel(25),meanchisq
	real*4 hi,hie,one
	integer nx,ny,nwt,loc
	external hi,hie

	integer*4 istat,lrecl
	logical flag,pflag,doing_solid,doing_cryo,doing_posi,doing_dummy_posi
	logical makedatahistos


	parameter(one=1.0)

	ytar = 0.0

	cryocor =1.0
	posicryocor =1.0

C       Open input file - read in file names and parameters
	open (unit=12,file='singles_ana.inp',status='old')
	read(12,'(a)') datafile     !data runs
	read(12,'(a)') dummyfile    !dummy runs
	read(12,'(a)') posifile     !positron runs
	read(12,'(a)') dumposifile  !dummy positron runs
	read(12,'(a)') simfile      !Monte Carlo file
	read(12,'(a)') radfile      !radiative corrections table file
	read(12,'(a)') extfile      !external radcor: dummy/walls ratio
	close (12)

	doing_solid=.FALSE.
	doing_cryo=.FALSE.
	doing_posi=.FALSE.
	doing_dummy_posi=.FALSE.



C       Open run list - read in run numbers
	do i=32,1,-1
	   if (datafile(i:i).eq.' ') filelen=i-1
	enddo
C       Open efficiency file..
	xfile = 'OUT/'//datafile(1:filelen-4)//'.eff'
	open (unit=11,file=xfile,status='unknown')
	write(11,1) 'RUN  ', 'HFID ','HFIDE','CO_LT','HELEC','HTRIG','HCERE',
	1  'HCALE','TrhoC','BO_Q'

C       Open output file
	xfile = 'OUT/'//datafile(1:filelen-4)//'.out'
	open (unit=33,file=xfile,status='unknown')

C       Open files for delta cross sections
	xfile = 'OUT/'//datafile(1:filelen-4)//'.deltabin.xsec'
	open (unit=44,file=xfile,status='unknown')

C       Open files for eprime cross sections
c	xfile = 'OUT/'//datafile(1:filelen-4)//'.eprimebin.xsec'
c	open (unit=45,file=xfile,status='unknown')

C       Open files for x-binned cross sections
	xfile = 'OUT/'//datafile(1:filelen-4)//'.xbin.xsec'
	open (unit=46,file=xfile,status='unknown')

C       Open files for xi cross sections
c	xfile = 'OUT/'//datafile(1:filelen-4)//'.xibin.xsec'
c	open (unit=48,file=xfile,status='unknown')

C       Open files for 2d binned cross sections
c	xfile = 'OUT/'//datafile(1:filelen-4)//'.2dbin.xsec'
c	open (unit=49,file=xfile,status='unknown')


C test if we are doing dummy stuff as target
	xfile = datafile(12:19)

	if(INDEX(xfile,'dummy_us').gt.0) then
	   ytar = -1.0
	elseif (INDEX(xfile,'dummy_ds').gt.0) then
	   ytar = 1.0
	endif


C       Define hbook file
	hbookfile = 'hbook/'//datafile(1:filelen-4)//'.hbook'

	xfile = 'RUNS/'//datafile(1:filelen)
	write(33,*) "opening ",xfile
	write(6,*) "opening ",xfile
	open (unit=10,file=xfile,status='old')
	do i = 1,200
	   read(10,*,end=100) datarun(i)
	   dummyflag(i) = 0
	enddo
 100	continue
	numdatruns = i-1
	close(10)

	if(numdatruns.lt.1) then
	   write(33,*) 'No runs for this target at this setting - abort'
	   write(6,*) 'No runs for this target at this setting-abort'
	   stop
	endif

C     Open dummy run list.
	do i=32,1,-1
	   if (dummyfile(i:i).eq.' ') filelen=i-1
	enddo
	if(filelen.gt.0) then
	   doing_cryo=.true.
	   xfile = 'RUNS/'//dummyfile(1:filelen)
	   write(33,*) "opening ",xfile
	   write(6,*) "opening ",xfile
	   open (unit=10,file=xfile,status='old')
	   do i = 1,200
	      read(10,*,end=200) dummyrun(i)
	      dummyflag(i) = 1
	   enddo
 200	   continue
	   numdumruns = i-1
	   close(10)
	   if(numdumruns.lt.1) then
	      write(33,*) 'No dummy runs even thought you gave me a runlist!'
	      write(6,*) 'No dummy runs even thought you gave me a runlist!'
	      stop
	   endif
	endif

C     Open positron run list.
	do i=32,1,-1
	   if (posifile(i:i).eq.' ') filelen=i-1
	enddo
	if(filelen.gt.0) then
	   doing_posi=.true.
	   xfile = 'RUNS/'//posifile(1:filelen)
	   write(33,*) "opening ",xfile
	   write(6,*) "opening ",xfile
	   open (unit=10,file=xfile,status='old')
	   do i = 1,200
	      read(10,*,end=300) posirun(i)
c	      dummyflag(i) = 1
	   enddo
 300	   continue
	   numposiruns = i-1
	   close(10)
	   if(numposiruns.lt.1) then
	      write(33,*) 'No positron runs even thought you gave me a runlist!'
	      write(6,*) 'No positron runs even thought you gave me a runlist!'
c	      stop
	   endif
	endif
	if(numposiruns.lt.1) then
	   doing_posi=.FALSE.
	endif
C set doing_posi to false for Nadia test
c	doing_posi=.FALSE.

C     Open dummy positron run list.
	do i=32,1,-1
	   if (dumposifile(i:i).eq.' ') filelen=i-1
	enddo
	if(filelen.gt.0) then
	   doing_dummy_posi = .true.
	   xfile = 'RUNS/'//dumposifile(1:filelen)
	   write(33,*) "opening ",xfile
	   write(6,*) "opening ",xfile
	   open (unit=10,file=xfile,status='old')
	   do i = 1,200
	      read(10,*,end=400) dumposirun(i)
c	      dummyflag(i) = 1
	   enddo
 400	   continue
	   numdumposiruns = i-1
	   close(10)
	   if(numdumposiruns.lt.1) then
	      write(33,*) 'No dummy positron runs even thought you gave me a runlist!'
	      write(6,*) 'No dummy psoitron runs even thought you gave me a runlist!'
c	      stop
	   endif
	endif
	if(numdumposiruns.lt.1) then
	   doing_dummy_posi=.FALSE.
	endif


CCCCC RC FILE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC USE this for bin centering (2D grid)
C       Open radiative corrections table
	do i=32,1,-1
	   if (radfile(i:i).eq.' ') filelen=i-1
	enddo

	if(filelen.gt.0) then
	   doing_rad=.true.
	   write(33,*) 'Radiative corrections file in input file'
	   xfile = 'RAD/'//radfile(1:filelen)
	   write(33,*) "opening ",xfile
	   open (unit=10,file=xfile,status='old')
	   if(abs(ytar).gt.0.01) then
	      ncomment=30
	   else
	      ncomment=28
	   endif
	   do i=1,ncomment
	      read(10,'(a)',end=99) line
	      if(line(1:3).ne.' **') then
		 write(6,*) 'Warning - expected comment line!!'
	      endif
	   enddo  !loop over comment lines

	   read(10,*) kmax  !number of p bins
	   write(33,*) 'number of momentum bins - RAD1',kmax

c	   ntheta = 91
c	   ntheta = 81
	   ntheta = 69
	   do i = 1,ntheta	!theta bins
	      do j=1,kmax	!loop over eprime dbins
		 read(10,1005,end=99) edum,eprad(j),thrad(i),q2dum,xdum,sig1dum,
	1	      sig2dum,sigv(i,j),sigrad(i,j),radcor(i,j),coulcor(i,j),dum1,
	2	      dum2,dum3
		 if(radcor(i,j).eq.10.0) then
		    radcor(i,j) = 1.0E30
		 endif
                      write(9,*) 'radcor read in:',i,j,eprad(j),thrad(i),radcor(i,j),sigrad(i,j),sigv(i,j),coulcor(i,j)
	      enddo		! loop over eprime
	   enddo !loop over theta
 99	   continue
	   close(10)
	endif
c	write(6,*) 'done w/first rc file'
CCCC END RC FILECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCC 2ND RC FILE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCC This reads in Peter's "externals_all" output
cc	write(6,*) 'opning 2nd rc file'
c	if(filelen.gt.0) then
c	   doing_rad=.true.
c	   write(33,*) '2nd Radiative corrections file in input file'
c	   write(6,*) '2nd Radiative corrections file in input file'
c	   xfile = 'RAD2/'//radfile(1:filelen)
c	   write(33,*) "opening ",xfile
c	   open (unit=10,file=xfile,status='old')
c	   do i=1,26
c	      read(10,'(a)',end=499) line
c	      if(line(1:3).ne.' **') then
c		 write(6,*) 'Warning - expected comment line!!'
c	      endif
c	   enddo  !loop over comment lines
c	   
c	   read(10,*) kmax2  !number of p bins
c	   ntheta = 69
c	   do j=1,kmax2		!loop over eprime dbins
cCDJG Nadia/Aji edited
c	      read(10,1004,end=499) edum,eprad2(j),thrad2,xdum,q2dum,
c	1	   sigv2(j),sigvi2(j),sigvqe2(j),sig3dum,sigrad2(j),sigrel2(j),sigrqe2(j),sigrdis2(j),coulcor2(j)
c	      if(sigrad2(j).gt.0.0) then
c		 radcor2(j) = sigv2(j)/sigrad2(j)
c	      else
c		 radcor2(j) = 0.0
c	      endif
cc             write(6,*) 'radcor read in 2:',i,j,eprad2(j),thrad2,radcor2(j),sigrad2(j),sigv2(j)
c	   enddo		! loop over eprime
c 499	   continue
c	   close(10)
c	endif
CCCCC END 2ND RC FILECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CC 3RD RC FILE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC XEM model - improved with no CC on radiated cross section
cdg	if(filelen.gt.0) then
cdgc	   doing_rad=.true.
cdg	   write(33,*) '3rd Radiative corrections file in input file'
cdg	   xfile = 'RAD3/'//radfile(1:filelen)
cdg	   write(33,*) "opening ",xfile
cdg	   open (unit=10,file=xfile,status='old')
cdg	   do i=1,26
cdg	      read(10,'(a)',end=599) line
cdg	      if(line(1:3).ne.' **') then
cdgc	      if(line(1:26).ne.' **') then
cdg		 write(6,*) 'Warning - expected comment line!!'
cdg	      endif
cdg	   enddo  !loop over comment lines
cdg	   
cdg	   read(10,*) kmax3  !number of p bins
cdgc	   write(6,*) 'number of momentum bins - RAD2',kmax
cdg
cdgc	   ntheta = 91
cdgc	   ntheta = 81
cdg	   ntheta = 69
cdgc	   do i = 1,ntheta	!theta bins
cdg	      do j=1,kmax3	!loop over eprime dbins
cdgCDJG Nadia/Aji edited
cdgc		 read(10,1004,end=499) edum,eprad2(j),thrad2,xdum,q2dum,
cdgc	1	      sigv2(j),sigrad2(j),sigrel2(j),sigrqe2(j),sigrdis2(j)
cdg		 read(10,1005,end=499) edum,eprad3(j),thrad3,q2dum,xdum,sig1dum,
cdg	1	      sig2dum,sigv3(j),sigrad3(j),radcor3(j),coulcor3(j),dum1,
cdg	2	      dum2,dum3
cdg		 if(sigrad3(j).gt.0.0) then
cdg		    radcor3(j) = sigv3(j)/sigrad3(j)
cdg		 else
cdg		    radcor3(j) = 0.0
cdg		 endif
cdgc		 write(6,*) 'radcor read in 3:',i,j,eprad3(j),thrad3,radcor3(j),sigrad3(j),sigv3(j)
cdg	      enddo		! loop over eprime
cdgc	   enddo !loop over theta
cdg 599	   continue
cdg	   close(10)
cdg	endif
cdgCCCCC END 3RD RC FILECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC 4TH RC FILE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC XEM model - improved with no CC on radiated cross section
cdg	if(filelen.gt.0) then
cdg	   doing_rad=.true.
cdg	   write(33,*) '4th Radiative corrections file in input file'
cdg	   xfile = 'RAD4/'//radfile(1:filelen)
cdg	   write(33,*) "opening ",xfile
cdg	   open (unit=10,file=xfile,status='old')
cdg	   do i=1,26
cdg	      read(10,'(a)',end=699) line
cdg	      if(line(1:3).ne.' **') then
cdgc	      if(line(1:26).ne.' **') then
cdg		 write(6,*) 'Warning - expected comment line!!'
cdg	      endif
cdg	   enddo  !loop over comment lines
cdg	   
cdg	   read(10,*) kmax4  !number of p bins
cdgc	   write(6,*) 'number of momentum bins - RAD2',kmax
cdg
cdgc	   ntheta = 91
cdgc	   ntheta = 81
cdg	   ntheta = 69
cdgc	   do i = 1,ntheta	!theta bins
cdg	      do j=1,kmax4	!loop over eprime dbins
cdgCDJG Nadia/Aji edited
cdgc		 read(10,1004,end=499) edum,eprad2(j),thrad2,xdum,q2dum,
cdgc	1	      sigv2(j),sigrad2(j),sigrel2(j),sigrqe2(j),sigrdis2(j)
cdg		 read(10,1005,end=499) edum,eprad4(j),thrad4,q2dum,xdum,sig1dum,
cdg	1	      sig2dum,sigv4(j),sigrad4(j),radcor4(j),dum1,dum1,
cdg	2	      dum2,dum3
cdg		 if(sigrad4(j).gt.0.0) then
cdg		    radcor4(j) = sigv4(j)/sigrad4(j)
cdg		 else
cdg		    radcor4(j) = 0.0
cdg		 endif
cdgc		 write(6,*) 'radcor read in 4:',i,j,eprad4(j),thrad4,radcor4(j),sigrad4(j),sigv4(j)
cdg	      enddo		! loop over eprime
cdgc	   enddo !loop over theta
cdg 699	   continue
cdg	   close(10)
cdg	endif
CCCCC END 4TH RC FILECCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C       Open table for ratio of dummy/can external radiative corrections
	do i=32,1,-1
	   if (extfile(i:i).eq.' ') filelen=i-1
	enddo

	if(filelen.gt.0) then
	   doing_ext=.true.
	   write(33,*) 'Dummy/can external radcor file in input file'
	   xfile = 'EXT/'//extfile(1:filelen)
	   write(33,*) "opening ",xfile
	   open (unit=10,file=xfile,status='old')
	   read(10,*) npmax
	   write(6,*) 'number of EXT p bins',npmax
	   do i = 1,ntheta	!theta bins
	      do j=1,npmax	!loop over eprime dbins
		 read(10,*,end=199) edum,epext(j),thext(i),q2dum,xdum,extcan(i,j),
	1	      extdummy(i,j),extrat(i,j)
c		 write(6,*) 'radcor read in EXT:',i,j,epext(j),thext(i),extcan(i,j),extdummy(i,j),extrat(i,j)
	      enddo ! loop over eprime
	   enddo !loop over theta
 199	   continue
	   close(10)
	endif

C hardwire some things here for tests
c	doing_posi=.false.
c	doing_cryo=.false.
c	doing_dummy_posi=.false.

	makedatahistos=.true.

	if(makedatahistos) then

	   call getkine(datarun(1),.false.,.false.)

	   if(abs(th_deg-40.0).lt.1.0) then
	      thdegcent=40.0
	   elseif(abs(th_deg-50.0).lt.1.0) then
	      thdegcent=50.0
	   elseif(abs(th_deg-32.0).lt.1.0) then
	      thdegcent=32.0
	   elseif(abs(th_deg-26.0).lt.1.0) then
	      thdegcent=26.0
	   elseif(abs(th_deg-22.0).lt.1.0) then
	      thdegcent=22.0
	   elseif(abs(th_deg-18.0).lt.1.0) then
	      thdegcent=18.0
	   elseif(abs(th_deg-46.0).lt.1.0) then
	      thdegcent=46.0
	   elseif(abs(th_deg-36.0).lt.1.0) then
	      thdegcent=36.0
	   elseif(abs(th_deg-29.0).lt.1.0) then
	      thdegcent=29.0
	   elseif(abs(th_deg-24.0).lt.1.0) then
	      thdegcent=24.0
	   endif	   
	   thradcent=thdegcent*3.141592654/180.0

	   epmin = pcent*(1.0+hdeltacutlo/100.0)
	   q2min = 4.0*ebeam_cor*epmin*sin(thradcent/2.0)**2
	   xbjmin = q2min/2./m_p/(ebeam_cor-epmin) 
	   ximin = 2.*xbjmin/(1.0+sqrt(1.+4.0*m_p**2*xbjmin**2/q2min))

	   epmax = pcent*(1.0+hdeltacuthi/100.0)
	   q2max = 4.0*ebeam_cor*epmax*sin(thradcent/2.0)**2
	   xbjmax = q2max/2./m_p/(ebeam_cor-epmax) 
	   ximax = 2.*xbjmax/(1.0+sqrt(1.+4.0*m_p**2*xbjmax**2/q2max))

C       Initialize histograms.	
	   call hbook_init

C       Fill data histos.
	   totcharge = 0.0
	   cleancharge = 0.0
	   totcounts = 0.0
	   cdens = 0.0
	   ccurrent = 0.0
	   cbeampos = 0.0
	   cdummypos = 0.0 ! although this is a correction to the dummy yield, the size of
	                   ! the correction comes from the beam position on the cryo target
	   do run=1,numdatruns
	      flag = .FALSE.	!doing dummy?
	      pflag = .FALSE.	!doing positrons?
	      charge=0.0
	      charge_elcl=0.0
	      counts=0.0
	      call makedathistos(datarun(run),flag,pflag,counts,charge,charge_elcl)
	      totcounts = totcounts+counts
	      totcharge = totcharge+charge
	      cleancharge = cleancharge+charge_elcl
	      cdens = cdens + charge*denscor
	      ccurrent = ccurrent + charge*currentcor
	      cbeampos = cbeampos + charge*beamposcor
	      cdummypos = cdummypos + charge*dummyposcor
	   enddo
C       Charge weighted average correction
	   cdens = cdens/totcharge
	   ccurrent = ccurrent/totcharge
	   cbeampos = cbeampos/totcharge
	   cdummypos = cdummypos/totcharge
	   write(33,*) 'ELECTRONS'
	   write(33,*) 'Charge weighted average target density correction:',cdens
	   write(33,*) 'Charge weighted average target boiling correction:',ccurrent
	   write(33,*) 'Charge weighted average target effective length correction:',cbeampos
	   write(33,*) 'Charge weighted average correction to dummy from beam pos.:',cdummypos
	   
	   cryocor = cdens*ccurrent*cbeampos

cdg	   cryocor = ccurrent*cbeampos

	   
	   if(doing_posi) then
	      posicharge = 0.0
	      posicounts = 0.0
	      cdens = 0.0
	      ccurrent = 0.0
	      cbeampos = 0.0
	      posicdummypos = 0.0
	      do run=1,numposiruns
		 flag = .FALSE.	!doing dummy?
		 pflag = .TRUE.	!doing positrons?
		 charge=0.0
		 charge_elcl=0.0
		 counts=0.0
		 call makedathistos(posirun(run),flag,pflag,counts,charge,charge_elcl)
		 posicounts = posicounts+counts
		 posicharge = posicharge+charge_elcl
		 cdens = cdens + charge_elcl*denscor
		 ccurrent = ccurrent + charge_elcl*currentcor
		 cbeampos = cbeampos + charge_elcl*beamposcor
		 posicdummypos = posicdummypos + charge_elcl*dummyposcor
	      enddo
	      
	      cdens = cdens/posicharge
	      ccurrent = ccurrent/posicharge
	      cbeampos = cbeampos/posicharge
	      posicdummypos = posicdummypos/posicharge
	      
	      write(33,*) 'POSITRONS'
	      write(33,*) 'Charge weighted average target density correction:',cdens
	      write(33,*) 'Charge weighted average target boiling correction:',ccurrent
	      write(33,*) 'Charge weighted average target effective length correction:',cbeampos
	      write(33,*) 'Charge weighted average correction to dummy from beam pos.:',posicdummypos
	      
	      posicryocor = cdens*ccurrent*cbeampos	      
ctmp	      posicryocor =ccurrent
	   endif


C       Fill dummy histos if it's a cryo target
	   if(doing_cryo) then
	      
	      dummycharge = 0.
	      cleandummycharge=0.0
	      dummycounts =0
	      do run=1,numdumruns
		 flag = .TRUE.
		 pflag = .FALSE.
		 charge = 0.0
		 charge_elcl=0.0
		 counts = 0.0
		 call makedathistos(dummyrun(run),flag,pflag,counts,charge,charge_elcl)
		 dummycounts = dummycounts+counts
		 dummycharge=dummycharge+charge
		 cleandummycharge=cleandummycharge+charge_elcl
	      enddo
C       If we have positrons, do the dummy-positrons
	      if(doing_dummy_posi) then
		 dummyposicharge = 0.
		 dummyposicounts =0
		 do run=1,numdumposiruns
		    flag = .TRUE.
		    pflag = .TRUE.
		    charge = 0.0
		    charge_elcl = 0.0
		    counts = 0.0
		    call makedathistos(dumposirun(run),flag,pflag,counts,charge,charge_elcl)
		    dummyposicounts = dummyposicounts+counts
		    dummyposicharge=dummyposicharge+charge_elcl
		 enddo
	      endif		!doing positrons
	      
	   endif		! doing dummy
	   
C       Subtract positrons from data
	   do i = 1,nvar
	      if(totcharge.gt.0.0) then
		 dataweight = 1.0/(totcharge*cryocor)
	      else
		 dataweight = 0.0
	      endif
	      posiweight=0.0
	      if(doing_posi.and.(posicharge.gt.0.0)) then
		 posiweight = 1.0/(posicharge*posicryocor)*cleancharge/totcharge
	      endif
	      dataid = 3000+(i-1)*100+1 !e- elreal 
	      dumid = 3000+(i-1)*100+2  !e- elclean
	      newid = 3000+(i-1)*100+3  !e- elreal/elclean
C       First calculate ratio of DATA(ELREAL)/DATA(ELCLEAN) - use binomial errors for fun!
	      call hopera(dataid,'/ B',dumid,newid,one,one)
	      dataid = 3000+(i-1)*100+4	!positrons
	      dumid = 3000+(i-1)*100+3 !ratio of e- elreal/elclean
	      newid = 3000+(i-1)*100+4 !pipe it back into positrons
C       Now multiply positrons by the above ratio - it's like an efficiency! 
	      call hopera(dataid,'* E',dumid,newid,one,one)
	      dataid = 3000+(i-1)*100+1	!electrons
	      dumid = 3000+(i-1)*100+4 !positrons
	      newid = 3000+(i-1)*100+5 !positron subtracted yield
C       Now subtract positrons - we're in counts/mC now 
	      call hopera(dataid,'- E',dumid,newid,sngl(dataweight),sngl(posiweight))
	   enddo
	   
C       Subtract positrons from dummy
	   do i = 1,nvar
	      if(dummycharge.gt.0.0) then
		 dataweight = 1.0/(dummycharge*cryocor/cdummypos)
	      else
		 dataweight = 0.0
	      endif
	      posiweight=0.0
	      if(doing_dummy_posi.and.(dummyposicharge.gt.0.0)) then
		 if(.not.doing_posi) then
		    write(6,*) 'Danger Will Robinson! Subtracting positrons from dummy'
		    write(6,*) 'but NOT from the cryotarget!!!'
		    posiweight = 1.0/(dummyposicharge)*cleandummycharge/dummycharge
		 else
		    posiweight = 1.0/(dummyposicharge*posicryocor/posicdummypos)*cleandummycharge/dummycharge
		 endif
	      endif
	      dataid = 3000+(i-1)*100+6 !e- elreal 
	      dumid = 3000+(i-1)*100+7 !e- elclean
	      newid = 3000+(i-1)*100+8 !e- elreal/elclean
C       First calculate ratio of DATA(ELREAL)/DATA(ELCLEAN) - use binomial errors for fun!
	      call hopera(dataid,'/ B',dumid,newid,one,one)
	      dataid = 3000+(i-1)*100+9 !positrons
	      dumid = 3000+(i-1)*100+8  !ratio of e- elreal/elclean
	      newid = 3000+(i-1)*100+9  !pipe it back into positrons
C       Now multiply positrons by the above ratio - it's like an efficiency! 
	      call hopera(dataid,'* E',dumid,newid,one,one)
	      dataid = 3000+(i-1)*100+6 !electrons
	      dumid = 3000+(i-1)*100+9  !positrons
	      newid = 3000+(i-1)*100+10 !positron subtracted yield
C       Now subtract positrons - counts/mC now
	      call hopera(dataid,'- E',dumid,newid,sngl(dataweight),sngl(posiweight))
	   enddo
	   


C       Subtract dummy from data
	   do i = 1,nvar
c       dataweight = 1.0/totcharge
	      dataweight = 1.0
	      dummyweight=0.0
	      if(doing_cryo) then
c       dummyweight = 1/7.32/dummycharge
		 dummyweight = 1.0
	      endif
	      dataid = 3000+(i-1)*100+5
	      dumid = 3000+(i-1)*100+10
	      newid = 3000+(i-1)*100+11
	      call hopera(dataid,'- E',dumid,newid,sngl(dataweight),sngl(dummyweight))
	   enddo
	   
C       Calculate e+/e-
	   do i = 1,nvar
	      dataweight = 1.0/totcharge
	      dummyweight=0.0
	      if(doing_posi.and.posicharge.gt.0.0) then
c       dummyweight = 1/7.32/dummycharge
		 dummyweight = 1/posicharge*cleancharge/totcharge
	      endif
c	      dataid = 3000+(i-1)*100+2 !crap - this should be +1!!!
	      dataid = 3000+(i-1)*100+1 !crap - this should be +1!!!
	      dumid = 3000+(i-1)*100+4
	      newid = 3000+(i-1)*100+15
	      call hopera(dumid,'/ E',dataid,newid,sngl(dummyweight),sngl(dataweight))
	   enddo
	   
	
	   totalcounts = hsum(3005)
	   write(33,*) 'Total yield is',totalcounts,'counts/mC'
	   write(33,*) 'Total charge is',totcharge
	   write(33,*) 'Dummy charge is',dummycharge

	else  !test on making data histos
c	   write(6,*) 'opening existing hbook file: istat is',istat,hbookfile
	   lrecl=0
	   call hlimit(nwpawc)
	   call hropen(2,'ABC',hbookfile,' ',lrecl,istat)
c	   write(6,*) 'after hropen',istat
	   call hrin(0,999999,0)
	   call hrend('ABC')
	   close(2)
	endif

	call getkine(datarun(1),.false.,.false.)	
c	write(6,*) 'calling makesimhistos'
	call makesimhistos(simfile)

c	write(6,*) 'back from makesimhistos'


c	call redo_sim_histos(model)
c	open(unit=20,file='fit.out',status='unknown')
c	if(niterate.eq.0)then
c	   nitvar=1
c	   goto 666
c	endif


C       Now start iterating.
c	do i=1,niterate
C       Do W first.
C       Divide data/simc.
c	   do j=1,nvar
c	      dataid = 1000+j*100+7
c	      simid = 1000+j*100+8
c	      newid = 1000+j*100+9
c	      call hopera(dataid,'/ E',simid,newid,1.,1.)
c	   enddo
c	   call hfithn(1109,'P0','Q',1,dumpar,' ',' ',' ',sigwpar,dumchi2(1))
c	   if(nparW.ne.0) then
c	      totalsim = 0.0
c	      call hbfun1 (1111,'wfunc',nchw,wmin,wmax,wfunc) !Fill func histo
c	      call hpake(1111,wferr) !set errors to zero
c	      call hopera (1111,'* E',1109,1110,1.0,1.0)
c	      call hfith (1110,wfunc,'Q',nparW,wpar,' ',' ',' ',sigwpar,
c	1	   wchi2)
c	      call redo_sim_histos(model)
c	      if(.not.doing_hydpi) call redo_central(model)
c	   endif
c	   totalsim = hsum(1608)
c	   simtotalw = hsum(1108)
c	   if(.not.doing_hydpi) then
c	      centralsum = hsum(999)
c	      cent_fac = centralsum/central_orig
c	   else
c	      cent_fac = corfunc(xdummy)
c	   endif
c	   xsec(1) = fac*sigcentral*totalcounts/totalsim*cent_fac
c	   xsec_cut(1) = fac*sigcentral*dattotalw/simtotalw*cent_fac
c	   if(doing_hepall) xsec_cut(1) = fac*hsum(998)*totalcounts/totalsim
c	   write(20,2) 'W',xsec_cut(1),xsec(1),dumchi2(1),totalcounts/totalsim,
c	1	sigcentral*cent_fac,wpar(1),wpar(2),wpar(3),
c	2	cent_fac

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       DJG Dump some final thoughts ...
        meanchi2 = 0.
        cutratio = 0.
 666    continue
c	call makefinal(simfile,model)
	meanchisq = 0.
	do i=1,nvar
	   dataid = 3000 + (i-1)*100 + 11
	   simid = 3000 + (i-1)*100 + 12
	   newid = 3000 + (i-1)*100 + 13
c	   write(6,*) dataid,simid,newid
	   call hopera(dataid,'/ E',simid,newid,one,one)
	   chisquare(i) = 0.0
	   nchannel(i) = 0.0
	   call hgive(dataid,title,nx,xmin,xmax,ny,ymin,ymax,nwt,loc)
	   do j=1,nx
	      x1 = hi(dataid,j)
	      x2 = hi(simid,j)
	      x1err = hie(dataid,j)
	      x2err = hie(simid,j)
	      if((x1.gt.0.0.or.x2.gt.0.0).and.(x1err.gt.0.0.or.x2err.gt.0.0)) then
		 chisquare(i) = chisquare(i) + 
	1	      (x1-x2)*(x1-x2)/(x1err**2+x2err**2)
		 nchannel(i) = nchannel(i) + 1.0
	      endif
c	      write(6,*) 'x1,x2',x1,x2
	   enddo
	   chisquare(i) = chisquare(i)/(nchannel(i)-1.0)
	   write(33,*) 'Histo, Chisq, nchannel',dataid,chisquare(i),nchannel(i)
	   meanchisq = meanchisq + chisquare(i)
	enddo
	meanchisq = meanchisq/13.0

	write(33,*) 'Mean chisq', meanchisq

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	call hrput(0,hbookfile,'N')

C OK - now try to extract some cross sections
	call calc_xsec_from_delta
c	call calc_xsec_from_eprime
	call calc_xsec_from_x
c	call calc_xsec_from_xi
c	call calc_xsec_2d
C close output files
	close(33)




 1	format(a5,2x,8(a5,2x),a5)
 2	format(a5,2x,e9.5e1,2x,e9.5e1,2x,f5.3,2x,f6.4,1x,e8.4e1,1x,
	1  f6.2,1x
	2 ,f7.2,1x,f6.2,1x,f6.3)
 3	format(a5,2x,a9,2x,a9,2x,a5,2x,a6,1x,a8,1x,a6,1x,a7,1x,a6,1x,a6)
 1003	format(f6.3,4f7.3,9e12.5)      
 1004   format(1x,5f9.4,9e13.5)
 1005   format(1f6.3,4f9.4,9e13.5)
 1009   format(1x,5f9.4,5e13.5)
      end

















