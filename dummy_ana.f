
	program dummy_ana

C+______________________________________________________________________________
!
! DESCRIPTION:
!
!   The analysis program to analyze dummy target data from XEM experiment.
!   To compile, do gmak.
!
! ARGUMENTS: (All R*4 data type)
!
!   RUNNNO   :   Run number
!   DUMMYFLAG:   dummy target analysis YES/NO
!   POSIFLAG :   Positron analysis YES/NO
!
! OUTPUTS
!   COUNTS   : NO of events in real peak
!   CHARGE   : CHARGE yield
!   CHARGE_ELCL : ?
!
!
! 02-03-2015 Buddhini Waidyawansa
!            Coppied from singles_ana.f to analyze dummy target data
!
C-______________________________________________________________________________

	
	implicit none

	include 'hbook.inc'
	include 'cuts.inc'
	include 'radcor.cmn'
	include 'histos.cmn'
	include 'cryocor.cmn'
	include 'math_physics.inc'
	include 'kinematics.cmn'
	
	character*80 datafile,dummyfile,posifile,dumposifile
	character*80 xfile,simfile,hbookfile,radfile
	character*80 line,title
	integer i,j,datarun(200),dummyrun(200),dummyflag(200)
	integer posirun(200)
	integer numdatruns,numdumruns,numposiruns,run,filelen
	integer dataid,dumid,newid,simid
	integer ncomment
	real*8 dataweight, dummyweight, posiweight
	real*8 counts,charge,charge_elcl
	real*8 totcounts,totcharge,cleancharge
	real*8 posicharge,posicounts
	real*8 totalcounts
	real*8 cdummypos,posicdummypos
	real*8 cdens,ccurrent,cbeampos,cryocor,posicryocor
	real*8 hsum
	real*8 edum,q2dum,xdum,sig1dum,sig2dum,dum1,dum2,dum3
	real*8 q2min,q2max,epmin,epmax

	
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C       Open input file - read in file names and parameters
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        print *, "###################################"  
        print *, "### Open input file dummy_ana.inp"
        print *, "###################################"  
	open (unit=12,file='dummy_ana.inp',status='old')
	read(12,'(a)') datafile     !data runs
	write(6,*)" - datafile = ",datafile
	read(12,'(a)') dummyfile    !dummy runs
	write(6,*)" - dummy file = ",dummyfile
	read(12,'(a)') posifile     !positron runs
	write(6,*)" - positron file = ",posifile
	read(12,'(a)') dumposifile !dummy positron runs
	write(6,*)" - dummy positron file = ",dumposifile
	read(12,'(a)') simfile	!Monte Carlo file
	write(6,*)" - simfile = ",simfile
	read(12,'(a)') radfile      !radiative corrections table file
	write(6,*)" - radfile = ",radfile
	close (12)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Flags to identify the analysis type 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	doing_solid=.FALSE.
	doing_cryo=.FALSE.
	doing_posi=.FALSE.
	doing_dummy_posi=.FALSE.


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C       Open run list - read in run numbers
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        print *, "Get datafile name (after removing the trailing characters)"
	do i=32,1,-1
	   if (datafile(i:i).eq.' ') filelen=i-1
	enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Open efficiency file..
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        print *, "Open detector efficiency file.."
	xfile = 'OUT/'//datafile(1:filelen-4)//'.eff'
	write(6,*) " - opening ",xfile
	open (unit=11,file=xfile,status='unknown')
	write(11,1) 'RUN  ', 'HFID ','HFIDE','CO_LT','HELEC','HTRIG','HCERE',
	1  'HCALE','TrhoC','BO_Q'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Open output file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        print *, "Open an output file to store the analyisis script print output.."
	xfile = 'OUT/'//datafile(1:filelen-4)//'.out'
	write(6,*) " - opening ",xfile
	open (unit=33,file=xfile,status='unknown')

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Open files for delta-binned cross sections
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        print *, "Open files for delta cross sections.."
	xfile = 'OUT/'//datafile(1:filelen-4)//'.deltabin.xsec'
	write(6,*) " - opening ",xfile
	open (unit=44,file=xfile,status='unknown')


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Open files for x-binned cross sections
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        print *, "Open files for x-binned cross sections.."
	xfile = 'OUT/'//datafile(1:filelen-4)//'.xbin.xsec'
	write(6,*) " - opening ",xfile
	open (unit=46,file=xfile,status='unknown')


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c test if we are doing dummy stuff as target
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	xfile = datafile(12:19)

	if(INDEX(xfile,'dummy_us').gt.0) then
	   ytar = -1.0
	elseif (INDEX(xfile,'dummy_ds').gt.0) then
	   ytar = 1.0
	endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Define hbook file
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        print *, "Open an hbook file to store results.."
	hbookfile = 'hbook/'//datafile(1:filelen-4)//'.hbook'
	write(6,*) " - opening ",hbookfile
	write(33,*) " - opening ",hbookfile

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Load the run numbers of the data files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        print *, "###################################" 
        print *, "Now get the list of data runs"
        print *, "###################################"
	xfile = 'RUNS/'//datafile(1:filelen)
	write(6,*) "opening ",xfile
	open (unit=10,file=xfile,status='old')
	do i = 1,200
	   read(10,*,end=100) datarun(i)
	   write(6,*) " -- ",datarun(i)
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Open dummy run list (if there is a dummy file given)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	do i=32,1,-1
	   if (dummyfile(i:i).eq.' ') filelen=i-1
	enddo
	if(filelen.gt.0) then
	   doing_cryo=.true.
	   xfile = 'RUNS/'//dummyfile(1:filelen)
	   write(33,*) "opening ",xfile
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Open positron run list.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	do i=32,1,-1
	   if (posifile(i:i).eq.' ') filelen=i-1
	enddo
	if(filelen.gt.0) then
	   print *, "###################################" 
	   print *, "Get the list of positron runs"
	   print *, "###################################" 
	   doing_posi=.true.
	   xfile = 'RUNS/'//posifile(1:filelen)
	   write(6,*) "opening ",xfile
	   write(33,*) "opening ",xfile
	   open (unit=10,file=xfile,status='old')
	   do i = 1,200
	      read(10,*,end=300) posirun(i)
	      dummyflag(i) = 1
	      write(6,*) " -- ",posirun(i)
	   enddo
 300	   continue
	   numposiruns = i-1
	   close(10)
	   if(numposiruns.lt.1) then
	      write(33,*) 'No positron runs even though you gave me a runlist!'
	      write(6,*) 'No positron runs even though you gave me a runlist!'
	      stop
	   endif
	endif
	if(numposiruns.lt.1) then
	   doing_posi=.FALSE.
	endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                         RC FILE 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c USE this for bin centering (2D grid)
c       Open radiative corrections table
	do i=32,1,-1
	   if (radfile(i:i).eq.' ') filelen=i-1
	enddo

	if(filelen.gt.0) then
	   doing_rad=.true.
	   print *, "###################################" 
	   print *, "Read in radiative corrections file"
	   print *, "###################################" 

	   write(33,*) 'Radiative corrections file in input file'
	   xfile = 'RAD/'//radfile(1:filelen)
	   write(33,*) "opening ",xfile
	   write(6,*) " - opening 1st rc file ",xfile
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
	   ntheta = 69          ! I counted and checked this in the radcor file-buddhini
	   do i = 1,ntheta	!theta bins
	      do j=1,kmax	!loop over eprime dbins
		 read(10,1005,end=99) edum,eprad(j),thrad(i),q2dum,xdum,sigdis(i,j),
	1	      sig2dum,sigv(i,j),sigrad(i,j),radcor(i,j),coulcor(i,j),dum1,
	2	      dum2,dum3

		 if(radcor(i,j).eq.10.0) then
		    radcor(i,j) = 1.0E30
		 endif
c		 write(6,*) 'radcor read in:',i,j,eprad(j),thrad(i),radcor(i,j),sigrad(i,j),sigdis(i,j),coulcor(i,j)
	      enddo		! loop over eprime
	   enddo !loop over theta
 99	   continue
	   close(10)
	endif
	write(6,*) 'done w/first rc file'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                     END RC FILE 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Calculate kinematics and store in histograms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	print *, "###################################" 
	print *, "Calculate kinematics and store in histograms"
	print *, "###################################" 

	makedatahistos=.true.

	if(makedatahistos) then
	   print *, "Make data histograms"
	   print *, "###################################"
	   print *, "First calculate the basic kinematics using the subroutine getkine() and the first run in the list"
	   write(6,*) " - run number ",datarun(1)
	   call getkine(datarun(1),.false.,.false.)

c Subroutine getkine(runno,dummyflag,posiflag) is defined in xem/getkine.f file
c For dummy target analysis set the dummyflag to false and positron flag to false

c Assign the central angle in theta
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


c The different inputs here were all calculated in the getkine() subroutine

	   epmin = pcent*(1.0+hdeltacutlo/100.0)
	   q2min = 4.0*ebeam_cor*epmin*sin(thradcent/2.0)**2
	   xbjmin = q2min/2./m_p/(ebeam_cor-epmin) 
	   ximin = 2.*xbjmin/(1.0+sqrt(1.+4.0*m_p**2*xbjmin**2/q2min))

	   epmax = pcent*(1.0+hdeltacuthi/100.0)
	   q2max = 4.0*ebeam_cor*epmax*sin(thradcent/2.0)**2
	   xbjmax = q2max/2./m_p/(ebeam_cor-epmax) 
	   ximax = 2.*xbjmax/(1.0+sqrt(1.+4.0*m_p**2*xbjmax**2/q2max))

	   write(6,*) " - central theta = ",thdegcent
	   write(6,*) " - epmin = ",epmin
	   write(6,*) " - epmax = ",epmax
	   write(6,*) " - q2min = ",q2min
	   write(6,*) " - q2max = ",q2max
	   write(6,*) " - xbjmin = ",xbjmin
	   write(6,*) " - xbjmax = ",xbjmax
	   write(6,*) " - ximin = ",ximin
	   write(6,*) " - ximax = ",ximax


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C       Initialize histograms.	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	print *, "####"
        print *, "Now initialize histograms.."	
	   call hbook_init

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C       Fill data histos.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        print *, "Fill data histos.."	
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
	      write(6,*) 'calculate and store histos for run ',datarun(run)
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
	print *, "###################################" 
        print *, "Charge weighted average correction.."
	   cdens = cdens/totcharge
	   ccurrent = ccurrent/totcharge
	   cbeampos = cbeampos/totcharge
	   cdummypos = cdummypos/totcharge
	   print *, "####"
	   write(33,*) 'ELECTRONS'
	   write(33,*) 'Charge weighted average target density correction:',cdens
	   write(33,*) 'Charge weighted average target boiling correction:',ccurrent
	   write(33,*) 'Charge weighted average target effective length correction:',cbeampos
	   write(33,*) 'Charge weighted average correction to dummy from beam pos.:',cdummypos
	   write(6,*) 'ELECTRONS'
	   write(6,*) 'Charge weighted average target density correction:',cdens
	   write(6,*) 'Charge weighted average target boiling correction:',ccurrent
	   write(6,*) 'Charge weighted average target effective length correction:',cbeampos
	   write(6,*) 'Charge weighted average correction to dummy from beam pos.:',cdummypos
	   print *, "###################################" 	

	   cryocor = cdens*ccurrent*cbeampos
	   
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
		 write(6,*) 'calculate and store histos for positron run ',posirun(run)
		 call makedathistos(posirun(run),flag,pflag,counts,charge,charge_elcl)
		 posicounts = posicounts+counts
		 posicharge = posicharge+charge_elcl
		 cdens = cdens + charge_elcl*denscor
		 ccurrent = ccurrent + charge_elcl*currentcor
		 cbeampos = cbeampos + charge_elcl*beamposcor
		 posicdummypos = posicdummypos + charge_elcl*dummyposcor
	      enddo
	      print *, "###################################" 
	      print *, "Charge weighted average correction.."
	      cdens = cdens/posicharge
	      ccurrent = ccurrent/posicharge
	      cbeampos = cbeampos/posicharge
	      posicdummypos = posicdummypos/posicharge
	      print *, "####"
	      write(33,*) 'POSITRONS'
	      write(33,*) 'Charge weighted average target density correction:',cdens
	      write(33,*) 'Charge weighted average target boiling correction:',ccurrent
	      write(33,*) 'Charge weighted average target effective length correction:',cbeampos
	      write(33,*) 'Charge weighted average correction to dummy from beam pos.:',posicdummypos
	      write(6,*) 'POSITRONS'
	      write(6,*) 'Charge weighted average target density correction:',cdens
	      write(6,*) 'Charge weighted average target boiling correction:',ccurrent
	      write(6,*) 'Charge weighted average target effective length correction:',cbeampos
	      write(6,*) 'Charge weighted average correction to dummy from beam pos.:',posicdummypos

	      posicryocor = cdens*ccurrent*cbeampos	      
	   endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Subtract positrons from data to elliminate charge symmetric background
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	print *, "###################################" 
	print *, "Subtract positrons from data to elliminate charge symmetric background"
	print *, "###################################" 
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
	   

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	   
C       Calculate e+/e-
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	print *, "First calculate e+/e-"
	print *, "###################################" 
	   do i = 1,nvar
	      dataweight = 1.0/totcharge
	      dummyweight=0.0
	      if(doing_posi.and.posicharge.gt.0.0) then
		 dummyweight = 1/posicharge*cleancharge/totcharge
	      endif
	      dataid = 3000+(i-1)*100+1
	      dumid = 3000+(i-1)*100+4
	      newid = 3000+(i-1)*100+15
	      call hopera(dumid,'/ E',dataid,newid,sngl(dummyweight),sngl(dataweight))
	   enddo
	   
	
	   totalcounts = hsum(3005)
	   write(33,*) 'Total yield is',totalcounts,'counts/mC'
	   write(33,*) 'Total charge is',totcharge,'C'
	   write(6,*) 'Total yield is',totalcounts,'counts/mC'
	   write(6,*) 'Total charge is',totcharge,'C'
c end of if started with (makehistos)
	else !test on making data histos
	   write(6,*) 'opening existing hbook file:',hbookfile 
	   lrecl=0 !Record length in machine words. If LREC=0 the actual record length is returned on exit.
	   call hlimit(nwpawc)
c Open a direct access HBOOK file. If several direct access files are opened, they are identified by the top directory only. 
	   call hropen(2,'ABC',hbookfile,'',lrecl,istat)
c ISTAT = Return code. ISTAT=0 indicates success.
	   if(istat.eq.0) then  
	      write(6,*) 'hbook file is open!'
	   else
	      write(6,*) 'Failed to open hbook file!'
	      stop
	   endif
c Read a histogram from the current directory on the direct access file into the current directory in memory.
c HRIN (ID,ICYCLE,IOFSET)
c ID=0 means scratch all histograms in the current directory.
c 999999 is the cycle number
c Character variable specifying options selected. 
	   call hrin(0,999999,0)
c Closes a direct access file if is not needed anymore or when the options 'N' or 'U' are specified in HRFILE.
	   call hrend('ABC')
	   close(2)
	endif
	print *, "Finished making data histograms!"
	print *, "###################################"
	
	call getkine(datarun(1),.false.,.false.)	
	write(6,*) 'calling makesimhistos'
c Looks inside the simtuples directory for the 'simfile'
	call makesimhistos(simfile)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Calculate the DATA/SIM ratio and store them in new histograms
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        meanchi2 = 0.
        cutratio = 0.
 666    continue
	meanchisq = 0.
	do i=1,nvar

c For the dummy Al target, since we are not doing any dummy bkg subtractions,
c what we should access to get corrected data is
c id = 1x05 ( DATA*(1-POSITRONS/DATA(ELCLEAN)) = DATA_POSCOR OR charge normalized yeild.)

	   dataid = 3000 + (i-1)*100 + 5
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

	call hrput(0,hbookfile,'N')

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c OK - now try to extract some cross sections
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	print *, "###################################" 
	write(6,*) 'calculate xsections from delta!'
	print *, "###################################" 
	call calc_xsec_from_delta
	print *, "###################################" 
	write(6,*) 'calculate xsections from delta!'
	print *, "###################################" 
	call calc_xsec_from_x

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C close output files
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

	print *, "###################################" 
	print *, "Done!"
	print *, "###################################" 
      end

















