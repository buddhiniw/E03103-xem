       subroutine total_eloss(angle,z,a,dens,hsp,tgthick,e_loss)

*------------------------------------------------------------------------------
*-         Prototype C routine
*- 
*-
*-    Purpose and Method :  In separate calls, calculate the energy loss for 
*-                          the incident electron in the target OR the energy  
*-                          loss for exiting particles in the target and 
*-                          other materials like windows. Cryogenic targets
*-                          must be beer-can cells. Solid targets are okay too.
*-                          Ytarget information is NOT used; all calculations
*-                          assume the reaction vertex is at the target center.
*-
*-    Output: loss            -   energy loss for the arm requested
*-    Created   1-Dec-1995  Rolf Ent
*------------------------------------------------------------------------------
**********************
* LH2 and LD2 targets
**********************
*
* Incoming beam sees the following materials to target center:
*    1.  a 3.0 mil Al-foil (upstream endcap of target) J. Dunne Dec 96
*    2.  half the target thickness
*
* Any particle exiting target center sees the following materials:  
*    3. Particle leaves thru side-walls:1.325 inch of target material corrected
* for the spectrometer angle, OR    
*       Particle leaves thru downstream window: half the target length, correc-
* ted for the spectrometer angle.
*
*    4.  A 5.0 mil Al-foil target wall thickness (J. Dunne Dec 96), corrected 
* for spectrometer angle.
*
****************** 
* Solid targets:
****************** 
*
* Incoming beam sees the following materials to target center:
*     1.  half the target thickness, corrected for the spectrometer angle.
* 
* Any particle exiting target center sees the following materials:
*     2.  half the target thickness, corrected for the spectrometer angle
*
***************************************************
*     Additional materials (irregardless of target):
***************************************************
*     * effective density for kevlar is 0.74 
*     * effective z for CH2 is 2.67, effective a is 4.67
*                                (values confirmed by T. Keppel Mar. 98)
*
*	HMS particles only: 
*	1. 16 mil Aluminum scattering chamber window (J. Mitchell Feb. 98)
*	2. 15 cm of air between chamber window and HMS entrance window.
*          *effective a for air is 14.68, effective z is 7.32, dens is .00121
*     	3. HMS entrance window, 17 mil kevlar and 5 mil mylar. 
*                                 (values confirmed by T. Keppel Mar. 98)
*
*	SOS particles only: 
*	1. 8.0 mil Al-foil scattering chamber window (J.Mitchell Feb 98)
*	2. 15 cm of air between chamber window and HMS entrance window.
*          *effective a for air is 14.68, effective z is 7.32, dens is .00121
*  	2.  SOS entrance window, 6 mil kevlar and 1.5 mil mylar.
*                                  (values confirmed by T. Keppel Mar. 98)
     
      IMPLICIT NONE
      SAVE
*

      LOGICAL liquid

      REAL*4 crit_angle,tg_spect_angle
      REAL*4 z,a,tgthick,dens,angle,tgangle,beta,type
      REAL*4 thick,thick_side,thick_front,e_loss,total_loss
      REAL*4 targ_win_loss,front_loss,back_loss,cell_wall_loss
      REAL*4 scat_win_loss,air_loss,h_win_loss,s_win_loss
      REAL*4 gcell_radius,gz_cell,ga_cell,gcell_den,gwall_thk,gend_thk
      REAL*4 hscat_win_thk,hscat_win_den,hscat_win_a,hscat_win_z
      REAL*4 hdet_ent_thk,hdet_ent_den,hdet_ent_a,hdet_ent_z
      REAL*4 gfront_thk,gair_dens,gair_thk,gair_a,gair_z
      REAL*4 hsp
      REAL*8 beta_temp,gamma_temp,X_temp,frac_temp,p_temp
      REAL*4 velocity
      real*4 mass_electron
      parameter (mass_electron = 0.000510999)

c      write(6,*) 'cheesy poofs',angle,z,a,dens,hsp,tgthick
********************INITIALIZE ENERGY LOSS VARIABLES*****************
      e_loss         = 0.0
      total_loss     = 0.0
      targ_win_loss  = 0.0
      front_loss     = 0.0
      back_loss      = 0.0
      cell_wall_loss = 0.0
      scat_win_loss  = 0.0
      air_loss       = 0.0
      h_win_loss     = 0.0
      s_win_loss     = 0.0
      liquid =.FALSE.
***********************SETUP OF PARAMETERS****************************
* Parameters should be accessed via the various common blocks
* No more hardwired #s in the code!!! (I.N. 2001)
*
*     
* z,a,tgthick,dens come via a common block
*
!      tgthick=gtarg_thick(gtarg_num)
!      dens=gtarg_dens(gtarg_num)
!      write(*,*) 'got passed tgt thick and dens', tgthick,dens
      tgangle=1.570796327 
      type=22
      if(z.lt.5) type=1.

      gair_dens         =   0.00121
      gair_thk  =   0.018
      gair_a    =   14.68
      gair _z   =   7.32

      gcell_radius = 2.008
      gz_cell    =  13.0
      ga_cell    =  27.0
      gcell_den  =   2.7
      gwall_thk  =  0.03477
      gend_thk   =  0.03477
      gfront_thk =  0.03477   

      hscat_win_thk  = 0.109728
      hscat_win_den =  2.70
      hscat_win_z   =  13.0
      hscat_win_a   =  27.0
      hdet_ent_thk  =  0.049098
      hdet_ent_den  =  0.878636
      hdet_ent_z    =  2.67
      hdet_ent_a    =  4.67


*******DIVIDE BY ZERO CHECK**************************************
      if ((gcell_radius.eq.0.0).or.(gz_cell.eq.0.0).or.(ga_cell.eq.0.0)
     &  .or.(gcell_den.eq.0.0).or.(gwall_thk.eq.0.0).or.(gend_thk.eq.0.0
     +  ).or.(gfront_thk.eq.0.0)) then
!        write(6,*)'Total_eloss: Uninitialized target variable(s)!!!'
!        write(6,*)'gcell_radius = ',gcell_radius
!        write(6,*)'gz_cell      = ',gz_cell
!        write(6,*)'ga_cell      = ',ga_cell
!        write(6,*)'gcell_den    = ',gcell_den
!        write(6,*)'gwall_thk    = ',gwall_thk
!        write(6,*)'gend_thk     = ',gend_thk
!        write(6,*)'gfront_thk   = ',gfront_thk
        stop
      elseif ((hscat_win_den.eq.0.0).or.
     &    (hscat_win_thk.eq.0.0).or.(hscat_win_z.eq.0.0).or.
     &    (hscat_win_a.eq.0.0).or.(hdet_ent_z.eq.0.0).or.
     &    (hdet_ent_a.eq.0.0)) then
        write(6,*)'Total_eloss: Uninitialized HMS window specs!!!'
        stop
      else
      endif
      
      if ((z*a*tgthick*dens*tgangle).eq.0.0) then
         write(6,*)'Total_eloss: Uninitialized target material!!!'
         write(6,*)
         write(6,*)'target angle = ',angle
         write(6,*)'target type  = ',type
         write(6,*)'thickness    = ',tgthick
         write(6,*)'Z            = ',z
         write(6,*)'A            = ',A
         write(6,*)'density      = ',dens

         stop
         
      endif
*
* If an angle is provided, use it, otherwise use the central
* spectrometer angle
*


 10            format(7(2x,A10))
 20            format(12x,6(2x,f10.9))
 30            format(5(2x,A10))
 40            format(12x,4(2x,f10.9))
 50            format(12x,3(2x,f10.9))
 60            format(4(A12))
 70            format(10(A11))
 80            format(2x,I9,9(2x,f9.6))
***********************END SETUP******************************

*******************************************************************************
* With the adoption of a new, beta-dependent energy loss correction formula
* for electrons, it became necessary to give the velocity of electrons in terms
* of log_10(beta*gamma), since REAL*4 was not good enough to distinguish the
* beta of electrons from 1. For hadrons, nothing will change.
*******************************************************************************

      velocity=0.

      p_temp=hsp
      p_temp=max(p_temp,.1D0)
      frac_temp=mass_electron/p_temp

      beta_temp=1./sqrt(1.+frac_temp**2)
      gamma_temp=sqrt(1.+frac_temp**2)/frac_temp
      X_temp=log(beta_temp*gamma_temp)/log(10.)   
      velocity=X_temp

**************************************************************************
* Calculate the angle at which the ejectile passes through the side of the
* target cell rather than the end.
**************************************************************************

       if ((type.eq.2).and.(tgthick.ne.0.).and.(dens.ne.0.)) then
          crit_angle= atan(gcell_radius/(tgthick/dens/2))
       else
          crit_angle= 0.45
       endif

**************************************************************************
* Define hydrogen, deuterium and 3,4He as liquid targets: z<=2
**************************************************************************

       if (a.le.4) liquid =.TRUE. 
       
       
       if (liquid) then         ! cryo target !
         thick = gcell_radius*dens !!! tuna-can !!!
       else
*     
*     Assume that tgangle = 90 deg 
*     corresponds to a target normal to the beam direction
*
         if (abs(sin(tgangle)).ge.0.01) then
           thick = tgthick/2./abs(sin(tgangle))
         else
           thick = tgthick/2./0.01
         endif
       endif
*********************************************************************
*Calculate the energy loss of ejectile after the target center.
*********************************************************************
!       write(*,*) 'after target, cryo and al'
      if (liquid ) then
*     Liquid target*********
         if (type.eq.1) then   
            call loss(z,a,thick,dens,velocity,back_loss) !liquid
            total_loss = total_loss + back_loss
            call loss(gz_cell,ga_cell,gfront_thk,gcell_den !aluminum 
     >           ,velocity,cell_wall_loss)                          
            total_loss = total_loss + cell_wall_loss
            
         else
           write(*,*) 'liquid not 1, should never be here'
*     write(6,*)'********************I am HERE*****************(VT)'
            thick=0.0
            thick_front=0.0
            thick_side=0.0
*     Through the end of the cell.
            if (angle.le.crit_angle) then        
               if (cos(angle).ge.0.01) then
                  thick = abs(gend_thk/cos(angle))
                  thick_front= abs(tgthick/2./cos(angle))
               else
                  thick = abs(gend_thk/0.01)
                  thick_front= abs(tgthick/2./0.01)
               endif
               call loss(z,a,thick_front,dens,velocity,back_loss) !liquid
               total_loss = total_loss + back_loss
               call loss(gz_cell,ga_cell,thick,gcell_den,velocity !aluminum
     >              ,cell_wall_loss)                          
               total_loss = total_loss + cell_wall_loss
*     Through the side of the cell. 
            else					
               if (abs(sin(angle)).ge.0.01) then
                  thick = abs(gwall_thk/abs(sin(angle)))
                  thick_side  = abs(gcell_radius*dens/abs(sin(angle)))
               else
                  thick = abs(gwall_thk/0.01)
                  thick_side  = abs(gcell_radius*dens/0.01)
               endif
               call loss(z,a,thick_side,dens,velocity,back_loss) !liquid
               total_loss = total_loss + back_loss
               call loss(gz_cell,ga_cell,thick,gcell_den,velocity !aluminum
     >              ,cell_wall_loss)                        
               total_loss = total_loss + cell_wall_loss
             endif
!             if(total_loss.GE.1e-2)
!     &         write(*,*)  velocity, 
!     &         "total_loss ", total_loss, " back_loss ",back_loss
             
           endif
*     Solid target************
         else    
           
*     In any ordinary case, the solid target has angle of 90 degrees
*     with respect to the beam direction: tgangle=90.*degrad
           
*     csa 1/5/99 -- Here I define tgangle > 90 deg to mean that the
*     solid target is facing the SOS.
           
           tg_spect_angle = angle + tgangle
         if (abs(sin(tg_spect_angle)).ge.0.01) then
           thick = abs((tgthick/2.)/abs(sin(tg_spect_angle)))
         else
           thick = abs((tgthick/2.)/0.01)
         endif
!         write(*,*) 'calling for solid target'
         call loss(z,a,thick,dens,velocity,back_loss) !generic solid target
         total_loss = total_loss + back_loss
       endif         
************************************
*     Now calculate the HMS energy loss.  
************************************
       
*     16 mil aluminum scattering chamber window on HMS side
!       write(*,*) 'hms chamber window'
       call loss(hscat_win_z,hscat_win_a,hscat_win_thk,
     &   hscat_win_den,velocity,scat_win_loss) !aluminum
       total_loss = total_loss + scat_win_loss
       
*     ENERGY LOSS IN AIR GAP BEWTEEN THE CHAMBER AND THE ENTRANCE WINDOW
!       write(*,*) 'hms in air gap'
       call loss(gair_z,gair_a,gair_thk,gair_dens,velocity,air_loss)
     +   
       total_loss = total_loss + air_loss
!       write(*,*) 'hms entrance window'
       
*     HMS Det. entrance window loss
       call loss(hdet_ent_z,hdet_ent_a,hdet_ent_thk,
     &   hdet_ent_den,velocity,h_win_loss) !HMS window
       total_loss = total_loss + h_win_loss
       
       e_loss = total_loss
       
*     eloss debug HMS
         if(liquid) then
!           write(6,10)'liquid',
!     &       'back','cell_wall','scat_win','air','HMS_win',
!     &       'total'
!           write(6,20) back_loss,cell_wall_loss,scat_win_loss,air_loss
!     +       ,h_win_loss,total_loss
         else
!          write(6,30)'solid', 'scat_win','air','HMS_win','total'
!          write(6,40) scat_win_loss,air_loss,h_win_loss,total_loss
         endif
!         write(6,*)
       
       
       
 100   continue
       
       RETURN
       END
      
*-------------------------------------------------------------
      subroutine loss(z,a,thick,dens,velocity,e_loss)
*-------------------------------------------------------------
*-         Prototype C function 
*- 
*-
*-    Purpose and Method :  Calculate energy loss 
*-    
*-    Output: -
*-    Created   1-Dec-1995  Rolf Ent
*-   
*-    Verification:  The non-electron portion on this subr. is Bethe_Bloch
*-                   equation (Physial Review D vol.50 (1994) 1251 with full
*-		     calculation of Tmax and the density correction. The electron
*-		     part has been switched from O'Brien, Phys. Rev. C9(1974)1418,
*-		     to Bethe-Bloch with relativistic corrections and density
*-		     density correction, Leo, Techniques for Nuclear and Particle 
*-		     Physics Experiments
*-                   J. Volmer 8/2/98 16:50
*------------------------------------------------------------------------------*
      IMPLICIT NONE
      SAVE
*
*
      LOGICAL electron
      REAL*4 eloss,z,a,thick,dens,beta,e_loss
      REAL*4 icon_ev,me_ev
      REAL*4 icon_gev,me_gev
      REAL*4 particle
      REAL*4 denscorr,hnup,c0,log10bg,pmass,tmax,gamma,velocity
      REAL*4 tau,betagamma
      real*4 mass_electron
      parameter (mass_electron = 0.000510999)
      parameter (me_ev = 510999.)
      parameter (me_gev = 0.000510999)
*

 91   format(7(A10))
 90   format(7(2x,f8.5))
      e_loss = 0.0
      eloss  = 0.0
      electron=.true.
!      write (*,*) 'subroutine got thickness ', thick

*************************************************************************
* for debugging print out all variables that have been passed on tol loss
*************************************************************************

*****************************************************************************
* calculate the mean excitation potential I in a newer parametrization 
* given in W.R. Leo's Techniques for Nuclear and Particle Physics Experiments
*****************************************************************************

*     csa 1/99 -- Note that this code calculates the mean energy loss,
*     not the most probable. This is appropriate for the case (as in
*     Hall C) where the resolution of the measurement is significantly
*     greater than the energy loss.

      if (z.lt.1.5) then
         icon_ev = 21.8
      elseif (z.lt.13) then
         icon_ev = 12.*z+7.
      elseif (z.ge.13) then
         icon_ev = z*(9.76+58.8*z**(-1.19))
      endif
      icon_gev = icon_ev*1.0e-9

**********************************************
* extract the velocity of the particle:
*     hadrons:   velocity = beta
*     electrons: velocity = log_10(beta*gamma)
**********************************************

         log10bg=velocity
         betagamma=exp(velocity*log(10.))
         beta=betagamma/(sqrt(1.+betagamma**2))
         gamma=sqrt(1.+betagamma**2)
         tau=gamma-1.

******************************************************
* calculate the density correction, as given in Leo,
* with Sternheimer's parametrization
* I is the mean excitation potential of the material
* hnup= h*nu_p is the plasma frequency of the material
******************************************************

      denscorr=0.
      if(A.gt.0.) then
         HNUP=28.816E-9*sqrt(abs(DENS*Z/A))
      else
         HNUP=28.816E-9*sqrt(abs(DENS*Z/1.))
      endif

* log(icon_gev/hnup)=log(icon_gev)-log(hnup)
      C0=-2*(log(icon_gev)-log(hnup)+.5)

      if(log10bg.lt.0.) then
         denscorr=0.
      elseif(log10bg.lt.3.) then
         denscorr=C0+2*log(10.)*log10bg+abs(C0/27.)*(3.-log10bg)**3
      elseif(log10bg.lt.4.7) then
         denscorr=C0+2*log(10.)*log10bg
      else
         denscorr=C0+2*log(10.)*4.7
      endif


**********************************************************************       
* now calculate the energy loss for electrons 
**********************************************************************
* electron
      if((thick.gt.0.0).and.(dens.gt.0.0).and.(a.gt.0.).and.(beta.gt.0.)
     >  .and.(tau.gt.0).and.(betagamma.gt.0))then
        eloss=0.1535e-03*z/a*thick/beta**2*(
     >    2*log(tau)+log((tau+2.)/2.)-2*(log(icon_gev)-log(me_gev))
     >    +1-beta**2+(tau**2/8-(2*tau+1)*log(2.))/(tau+1)**2
     >    -(-(2*(log(icon_gev)-log(hnup))+1)+2*log(betagamma)))
      endif
      

      if (eloss.le.0.) write(6,*)'loss: eloss<=0!'
* units should be in GeV
      e_loss = eloss

      if ((eloss.le.0)) then
         particle=0.0
         if (electron) particle=1.0
         write(6,91) 'electron?','ztgt','atgt','thick','dens','velocity','e_loss'
         write(6,90) particle,z,a,thick,dens,velocity,e_loss
         write(6,'(4A10)') 'velocity','beta','pmass','denscorr'
         write(6,'(6(2x,f8.5))') velocity,beta,pmass,denscorr
         write(6,'(6A10)') 'betagamma','log10bg','tau','gamma','icon_ev','hnup (eV)'
         write(6,'(6(2x,F8.3))') betagamma,log10bg,tau,gamma,icon_ev,hnup*1e9
      endif

      RETURN
      END
