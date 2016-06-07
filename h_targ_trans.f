      SUBROUTINE H_TARG_TRANS(hx_fp,hxp_fp,hy_fp,hyp_fp,gbeam_y,
     >         hdelta_tar,hxp_tar,hyp_tar,hy_tar)
*--------------------------------------------------------
*-
*-   Purpose and Methods :  Transforms tracks from HMS focal plane to 
*-                          target.
*-
*-      Required Input BANKS     HMS_FOCAL_PLANE
*-
*-      Output BANKS             HMS_TARGET
*-
*-   Output: ABORT           - success or failure
*-         : err             - reason for failure, if any
*-
*-  istat   (integer) Status flag. Value returned indicates the following:
*-           = 1      Normal return.
*-           = 2      Matrix elements not initted correctly.
*
* Abstract: Reconstruct target scattering variables from track variables in
*           the detectors, using a polynomial (Taylor series) map. The track,
*           target, and map data are all maintained in common blocks.
*
* NOTE:     This version assumes that the beam is not rastered.
*           Also, there is no treatment of error matrices, yet.
*
* Output arguments:
*
*
* Right-handed coordinates are assumed: X=down, Z=downstream, Y = (Z cross X)
*
* Author:   David H. Potterveld, Argonne National Lab, Nov. 1993
*--------------------------------------------------------
      IMPLICIT NONE

*
      logical ABORT
      integer*4   istat
*
      include 'hms_recon_elements.cmn'
*
*--------------------------------------------------------
*
* Misc. variables.

      integer*4        i,j

      real*8           hx_fp,hxp_fp,hy_fp,hyp_fp,gbeam_y,hdelta_tar
      real*8           hxp_tar,hyp_tar,hy_tar,hx_tar,hz_tar
      real*8           term
      real*8           sum(4),hut(5),hut_rot(5)

      SAVE
*=============================Executable Code =============================
      ABORT= .FALSE.
c      err= ' '
* Check for correct initialization.

      if (h_recon_initted.ne.1) then
         istat = 2
         return
      endif
      istat = 1
      

      hdelta_tar=0.0
      hxp_tar=0.0
      hyp_tar=0.0
      hy_tar=0.0
*     Reset COSY sums.
      do i = 1,4
         sum(i) = 0.
      enddo
      
      
*     Load track data into local array, Converting to COSY units.
*     It is assumed that the track coordinates are reported at
*     the same focal plane as the COSY matrix elements were calculated.
*     Also note that the COSY track slopes HUT(2) and HUT(4) are actually
*     the SINE of the track angle in the XZ and YZ planes.
      
      hut(1) = hx_fp/100. + h_z_true_focus*hxp_fp
     >     + h_det_offset_x     ! include detector offset  (m)
!     includes transformation to actual focus if not at Z=0.
      
      hut(2) = hxp_fp + h_ang_offset_x !radians
      
      hut(3) = hy_fp/100. + h_z_true_focus*hyp_fp
     >     + h_det_offset_y     !m
!     again icludes transformation to true focus.
      
      hut(4) = hyp_fp + h_ang_offset_y !radians
      
      hut(5)= -gbeam_y/100.     ! spectrometer target X in meter!
                                ! note that pos. spect. X = neg. beam Y
      
!     now transform 
*     hx_fp_rot(itrk)=  hut(1) + h_det_offset_x    ! include detector offset
*     hy_fp_rot(itrk)=  hut(3) + h_det_offset_y 
*     hxp_fp_rot(itrk)= hut(2) + hut(1)*h_ang_slope_x
*     hyp_fp_rot(itrk)= hut(4) + hut(3)*h_ang_slope_y
*     hut_rot(1)= hx_fp_rot(itrk)
*     hut_rot(2)= hxp_fp_rot(itrk)
*     hut_rot(3)= hy_fp_rot(itrk)
*     hut_rot(4)= hyp_fp_rot(itrk)
*     h*_fp_rot never used except here, so remove the intermediate step.
      
      hut_rot(1) = hut(1)
      hut_rot(2) = hut(2) + hut(1)*h_ang_slope_x
      hut_rot(3) = hut(3)
      hut_rot(4) = hut(4) + hut(3)*h_ang_slope_y
      hut_rot(5) = hut(5)
*     Compute COSY sums.
      do i = 1,h_num_recon_terms
         term = 1.
         do j = 1,5
            if (h_recon_expon(j,i).ne.0.) then
               term = term*hut_rot(j)**h_recon_expon(j,i)
            endif
         enddo
         sum(1) = sum(1) + term*h_recon_coeff(1,i)
         sum(2) = sum(2) + term*h_recon_coeff(2,i)
         sum(3) = sum(3) + term*h_recon_coeff(3,i)
         sum(4) = sum(4) + term*h_recon_coeff(4,i)
      enddo
*     Protext against asin argument > 1.
c     if(sum(1).gt. 1.0) sum(1)=  0.99
c     if(sum(1).lt. -1.0) sum(1)= -.99
c     if(sum(3).gt. 1.0) sum(3)=  0.99
c     if(sum(3).lt. -1.0) sum(3)= -.99
      
      
*     Load output values.
      
      hx_tar = 0.               ! ** No beam raster yet **
      hy_tar = sum(2)*100.      !cm.
      hxp_tar = sum(1)          !Slope xp
      hyp_tar = sum(3)          !Slope yp
      
      hz_tar = 0.0              !Track is at origin
      hdelta_tar = sum(4)*100.  !percent.

* Apply offsets to reconstruction.
c         hdelta_tar = hdelta_tar + hdelta_offset
c         hyp_tar = hyp_tar + htheta_offset
c         hxp_tar = hxp_tar + hphi_offset

C for now use hardwired values

         hdelta_tar = hdelta_tar + 0.0
         hyp_tar = hyp_tar + 4.83d-5
         hxp_tar = hxp_tar - 4.73d-3


         return
         end
