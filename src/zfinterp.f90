SUBROUTINE ZFINTERP(zpos,fbhbpos,sbsspos,dellpos,deltpos,spec,lbol,mass,tpos)

  !Linearly interpolate a grid of SSPs over 
  ! metallicity (zpos)
  ! fbhb (fbhbpos)
  ! sbss (sbsspos)
  ! dell (dellpos)
  ! delt (deltpos)

  !With the option to:
  !1) also interpolate age (tpos) for a single age
  !2) output over a grid of ages

  USE sps_vars
  USE sps_utils, ONLY : locate, tsum, ndinterpolate
  IMPLICIT NONE

  REAL(SP),INTENT(in) :: zpos, fbhbpos, sbsspos, dellpos, deltpos 
  REAL(SP),INTENT(in),OPTIONAL :: tpos
  REAL(SP),INTENT(inout),DIMENSION(:) :: mass, lbol 
  REAL(SP),INTENT(inout),DIMENSION(:,:) :: spec

  INTEGER  :: zlo,fbhblo,sbsslo,delllo,deltlo,tlo,i,j
  REAL(SP) :: dz,dfbhb,dsbss,ddell,ddelt,dt,z0,imdf,w1=0.25,w2=0.5,w3=0.25
  REAL(SP), DIMENSION(nz) :: mdf

  INTEGER, DIMENSION(6), PARAMETER :: nent6d=(/ntfull,nz,nfbhb,nsbss,ndell,ndelt/)
  REAL(SP), DIMENSION(6) :: pos6d
  REAL(SP), DIMENSION(SUM(nent6d)) :: grid6d

  INTEGER, DIMENSION(5), PARAMETER :: nent5d=(/nz,nfbhb,nsbss,ndell,ndelt/)
  REAL(SP), DIMENSION(5) :: pos5d
  REAL(SP), DIMENSION(SUM(nent5d)) :: grid5d

  !------------------------------------------------------------!


  !interpolate to a single metallicity and a single time
  IF (PRESENT(tpos)) THEN

     IF (SIZE(mass).GT.1.OR.SIZE(lbol).GT.1.OR.SIZE(spec(:,1)).GT.nspec.OR.SIZE(spec(1,:)).GT.1) THEN
        WRITE(*,*) 'ZFINERP ERROR: you specified an age but are '//&
             'asking for the full age array as output!'
        STOP
     ENDIF

     pos6d  = [tpos, zpos, fbhbpos, sbsspos, dellpos, deltpos]
     grid6d = [time_full, LOG10(zlegend/zsol), fbhb_legend, sbss_legend, dell_legend, delt_legend]

     ! 1) age index & differential
     tlo = MAX(MIN(locate(time_full,tpos),ntfull-1),1)
     dt  = (tpos - time_full(tlo)) / (time_full(tlo+1) - time_full(tlo))

     ! 2) metallicity index & differential
     zlo = MAX(MIN(locate(LOG10(zlegend/zsol),zpos),nz-1),1)
     dz  = (zpos-LOG10(zlegend(zlo)/zsol)) / &
          ( LOG10(zlegend(zlo+1)/zsol) - LOG10(zlegend(zlo)/zsol) )

     ! 3) fbhb index & differential
     fbhblo = MAX(MIN(locate(fbhb_legend,fbhbpos),nfbhb-1),1)
     dfbhb  = (fbhbpos - fbhb_legend(fbhblo)) / (fbhb_legend(fbhblo+1) - fbhb_legend(fbhblo))

     ! 4) sbss index & differential
     sbsslo = MAX(MIN(locate(sbss_legend,sbsspos),nsbss-1),1)
     dsbss  = (sbsspos - sbss_legend(sbsslo)) / (sbss_legend(sbsslo+1) - sbss_legend(sbsslo))

     ! 5) dell index & differential
     delllo = MAX(MIN(locate(dell_legend,dellpos),ndell-1),1)
     ddell  = (dellpos - dell_legend(delllo)) / (dell_legend(delllo+1) - dell_legend(delllo))

     ! 6) dell index & differential
     deltlo = MAX(MIN(locate(delt_legend,deltpos),ndelt-1),1)
     ddelt  = (deltpos - delt_legend(deltlo)) / (delt_legend(deltlo+1) - delt_legend(deltlo))

     ! do the n-dimensional interpolation with the fint function (must first reshape arrays to be 1D)
     mass = ndinterpolate(6, pos6d, nent6d, grid6d, mass_ssp_itp_flat)
     lbol = ndinterpolate(6, pos6d, nent6d, grid6d, lbol_ssp_itp_flat)
      
     DO i=1,nspec 
       spec(i,1) = ndinterpolate(6, pos6d, nent6d, grid6d, spec_ssp_itp_flat(i,:))
     ENDDO
     
  ELSE

     !interpolate to a single metallicity and return a grid of ages
     pos5d  = [zpos, fbhbpos, sbsspos, dellpos, deltpos]
     grid5d = [LOG10(zlegend/zsol), fbhb_legend, sbss_legend, dell_legend, delt_legend]

      ! 1) metallicity index & differential
      zlo = MAX(MIN(locate(LOG10(zlegend/zsol),zpos),nz-1),1)
      dz  = (zpos-LOG10(zlegend(zlo)/zsol)) / &
            ( LOG10(zlegend(zlo+1)/zsol) - LOG10(zlegend(zlo)/zsol) )

      ! 2) fbhb index & differential
      fbhblo = MAX(MIN(locate(fbhb_legend,fbhbpos),nfbhb-1),1)
      dfbhb  = (fbhbpos - fbhb_legend(fbhblo)) / (fbhb_legend(fbhblo+1) - fbhb_legend(fbhblo))

      ! 3) sbss index & differential
      sbsslo = MAX(MIN(locate(sbss_legend,sbsspos),nsbss-1),1)
      dsbss  = (sbsspos - sbss_legend(sbsslo)) / (sbss_legend(sbsslo+1) - sbss_legend(sbsslo))

      ! 4) dell index & differential
      delllo = MAX(MIN(locate(dell_legend,dellpos),ndell-1),1)
      ddell  = (dellpos - dell_legend(delllo)) / (dell_legend(delllo+1) - dell_legend(delllo))

      ! 5) dell index & differential
      deltlo = MAX(MIN(locate(delt_legend,deltpos),ndelt-1),1)
      ddelt  = (deltpos - delt_legend(deltlo)) / (delt_legend(deltlo+1) - delt_legend(deltlo))
   
      DO j=1,ntfull
         mass(j) = ndinterpolate(5, pos5d, nent5d, grid5d, mass_ssp_itp(j,:))
         lbol(j) = ndinterpolate(5, pos5d, nent5d, grid5d, lbol_ssp_itp(j,:))
         DO i=1,nspec
            spec(i,j) = ndinterpolate(5, pos5d, nent5d, grid5d, spec_ssp_itp(i,j,:))
         ENDDO
      ENDDO

  ENDIF


END SUBROUTINE ZFINTERP
