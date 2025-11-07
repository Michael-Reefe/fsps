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

  INTEGER  :: zlo,fbhblo,sbsslo,delllo,deltlo,tlo
  INTEGER  :: zi,fbhbi,sbssi,delli,delti,ti
  REAL(SP) :: dz,dfbhb,dsbss,ddell,ddelt,dt,wght,w1,w2,w3,w4,w5,w6

  !------------------------------------------------------------!


  !interpolate to a single metallicity and a single time
  IF (PRESENT(tpos)) THEN

     IF (SIZE(mass).GT.1.OR.SIZE(lbol).GT.1.OR.SIZE(spec(:,1)).GT.nspec.OR.SIZE(spec(1,:)).GT.1) THEN
        WRITE(*,*) 'ZFINERP ERROR: you specified an age but are '//&
             'asking for the full age array as output!'
        STOP
     ENDIF

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

     ! do the n-dimensional interpolation over 6 axes
     mass = 0.
     lbol = 0.
     spec(:,1) = 0.
     DO ti=tlo,tlo+1
       ! weight for time axis, will be 1-dt when ti=tlo, and dt when ti=tlo+1
       w1 = merge(1.-dt, dt, ti==tlo)
       DO zi=zlo,zlo+1
         ! weight for metallicity axis
         w2 = merge(1.-dz, dz, zi==zlo)
         DO fbhbi=fbhblo,fbhblo+1
            ! weight for fbhb axis
            w3 = merge(1.-dfbhb, dfbhb, fbhbi==fbhblo)
            DO sbssi=sbsslo,sbsslo+1
               ! weight for sbss axis
               w4 = merge(1.-dsbss, dsbss, sbssi==sbsslo)
               DO delli=delllo,delllo+1
                  ! weight for dell axis
                  w5 = merge(1.-ddell, ddell, delli==delllo)
                  DO delti=deltlo,deltlo+1
                     ! weight for delt axis
                     w6 = merge(1.-ddelt, ddelt, delti==deltlo)
                     ! overall weight 
                     wght = w1*w2*w3*w4*w5*w6
                     ! add the contribution from this point, multiplied by the weight
                     mass = mass + wght * mass_ssp_itp(ti, zi, fbhbi, sbssi, delli, delti)
                     lbol = lbol + wght * lbol_ssp_itp(ti, zi, fbhbi, sbssi, delli, delti)
                     spec(:,1) = spec(:,1) + wght * spec_ssp_itp(:, ti, zi, fbhbi, sbssi, delli, delti)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
       ENDDO
     ENDDO
     
  ELSE

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

      ! do the n-dimensional interpolation over 5 axes
      mass = 0.
      lbol = 0.
      spec = 0.
      DO zi=zlo,zlo+1
        ! weight for metallicity axis
        w1 = merge(1.-dz, dz, zi==zlo)
        DO fbhbi=fbhblo,fbhblo+1
          ! weight for fbhb axis
          w2 = merge(1.-dfbhb, dfbhb, fbhbi==fbhblo)
          DO sbssi=sbsslo,sbsslo+1
            ! weight for sbss axis
            w3 = merge(1.-dsbss, dsbss, sbssi==sbsslo)
            DO delli=delllo,delllo+1
               ! weight for dell axis
               w4 = merge(1.-ddell, ddell, delli==delllo)
               DO delti=deltlo,deltlo+1
                  ! weight for delt axis
                  w5 = merge(1.-ddelt, ddelt, delti==deltlo)
                  ! overall weight
                  wght = w1*w2*w3*w4*w5
                  ! add the contribution from this point, multiplied by the weight
                  mass(:) = mass(:) + wght * mass_ssp_itp(:, zi, fbhbi, sbssi, delli, delti)
                  lbol(:) = lbol(:) + wght * lbol_ssp_itp(:, zi, fbhbi, sbssi, delli, delti)
                  spec(:,:) = spec(:,:) + wght * spec_ssp_itp(:, :, zi, fbhbi, sbssi, delli, delti)
               ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

  ENDIF

END SUBROUTINE ZFINTERP


SUBROUTINE ZFINTERP2D(fbhbpos,sbsspos,dellpos,deltpos,spec,lbol,mass)

  !Linearly interpolate a grid of SSPs over 
  ! fbhb (fbhbpos)
  ! sbss (sbsspos)
  ! dell (dellpos)
  ! delt (deltpos)

  USE sps_vars
  USE sps_utils, ONLY : locate, tsum, ndinterpolate
  IMPLICIT NONE

  REAL(SP),INTENT(in) :: fbhbpos, sbsspos, dellpos, deltpos 
  REAL(SP),INTENT(inout),DIMENSION(:,:) :: mass, lbol 
  REAL(SP),INTENT(inout),DIMENSION(:,:,:) :: spec

  INTEGER  :: fbhblo,sbsslo,delllo,deltlo
  INTEGER  :: fbhbi,sbssi,delli,delti
  REAL(SP) :: dfbhb,dsbss,ddell,ddelt,wght,w1,w2,w3,w4

  !------------------------------------------------------------!

   ! 1) fbhb index & differential
   fbhblo = MAX(MIN(locate(fbhb_legend,fbhbpos),nfbhb-1),1)
   dfbhb  = (fbhbpos - fbhb_legend(fbhblo)) / (fbhb_legend(fbhblo+1) - fbhb_legend(fbhblo))

   ! 2) sbss index & differential
   sbsslo = MAX(MIN(locate(sbss_legend,sbsspos),nsbss-1),1)
   dsbss  = (sbsspos - sbss_legend(sbsslo)) / (sbss_legend(sbsslo+1) - sbss_legend(sbsslo))

   ! 3) dell index & differential
   delllo = MAX(MIN(locate(dell_legend,dellpos),ndell-1),1)
   ddell  = (dellpos - dell_legend(delllo)) / (dell_legend(delllo+1) - dell_legend(delllo))

   ! 4) dell index & differential
   deltlo = MAX(MIN(locate(delt_legend,deltpos),ndelt-1),1)
   ddelt  = (deltpos - delt_legend(deltlo)) / (delt_legend(deltlo+1) - delt_legend(deltlo))

   ! do the n-dimensional interpolation over 5 axes
   mass = 0.
   lbol = 0.
   DO fbhbi=fbhblo,fbhblo+1
      ! weight for fbhb axis
      w1 = merge(1.-dfbhb, dfbhb, fbhbi==fbhblo)
      DO sbssi=sbsslo,sbsslo+1
         ! weight for sbss axis
         w2 = merge(1.-dsbss, dsbss, sbssi==sbsslo)
         DO delli=delllo,delllo+1
            ! weight for dell axis
            w3 = merge(1.-ddell, ddell, delli==delllo)
            DO delti=deltlo,deltlo+1
               ! weight for delt axis
               w4 = merge(1.-ddelt, ddelt, delti==deltlo)
               ! overall weight
               wght = w1*w2*w3*w4
               ! add the contribution from this point, multiplied by the weight
               mass(:,:) = mass(:,:) + wght * mass_ssp_itp(:, :, fbhbi, sbssi, delli, delti)
               lbol(:,:) = lbol(:,:) + wght * lbol_ssp_itp(:, :, fbhbi, sbssi, delli, delti)
               spec(:,:,:) = spec(:,:,:) + wght * spec_ssp_itp(:, :, :, fbhbi, sbssi, delli, delti)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE ZFINTERP2D

