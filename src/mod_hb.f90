SUBROUTINE MOD_HB(f_bhb,t,mini,mact,logl,logt,logg,phase, &
     wght,hb_wght,nmass,hbtime,macthb,loglhb,logthb,logghb,mchb, &
     yhb,zhb,timehb)

  !routine to modify the horizontal branch to include bluer 
  !stars.  The fiducial HB stars are identified and then a fraction
  !f_bhb uniformly redistributed from the red clump to 10^4 K 
  !We also want to be able to call this routine when f_bhb=0 
  !in order to count the weight of stars on the HB.

  !see e.g. Sarajedini et al. 2007.  This treatment of BHB 
  !produces spectra similar to the Maraston 2005 models
  !see also Jimenez et al. 2004

  !Note that the parameter bhb_sbs_time, set in sps_vars.f90,
  !sets the turn-on time for this modification.

  USE sps_vars
  USE sps_utils, ONLY : imf, funcint
  IMPLICIT NONE

  REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: mini,mact,&
       logl,logt,logg,phase
  REAL(SP), INTENT(inout), DIMENSION(nt_dro,nm_dro) :: macthb,&
       loglhb,logthb,logghb,mchb,yhb,zhb
  REAL(SP), INTENT(inout), DIMENSION(nt_dro) :: timehb
  REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
  REAL(SP), DIMENSION(nm) :: tphase=0.0
  INTEGER, INTENT(inout), DIMENSION(nt) :: nmass
  REAL(SP), INTENT(inout) :: hb_wght
  INTEGER, INTENT(in) :: t
  REAL(SP), INTENT(in) :: f_bhb, hbtime

  !number of blue HB to add per HB star
  !(not important b/c their total weight remains fixed)
  INTEGER, PARAMETER :: nhb=10
  INTEGER :: j, i, k, jrel=0, flip=0, tnhb, w_ind, nhbj
  REAL(SP) :: tgrad=0.,hblum=-999.,minteff=1E6,mint,maxt,mzams,lzams,azams,tnow
  REAL(SP) :: hb_wght0, hb_wght1
  REAL(SP), DIMENSION(nhb) :: dumarr=0.
  LOGICAL :: do_mod

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  hblum    = -999.
  flip     = 0
  hb_wght  = 0.
  hb_wght0 = 0.
  hb_wght1 = 0.
  tphase   = phase(t,:)

  !we need to count the total number of HB stars in 
  !these isochrones.  Also the minimum Teff for the HB
  IF (isoc_type.EQ.'bsti'.OR.isoc_type.EQ.'mist') THEN
     tnhb = 0
     minteff=1E6
     DO i=1,nm
        IF (tphase(i).GE.3.AND.tphase(i).LE.6) THEN
           tnhb=tnhb+1
           hb_wght0=hb_wght0+wght(i)
           !need this delta(mass) cut to remove the stars 
           !descending from the TRGB in the MIST models
           IF (logt(t,i).LT.minteff.AND.(mini(t,i+1)-mini(t,i)).GT.1E-6)  &
                minteff=logt(t,i)
        ENDIF
     ENDDO
     i=1
  ENDIF

  mint = 10**(hbtime-0.025) / 1E6
  maxt = 10**(hbtime+0.025) / 1E6

  WRITE(*,*) '----------------------'
  WRITE(*,*) 't        = ', t
  WRITE(*,*) 'time     = ', 10**hbtime
  WRITE(*,*) 'mint    = ', mint
  WRITE(*,*) 'maxt    = ', maxt
  WRITE(*,*) 'hb_wght0 = ', hb_wght0

  IF (isoc_type.EQ.'pdva') THEN

      DO j=2,nm

        IF (flip.NE.-1) THEN 
        
           tgrad = (logl(t,j)-logl(t,j-1)) / (mini(t,j)-mini(t,j-1))
           
           !if the lum jump is negative and large, we're on the HB
           IF (tgrad.LE.-5E2.AND.logl(t,j-1).GT.2.5) THEN
              flip = 1
              hblum = logl(t,j)
           ENDIF
        
           IF (flip.EQ.1) THEN

              !modify the HB
              !once lum increases by 0.1 dex we've left the HB
              IF (ABS(logl(t,j)-hblum).LT.0.1) THEN 

                 !keep track of total HB weight
                 hb_wght  = hb_wght+wght(j)
                 
                 !Blue HB stars have to be old
                 IF (f_bhb.GT.1E-3.AND.hbtime.GE.bhb_sbs_time) THEN

                    !add blue HB stars (their mass and Lbol remain the same)
                    mini(t,nmass(t)+1:nmass(t)+nhb)  = dumarr+mini(t,j)
                    mact(t,nmass(t)+1:nmass(t)+nhb)  = dumarr+mact(t,j)
                    logl(t,nmass(t)+1:nmass(t)+nhb)  = dumarr+logl(t,j)
                    phase(t,nmass(t)+1:nmass(t)+nhb) = dumarr+8.
                    DO i=1,nhb
                       !distribute Teff uniformly to high T
                       logt(t,nmass(t)+i) = logt(t,j)+(4.2-logt(t,j))*i/REAL(nhb)
                       !logt(t,nmass(t)+i) = 4.0+(4.3-4.0)*i/REAL(nhb)
                       !compute logg
                       logg(t,nmass(t)+i) = LOG10( gsig4pi*mact(t,nmass(t)+i)/&
                            10**logl(t,nmass(t)+i) ) + 4*logt(t,nmass(t)+i) 

                    ENDDO
                    wght(nmass(t)+1:nmass(t)+nhb)   = dumarr + &
                         f_bhb*wght(j)/nhb
                    !modify the weight of the existing HB stars
                    wght(j)   = wght(j) * (1-f_bhb)                 
                    !update number of stars in the isochrone
                    nmass(t) = nmass(t)+nhb
                 ENDIF

              ELSE 
                 flip = -1
              ENDIF
           
           ENDIF
        
        ENDIF
      
      ENDDO

   ELSE IF (isoc_type.EQ.'bsti'.OR.isoc_type.EQ.'mist') THEN

      !Blue HB stars have to be old
      !here, we're adding one additional BHB stars per CHeB star
      !we're also putting the original BHB star to the RC, so that
      !if f_bhb is small, then the actual BHB contribution is small
      !if f_bhb<1E-4, then the default MIST BHB is used
      do_mod = f_bhb.GT.1E-4.AND.hbtime.GE.bhb_sbs_time

      IF (do_mod) THEN

         !The stars in the DRO3 isochrones aren't all going to have counterparts in the MIST isochrones that they can
         !be "linked" to, because they span different parameter spaces (by design).  The DRO3 isochrones ONLY cover stars with
         !envelope masses smaller than that required to reach the thermally pulsing AGB stage.  And since the helium core
         !mass is essentially a constant, this restricts the resultant total stellar masses (core+envelope) to < ~1 Msun.
         !In addition, the DRO3 isochrones are sorted by age since the zero age horizontal branch (ZAHB), rather than 
         !absolute age.  This is a problem for us because we need to know the absolute age to know where to insert them
         !into the existing MIST isochrones.

         !This means we have to use a (non-ideal) approach of analytically estimating the main-sequence lifetime of each 
         !evolutionary track using the information we do have.  I take the simple approach of using a scaling relation 
         !where age = 10 Gyr * (M/Msun)**-2.5

         !We also need to know how to weight the different evolutionary tracks. What we would like to do is inform our choices 
         !on the weights based on the initial mass function, just like FSPS does with the basic isochrones. However, this is not 
         !a simple procedure because the DRO3 isochrones contain no information about the initial masses of the stars.  We have to 
         !make some assumptions about the relationship between the ZAHB mass and the ZAMS mass.  According to the DRO3 paper, stars 
         !lose between 0.1-0.3 Msun over the full duration of the RGB phase.  Then, if we just make a simple assumption of an average 
         !mass loss of 0.2 solar masses, we can estimate the initial mass.

         !We also need to be a bit careful because we can no longer have a direct relationship between the number of stars in the
         !default isochrones and those in the DRO3 isochrones.

         !loop through all of the indices on the new DRO3 isochrones
         nhbj=0
         DO jrel=1,nm_dro
            DO k=1,nt_dro
               !do not add a new HB star if the isochrones have no data for this j index
               IF (macthb(k,jrel).EQ.0.0) CYCLE
               
               mzams = macthb(k,jrel) + 0.2     ! estimate for the ZAMS mass (see the assumptions described above)
               lzams = mzams**4                 ! main sequence mass-luminosity relation for stars with 0.43 < M/Msun < 2
               !estimate the main sequence lifetime in Myr
               azams = 0.1*(1-yhb(k,jrel)-zhb(k,jrel))*0.00685*(mzams*msun)*(clight/1E8)**2 / (lzams*lsun)
               azams = azams / (3600.*24.*365.25*1e6)
               azams = MAX(MIN(azams, 15000.), 3200.)
               tnow = azams + timehb(k)

               !do not add a new HB star if time is not within this age bin
               IF (tnow.LT.mint.OR.tnow.GT.maxt) CYCLE  
               !update number of stars in the isochrone
               nmass(t) = nmass(t)+1
               !add blue HB stars (their mass and Lbol are now linked to the DRO3 isochrones)
               mini(t,nmass(t))  = mzams               ! estimate the initial mass using the prescription described above
               mact(t,nmass(t))  = macthb(k,jrel)      ! take actual mass directly from the DRO3 isochrones
               logl(t,nmass(t))  = loglhb(k,jrel)      ! do the same for logL
               logt(t,nmass(t))  = logthb(k,jrel)      ! and for logTeff
               logg(t,nmass(t))  = logghb(k,jrel)      ! and for logg
               phase(t,nmass(t)) = 8                ! FSPS default
               !logt(t,j) = minteff                 ! this moves the original HB star to the red cloud; probably should be commented out so we dont affect the original HB stars in the new prescription
               nhbj = nhbj + 1                      ! keep track of how many new stars have been added

               !find the weight of the closest existing mass bin (as the bins are meant to represent the entire range between -deltaM to +deltaM)
               !this will only be used to consider the *relative* weighting of these stars with respect to each other; the *absolute* scale will
               !be determined later
               w_ind = MINLOC(ABS(mini(t,1:nmass(t)-nhbj)-mini(t,nmass(t))),1)
               wght(nmass(t)) = wght(w_ind)
               hb_wght1 = hb_wght1 + wght(nmass(t))
            ENDDO
         ENDDO

         WRITE(*,*) 'hb_wght1 = ', hb_wght1

      ENDIF

      !the weights now also need to be modified based on 2 factors:
      ! 1. the relative total weights of horizontal branch stars from the original isochrones (hb_wght0) and the DRO3 ones (hb_wght1)
      ! 2. the f_bhb parameter, which should decide the relative weighting of the two populations

      !loop through all of the indices on the MIST isochrones
      DO j=2,nm
         IF (phase(t,j).GE.3.AND.phase(t,j).LE.6) THEN 
            !modify the weights of the original HB/AGB/PAGB stars
            IF (do_mod) wght(j) = (1-f_bhb)*wght(j)
            hb_wght = hb_wght + wght(j)
         ELSE IF (phase(t,j).EQ.8) THEN
            IF (do_mod) wght(j) = f_bhb*hb_wght0/hb_wght1*wght(j)
            hb_wght = hb_wght + wght(j)
         ENDIF
      ENDDO

   ENDIF

  WRITE(*,*) '----------------------'

  IF (nmass(t).GT.nm) THEN
     WRITE(*,*) 'MOD_HB ERROR: number of mass points GT nm', nmass(t), nm
     STOP
  ENDIF

  RETURN

END SUBROUTINE MOD_HB
