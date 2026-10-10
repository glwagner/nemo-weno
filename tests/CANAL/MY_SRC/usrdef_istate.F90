MODULE usrdef_istate
   !!======================================================================
   !!                     ***  MODULE usrdef_istate   ***
   !!
   !!                      ===  CANAL configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2017-11  (J. Chanut) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce        
   USE phycst         ! physical constants
   USE eosbn2, ONLY: rn_a0
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran    ! to use sign with key_nosignedzero
   !   
   USE usrdef_nam
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate       ! called by istate.F90
   PUBLIC   usr_def_istate_ssh   ! called by sshwzv.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here CANAL configuration 
      !!
      !! ** Method  :   Set a gaussian anomaly of pressure and associated
      !!                geostrophic velocities
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      !
      INTEGER  :: ji, jj, jk, jl  ! dummy loop indices
      REAL(wp) :: zx, zy, zP0, zumax, zlambda, zr_lambda2, zn2, zf0, zH, zrho1, za, zf, zdzF
      REAL(wp) :: zpsurf, zdyPs, zdxPs
      REAL(wp) :: zdt, zdu, zdv
      REAL(wp) :: zjetx, zjety, zbeta
      REAL(wp), DIMENSION(:,:), ALLOCATABLE  ::   zrandom
      INTEGER                                ::   iseed    ! size of the random seed
      INTEGER , DIMENSION(:)  , ALLOCATABLE  ::   iseeds   ! random seed
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : CANAL configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   '
      !
      zjetx = ABS(rn_ujetszx)/2._wp
      zjety = ABS(rn_ujetszy)/2._wp
      !
      zf0   = 2._wp * omega * SIN( rad * rn_ppgphi0 )
      !
      SELECT CASE(nn_initcase)

      CASE(-1)    ! stratif at rest

         ! temperature:
         pts(:,:,1,jp_tem) = 25. !!30._wp
         pts(:,:,2:jpk,jp_tem) = 22. !!24._wp
         ! salinity:  
         pts(:,:,:,jp_sal) = 35._wp
         ! velocities:
         pu(:,:,:) = 0._wp
         pv(:,:,:) = 0._wp

      CASE(0)    ! rest
         !
         ! temperature:
         pts(:,:,:,jp_tem) = 10._wp
         ! salinity:  
         pts(:,:,:,jp_sal) = 35._wp
         ! velocities:
         pu(:,:,:) = 0._wp
         pv(:,:,:) = 0._wp
         
      CASE(1)    ! geostrophic zonal jet from -zjety to +zjety
         !
         ! temperature:
         pts(:,:,:,jp_tem) = 10._wp
         ! salinity:  
         pts(:,:,jpk,jp_sal) = 0._wp
         DO jk=1, jpkm1
            WHERE( ABS(gphit) <= zjety )
!!$            WHERE( ABS(gphit) <= zjety*0.5 .AND. ABS(glamt) <= zjety*0.5 ) ! for a square of salt
               pts(:,:,jk,jp_sal) = 35._wp
            ELSEWHERE
               pts(:,:,jk,jp_sal) = 30._wp
            END WHERE                    
         END DO
         ! velocities:
         pu(:,:,:) = 0._wp
         DO jk=1, jpkm1
            WHERE( ABS(gphit) <= zjety ) pu(:,:,jk) = rn_uzonal
         END DO
         pv(:,:,:) = 0._wp
         !                  
      CASE(2)    ! geostrophic zonal current shear
         !
         ! temperature:
         pts(:,:,:,jp_tem) = 10._wp
         ! salinity:  
         pts(:,:,:,jp_sal) = 30._wp
         DO jk=1, jpkm1
            WHERE( ABS(gphiv) <= zjety ) pts(:,:,jk,jp_sal) = 30._wp + SIGN(1._wp,gphiv(:,:))
         END DO
         ! velocities:
         pu(:,:,:) = 0._wp
         DO jk=1, jpkm1
            WHERE( ABS(gphiv) <= zjety ) pu(:,:,jk) = SIGN(rn_uzonal,gphit(:,:))*SIGN(1._wp,rn_uzonal)
            WHERE( ABS(gphiv) == 0._wp ) pu(:,:,jk) = 0._wp
         END DO
         pv(:,:,:) = 0._wp
          !                  
      CASE(3)    ! gaussian zonal currant
         !
         ! zonal current
         pu(:,:,jpk) = 0._wp
         DO jk=1, jpkm1
            ! gphit and lambda are both in km
            pu(:,:,jk) = rn_uzonal * EXP( - 0.5_wp * gphit(:,:)**2 / rn_lambda**2 )
         END DO
         ! temperature:
         pts(:,:,:,jp_tem) = 10._wp
         ! salinity:  
         pts(:,:,jpk,jp_sal) = 0._wp
         DO jk=1, jpkm1
            pts(:,:,jk,jp_sal) = gphit(:,:)
         END DO
         ! velocities:
         pv(:,:,:) = 0._wp
         !            
      CASE(4)    ! geostrophic zonal pulse
         !
         DO_2D( 1, 1, 1, 1 )
            IF ( ABS(glamt(ji,jj)) <= zjetx ) THEN
               zdu = rn_uzonal
            ELSEIF ( ABS(glamt(ji,jj)) <= zjetx + 100._wp ) THEN
               zdu = rn_uzonal * ( ( zjetx-ABS(glamt(ji,jj)) )/100._wp + 1._wp )
            ELSE
               zdu = 0._wp
            ENDIF
            IF ( ABS(gphit(ji,jj)) <= zjety ) THEN
               pu(ji,jj,:) = zdu
               pts(ji,jj,:,jp_sal) = zdu / rn_uzonal + 1._wp
            ELSE
               pu(ji,jj,:) = 0._wp
               pts(ji,jj,:,jp_sal) = 1._wp
            ENDIF
         END_2D
         !
         ! temperature:
         pts(:,:,:,jp_tem) = 10._wp * ptmask(:,:,:)        
         pv(:,:,:) = 0._wp
         !
      CASE(5)    ! vortex
         !
         zf0   = 2._wp * omega * SIN( rad * rn_ppgphi0 )
         zumax = rn_vtxmax * SIGN(1._wp, zf0)  ! Here Anticyclonic: set zumax=-1 for cyclonic
         zlambda = SQRT(2._wp)*rn_lambda*1.e3       ! Horizontal scale in meters 
         zn2 = 3.e-3**2
         zH = 0.5_wp * 5000._wp
         !
         zr_lambda2 = 1._wp / zlambda**2
         zP0 = rho0 * zf0 * zumax * zlambda * SQRT(EXP(1._wp)/2._wp)
         !
         DO_2D( 1, 1, 1, 1 )
            zx = glamt(ji,jj) * 1.e3
            zy = gphit(ji,jj) * 1.e3
            ! Surface pressure: P(x,y,z) = F(z) * Psurf(x,y)
            zpsurf = zP0 * EXP(-(zx**2+zy**2)*zr_lambda2) - rho0 * ff_t(ji,jj) * rn_uzonal * zy
            ! temperature:
            DO jk=1,jpk
               zdt =  pdept(ji,jj,jk) 
               zrho1 = rho0 * (1._wp + zn2*zdt/grav)
               IF (zdt < zH) THEN
                  zdzF = (1._wp-EXP(zdt-zH)) / (zH-1._wp + EXP(-zH))   ! F'(z)
                  zrho1 = zrho1 - zdzF * zpsurf / grav    ! -1/g Dz(P) = -1/g * F'(z) * Psurf(x,y)
               ENDIF
               !               pts(ji,jj,jk,jp_tem) = (20._wp + (rho0-zrho1) / rn_a0) * ptmask(ji,jj,jk)
               pts(ji,jj,jk,jp_tem) = (10._wp + (rho0-zrho1) / rn_a0) * ptmask(ji,jj,jk)
            END DO
         END_2D
         !
         ! salinity:  
         pts(:,:,:,jp_sal) = 35._wp * ptmask(:,:,:) 
         !
         ! velocities:
         za = 2._wp * zP0 / zlambda**2
         DO_2D( 0, 0, 0, 0 )
            zx = glamu(ji,jj) * 1.e3
            zy = gphiu(ji,jj) * 1.e3
            DO jk=1, jpk
               zdu = 0.5_wp * (pdept(ji,jj,jk) + pdept(ji+1,jj,jk))
               IF (zdu < zH) THEN
                  zf = (zH-1._wp-zdu+EXP(zdu-zH)) / (zH-1._wp+EXP(-zH))
                  zdyPs = - za * zy * EXP(-(zx**2+zy**2)*zr_lambda2) - rho0 * ff_t(ji,jj) * rn_uzonal
                  pu(ji,jj,jk) = - zf / ( rho0 * ff_t(ji,jj) ) * zdyPs * ptmask(ji,jj,jk) * ptmask(ji+1,jj,jk)
               ELSE
                  pu(ji,jj,jk) = 0._wp
               ENDIF
            END DO
         END_2D
         !
         DO_2D( 0, 0, 0, 0 )
            zx = glamv(ji,jj) * 1.e3
            zy = gphiv(ji,jj) * 1.e3
            DO jk=1, jpk
               zdv = 0.5_wp * (pdept(ji,jj,jk) + pdept(ji,jj+1,jk))
               IF (zdv < zH) THEN
                  zf = (zH-1._wp-zdv+EXP(zdv-zH)) / (zH-1._wp+EXP(-zH))
                  zdxPs = - za * zx * EXP(-(zx**2+zy**2)*zr_lambda2)
                  pv(ji,jj,jk) = zf / ( rho0 * ff_f(ji,jj) ) * zdxPs * ptmask(ji,jj,jk) * ptmask(ji,jj+1,jk)
               ELSE
                  pv(ji,jj,jk) = 0._wp
               ENDIF
            END DO
         END_2D
         !            
      END SELECT
      !
      !                          !==  add noise  ==!
      IF (ln_sshnoise) THEN
         ALLOCATE(zrandom(Ni0glo,Nj0glo))
         ! fixed seed, identical on all processes: the noise is reproducible and independent of the domain decomposition
         ! (the default seed is not fixed, e.g. gfortran >= 7 seeds from the OS). Comment out to get a new seed at each run.
         CALL RANDOM_SEED( size = iseed )
         ALLOCATE( iseeds(iseed) )
         iseeds(:) = 1234567
         CALL RANDOM_SEED( put = iseeds )
         DEALLOCATE( iseeds )
         CALL RANDOM_NUMBER(zrandom)
         DO_2D( 0, 0, 0, 0 )
            pv(ji,jj,1) = pv(ji,jj,1) + 0.5_wp * ( zrandom(mig(ji,0),mjg(jj,0)) - 0.5_wp )
         END_2D
         DEALLOCATE(zrandom)
      ENDIF
      !
      CALL lbc_lnk( 'usrdef_istate', pts , 'T',  1._wp )
      CALL lbc_lnk( 'usrdef_istate', pu, 'U', -1._wp, pv, 'V', -1._wp )

   END SUBROUTINE usr_def_istate

  
   SUBROUTINE usr_def_istate_ssh( ptmask, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate_ssh  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here CANAL configuration 
      !!
      !! ** Method  :   Set ssh 
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      !
      INTEGER  :: ji, jj, jk, jl  ! dummy loop indices
      REAL(wp) :: zx, zy, zP0, zumax, zlambda, zr_lambda2, zn2, zf0, zH, zrho1, za, zf, zdzF
      REAL(wp) :: zpsurf, zdyPs, zdxPs
      REAL(wp) :: zdt, zdu, zdv
      REAL(wp) :: zjetx, zjety, zbeta
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate_ssh : CANAL configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   '
      !
      zjetx = ABS(rn_ujetszx)/2._wp
      zjety = ABS(rn_ujetszy)/2._wp
      !
      SELECT CASE(nn_initcase)
      CASE(0)                      !==   rest  ==!
         !
         pssh(:,:) = 0._wp
         !
      CASE(1)                      !==  geostrophic zonal jet from -zjety to +zjety  ==!
         !
         SELECT CASE( nn_fcase )
         CASE(0)                          !* f = f0 : ssh = - fuy / g
            WHERE( ABS(gphit) <= zjety )
               pssh(:,:) = - ff_t(:,:) * rn_uzonal * gphit(:,:) * 1.e3_wp / grav
            ELSEWHERE
               pssh(:,:) = - ff_t(:,:) * rn_uzonal * SIGN(zjety, gphit(:,:)) * 1.e3_wp / grav
            END WHERE
         CASE(1)                          !* f = f0 + beta*y : ssh = - u / g * ( fy + 0.5 * beta * y^2 )
            zbeta = 2._wp * omega * COS( rad * rn_ppgphi0 ) / ra
            WHERE( ABS(gphit) <= zjety )
               pssh(:,:) = - rn_uzonal / grav * ( ff_t(:,:) * gphit(:,:) * 1.e3_wp   &
                  &                             + 0.5_wp * zbeta * gphit(:,:) * gphit(:,:) * 1.e6_wp )
            ELSEWHERE
               pssh(:,:) = - rn_uzonal / grav * ( ff_t(:,:) * SIGN(zjety, gphit(:,:)) * 1.e3_wp   &
                  &                             + 0.5_wp * zbeta * zjety * zjety * 1.e6_wp )
            END WHERE
         END SELECT
         !                  
      CASE(2)                      !==  geostrophic zonal current shear  ==!
         !
         SELECT CASE( nn_fcase )
         CASE(0)                          !* f = f0 : ssh = - fuy / g
            WHERE( ABS(gphit) <= zjety )
               pssh(:,:) = - ff_t(:,:) * rn_uzonal * ABS(gphit(:,:)) * 1.e3_wp / grav
            ELSEWHERE
               pssh(:,:) = - ff_t(:,:) * rn_uzonal * zjety * 1.e3_wp / grav
            END WHERE
         CASE(1)                          !* f = f0 + beta*y : ssh = - u / g * ( fy + 0.5 * beta * y^2 )
            zbeta = 2._wp * omega * COS( rad * rn_ppgphi0 ) / ra
            WHERE( ABS(gphit) <= zjety )
               pssh(:,:) = - SIGN(rn_uzonal, gphit(:,:)) / grav   &
                  &        * ( ff_t(:,:) * gphit(:,:) * 1.e3_wp   &
                  &          + 0.5_wp * zbeta * gphit(:,:) * gphit(:,:) * 1.e6_wp )
            ELSEWHERE
               pssh(:,:) = - SIGN(rn_uzonal, gphit(:,:)) / grav   &
                  &        * ( ff_t(:,:) * SIGN(zjety, gphit(:,:)) * 1.e3_wp   &
                  &          + 0.5_wp * zbeta * zjety * zjety * 1.e6_wp )
            END WHERE
         END SELECT
         !                  
      CASE(3)                      !==  gaussian zonal currant  ==!
         !
         pssh(:,1) = - ff_t(:,1) / grav * e2t(:,1) * rn_uzonal * EXP( - 0.5_wp * gphit(:,1)**2 / rn_lambda**2 )
         DO jl=1, jpnj
            DO_2D( 0, 0, 0, 0 )
               pssh(ji,jj) = pssh(ji,jj-1) - ff_t(ji,jj) / grav * rn_uzonal * EXP( - 0.5_wp * gphit(ji,jj)**2 / rn_lambda**2 ) * e2t(ji,jj)
            END_2D
            CALL lbc_lnk( 'usrdef_istate_ssh', pssh, 'T',  1._wp )
         END DO
         !            
      CASE(4)                      !==  geostrophic zonal pulse !!st need to implement a way to separate ssh properly  ==!
         !
         DO_2D( 1, 1, 1, 1 )
            IF ( ABS(glamt(ji,jj)) <= zjetx ) THEN
               zdu = rn_uzonal
            ELSEIF ( ABS(glamt(ji,jj)) <= zjetx + 100._wp ) THEN
               zdu = rn_uzonal * ( ( zjetx-ABS(glamt(ji,jj)) )/100._wp + 1._wp )
            ELSE
               zdu = 0._wp
            ENDIF
            IF ( ABS(gphit(ji,jj)) <= zjety ) THEN
               pssh(ji,jj) = - ff_t(ji,jj) * zdu * gphit(ji,jj) * 1.e3_wp / grav
            ELSE
               pssh(ji,jj) = - ff_t(ji,jj) * zdu * SIGN(zjety,gphit(ji,jj)) * 1.e3_wp / grav 
            ENDIF
         END_2D
         !
      CASE(5)                    !==  vortex  ==!
         !
         zf0   = 2._wp * omega * SIN( rad * rn_ppgphi0 )
         zumax = rn_vtxmax * SIGN(1._wp, zf0)   ! Here Anticyclonic: set zumax=-1 for cyclonic
         zlambda = SQRT(2._wp)*rn_lambda        ! Horizontal scale in meters 
         zn2 = 3.e-3**2
         zH = 0.5_wp * 5000._wp
         !
         zr_lambda2 = 1._wp / zlambda**2
         zP0 = rho0 * zf0 * zumax * zlambda * SQRT(EXP(1._wp)/2._wp)
         !
         DO_2D( 1, 1, 1, 1 )
            zx = glamt(ji,jj) * 1.e3
            zy = gphit(ji,jj) * 1.e3
            !                                   ! Surface pressure: P(x,y,z) = F(z) * Psurf(x,y)
            zpsurf = zP0 * EXP(-(zx**2+zy**2)*zr_lambda2) - rho0 * ff_t(ji,jj) * rn_uzonal * zy
            pssh(ji,jj) = 0.
            DO jl=1,5
               zdt = pssh(ji,jj)
               zdzF = (1._wp - EXP(zdt-zH)) / (zH - 1._wp + EXP(-zH))   ! F'(z)
               zrho1 = rho0 * (1._wp + zn2*zdt/grav) - zdzF * zpsurf / grav    ! -1/g Dz(P) = -1/g * F'(z) * Psurf(x,y)
               pssh(ji,jj) = zpsurf / (zrho1*grav) * ptmask(ji,jj,1)   ! ssh = Psurf / (Rho*g)
            END DO
         END_2D
         !            
      END SELECT
      !
      CALL lbc_lnk( 'usrdef_istate_ssh', pssh, 'T',  1._wp )
      !
   END SUBROUTINE usr_def_istate_ssh
   
   !!======================================================================
END MODULE usrdef_istate
