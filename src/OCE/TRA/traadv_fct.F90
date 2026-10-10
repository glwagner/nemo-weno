MODULE traadv_fct
   !!==============================================================================
   !!                       ***  MODULE  traadv_fct  ***
   !! Ocean  tracers:  horizontal & vertical advective trend (2nd/4th order Flux Corrected Transport method)
   !!==============================================================================
   !! History :  3.7  !  2015-09  (L. Debreu, G. Madec)   original code (inspired from traadv_tvd.F90)
   !!            5.x  !  2026-03  (S. Griffies, G. Madec) thickness weighted tracer tendency 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  tra_adv_fct    : update the tracer trend with a 3D advective trends using a 2nd or 4th order FCT scheme
   !!                   with sub-time-stepping in the vertical direction
   !!  nonosc         : compute monotonic tracer fluxes by a non-oscillatory algorithm
   !!  interp_4th_cpt : 4th order compact scheme for the vertical component of the advection
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   USE trc_oce        ! share passive tracers/Ocean variables
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! tracers trends
   USE diaptr         ! poleward transport diagnostics
   USE diaar5         ! AR5 diagnostics
   USE phycst  , ONLY : rho0_rcp
   USE zdf_oce , ONLY : ln_zad_Aimp
   !
   USE in_out_manager ! I/O manager
   USE iom            !
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary condition (or mpp link)
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_fct        ! called by traadv.F90
   PUBLIC   interp_4th_cpt     ! called by traadv_cen.F90

   LOGICAL  ::   ll_zAimp1 ! flag to apply optimized adaptive implicit vertical advection
   LOGICAL  ::   ll_zAimp2 ! flag to apply accurate  adaptive implicit vertical advection
   LOGICAL  ::   l_trd     ! flag to compute trends
   LOGICAL  ::   l_ptr     ! flag to compute poleward transport
   LOGICAL  ::   l_hst     ! flag to compute heat/salt transport
   REAL(wp) ::   r1_6 = 1._wp / 6._wp   ! =1/6

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztmp, zwdia, zwinf, zwsup
   !                                        ! tridiag solver associated indices:
   !! INTEGER, PARAMETER ::   np_NH   = 0   ! Neumann homogeneous boundary condition
   !! INTEGER, PARAMETER ::   np_CEN2 = 1   ! 2nd order centered  boundary condition

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_fct( kt, kit000, cdtype, pdt, pU, pV, pW,       &
      &                    Kbb, Kmm, Kaa, pt, kjpt, Krhs, kn_fct_h, kn_fct_v, kn_fct_imp )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_fct  ***
      !!
      !! **  Purpose :   Compute the now trend due to total advection of tracers
      !!               and add it to the general trend of tracer equations
      !!
      !! **  Method  : - 2nd or 4th FCT scheme on the horizontal direction
      !!               (choice through the value of kn_fct)
      !!               - on the vertical the 4th order is a compact scheme
      !!               - corrected flux (monotonic correction)
      !!
      !! ** Action : - update pt(:,:,:,:,Krhs)  with the now advective tracer trends
      !!             - send trends to trdtra module for further diagnostics (l_trdtra=T)
      !!             - poleward advective heat and salt transport (ln_diaptr=T)
      !!----------------------------------------------------------------------
      INTEGER                                  , INTENT(in   ) ::   kt                   ! ocean time-step index
      INTEGER                                  , INTENT(in   ) ::   Kbb, Kmm, Kaa, Krhs  ! ocean time level indices
      INTEGER                                  , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                         , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                                  , INTENT(in   ) ::   kjpt            ! number of tracers
      INTEGER                                  , INTENT(in   ) ::   kn_fct_h        ! order of the FCT scheme (=2 or 4)
      INTEGER                                  , INTENT(in   ) ::   kn_fct_v        ! order of the FCT scheme (=2 or 4)
      INTEGER                                  , INTENT(in   ) ::   kn_fct_imp      ! treatment of implicit optmized(1) or accurate(2)
      REAL(wp)                                 , INTENT(in   ) ::   pdt             ! tracer time-step
      REAL(wp), DIMENSION(T2D(nn_hls),jpk     ), INTENT(in   ) ::   pU, pV, pW      ! ocean volume flux components [m^3/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt,jpt), INTENT(inout) ::   pt              ! tracers and RHS of tracer equation
      !
      INTEGER  ::   ji, jj, jk, jn                 ! dummy loop indices
      INTEGER  ::   kbnd
      REAL(wp) ::   ztra                           ! local scalar
      REAL(wp) ::   zC2t_u, zC4t_u, zC2t_v, zC4t_v !   -      -
      !
      REAL(wp), DIMENSION(T2D(nn_hls),jpk)    ::   zta_up1, ztFu, ztFv, ztFw
#if ! defined key_PSYCLONE_2p5p0
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztw
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   ztu, ztv !!, zltu, zltv
#else
      REAL(wp), DIMENSION(T2D(1),jpk)         ::   ztw
      REAL(wp), DIMENSION(T2D(2))             ::   ztu, ztv !!, zltu, zltv
#endif
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdx, ztrdy, ztrdz, zptry
      !!----------------------------------------------------------------------
      !
      ! Tiling can't be used with this routine due to lbc_lnk calls. Tiling should be paused in tra_adv as a workaround.
      IF( l_istiled ) CALL ctl_stop( 'Cannot use FCT advection with tiling active- this should not have happened!' )
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == kit000 )  THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'tra_adv_fct : FCT advection scheme on ', cdtype
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         ENDIF
         !
         l_trd = .FALSE.            ! set local switches
         l_hst = .FALSE.
         l_ptr = .FALSE.
         ll_zAimp1 = .FALSE.
         ll_zAimp2 = .FALSE.         
         IF( ( cdtype == 'TRA' .AND. l_trdtra ) .OR. ( cdtype == 'TRC' .AND. l_trdtrc ) )                      l_trd = .TRUE.
         IF(  l_diaptr .AND. cdtype == 'TRA' .AND. ( iom_use( 'sophtadv' )  .OR. iom_use( 'sopstadv'  ) )  )   l_ptr = .TRUE.
         IF(  l_diaar5 .AND. cdtype == 'TRA' .AND. ( iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. &
            &                                        iom_use("uadv_salttr") .OR. iom_use("vadv_salttr") )  )   l_hst = .TRUE.
      ENDIF
      ! Initialization of the last level otherwise mpi communications at north fold crash
      zta_up1(:,:,jpk) = 0._wp
      ztFu   (:,:,jpk) = 0._wp
      ztFv   (:,:,jpk) = 0._wp
      ztFw   (:,:,jpk) = 0._wp
      !
      IF( l_trd .OR. l_hst )  THEN
         ALLOCATE( ztrdx(T2D(1),jpk), ztrdy(T2D(1),jpk), ztrdz(T2D(0),jpk) )
         ztrdx(:,:,jpk) = 0._wp   ;   ztrdy(:,:,jpk) = 0._wp
         ztrdz(:,:,1  ) = 0._wp   ;   ztrdz(:,:,jpk) = 0._wp
      ENDIF
      !
      IF( l_ptr ) THEN
         ALLOCATE( zptry(T2D(0),jpk) )
         zptry(:,:,jpk) = 0._wp
      ENDIF
      !
      ! If adaptive vertical advection, check if it is needed on this PE at this time
      IF( ln_zad_Aimp ) THEN
         IF( MAXVAL( ABS( wi(T2D(1),:) ) ) > 0._wp ) THEN
            IF( kn_fct_imp == 1 )   ll_zAimp1 = .TRUE.
            IF( kn_fct_imp == 2 )   ll_zAimp2 = .TRUE.
         ENDIF
      END IF
      !
      IF( ll_zAimp2 )   ALLOCATE( ztmp(T2D(nn_hls),jpk), zwdia(T2D(nn_hls),jpk), zwinf(T2D(nn_hls),jpk), zwsup(T2D(nn_hls),jpk) )
      !
      DO jn = 1, kjpt            !==  loop over the tracers  ==!
         !
         !
         !        !==  Upstream and high order fluxes  ==!
         !
         ! -- Upstream fluxes
         ! ------------------
         CALL fct_up1_2stp( Kbb, Kmm, Kaa, pdt, pt(:,:,:,jn,Kbb), pU, pV, pW, ztFu, ztFv, ztFw, zta_up1, pt(:,:,:,jn,Krhs) )
         ! output => ztFu(1,0,1,0), ztFv(1,0,1,0), zta_up1(0,0,0,0), ztFw(0,0,0,0) if Amip2 or ztFw(1,1,1,1) if .not.Aimp2
         !
         IF( l_trd .OR. l_hst )  THEN             ! trend diagnostics (contribution of upstream fluxes)
            DO_3D( 1, 0, 1, 0, 1, jpkm1 )
               ztrdx(ji,jj,jk) = ztFu(ji,jj,jk)
               ztrdy(ji,jj,jk) = ztFv(ji,jj,jk)
            END_3D
            DO_3D( 0, 0, 0, 0, 2, jpkm1 )
               ztrdz(ji,jj,jk) = ztFw(ji,jj,jk)
            END_3D
         END IF
         !                             ! "Poleward" heat and salt transports (contribution of upstream fluxes)
         IF( l_ptr ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               zptry(ji,jj,jk) = ztFv(ji,jj,jk)
            END_3D
         ENDIF
         !
         !
         ! -- Anti-diffusive fluxes : high order minus low order
         ! ------------------------
         SELECT CASE( kn_fct_h )    !* horizontal anti-diffusive fluxes
         !
         CASE(  2  )                   !- 2nd order centered
            DO_3D( 1, 0, 1, 0, 1, jpkm1 )
               ztFu(ji,jj,jk) = 0.5_wp * pU(ji,jj,jk) * ( pt(ji,jj,jk,jn,Kmm) + pt(ji+1,jj,jk,jn,Kmm) ) - ztFu(ji,jj,jk)
               ztFv(ji,jj,jk) = 0.5_wp * pV(ji,jj,jk) * ( pt(ji,jj,jk,jn,Kmm) + pt(ji,jj+1,jk,jn,Kmm) ) - ztFv(ji,jj,jk)
            END_3D
            !
!!$         CASE(  4  )                   !- 4th order centered
!!$            !
!!$            ALLOCATE( ztu(T2D(2)), ztv(T2D(2)), zltu(T2D(2)), zltv(T2D(2)) )
!!$            !
!!$            DO jk = 1, jpkm1                 ! Laplacian
!!$               DO_2D( 2, 1, 2, 1 )                 ! 1st derivative (gradient)
!!$                  ztu(ji,jj) = ( pt(ji+1,jj  ,jk,jn,Kmm) - pt(ji,jj,jk,jn,Kmm) ) * umask(ji,jj,jk)
!!$                  ztv(ji,jj) = ( pt(ji  ,jj+1,jk,jn,Kmm) - pt(ji,jj,jk,jn,Kmm) ) * vmask(ji,jj,jk)
!!$               END_2D
!!$               DO_2D( 1, 1, 1, 1 )                 ! 2nd derivative * 1/ 6
!!$                  zltu(ji,jj) = (  ztu(ji,jj) + ztu(ji-1,jj)  ) * r1_6
!!$                  zltv(ji,jj) = (  ztv(ji,jj) + ztv(ji,jj-1)  ) * r1_6
!!$               END_2D
!!$               ! NOTE [ comm_cleanup ] : need to change sign to ensure halo 1 - halo 2 compatibility
!!$               !
!!$               DO_2D( 1, 0, 1, 0 )
!!$                  zC2t_u = pt(ji,jj,jk,jn,Kmm) + pt(ji+1,jj  ,jk,jn,Kmm)   ! 2 x C2 interpolation of T at u- & v-points
!!$                  zC2t_v = pt(ji,jj,jk,jn,Kmm) + pt(ji  ,jj+1,jk,jn,Kmm)
!!$                  !                                                        ! C4 interpolation of T at u- & v-points (x2)
!!$                  zC4t_u = zC2t_u + ( zltu(ji,jj) - zltu(ji+1,jj  ) )      !    round brackets added to fix the order of floating point operations
!!$                  zC4t_v = zC2t_v + ( zltv(ji,jj) - zltv(ji  ,jj+1) )      !    needed to ensure the North Pole reproducibility
!!$                  !                                                        ! C4 minus upstream advective fluxes
!!$                  ztFu(ji,jj,jk) =  0.5_wp * pU(ji,jj,jk) * zC4t_u - ztFu(ji,jj,jk)
!!$                  ztFv(ji,jj,jk) =  0.5_wp * pV(ji,jj,jk) * zC4t_v - ztFv(ji,jj,jk)
!!$               END_2D
!!$            END DO
!!$            !
!!$            DEALLOCATE( ztu, ztv, zltu, zltv )  
            !
         CASE(  4  )                   !- 4th order centered
            !
#if ! defined key_PSYCLONE_2p5p0
            ALLOCATE( ztu(T2D(2)), ztv(T2D(2)) )  
#endif
            !
            DO jk = 1, jpkm1
               DO_2D( 2, 1, 2, 1 )    ! 1st derivative (gradient)
                  ztu(ji,jj) = ( pt(ji+1,jj  ,jk,jn,Kmm) - pt(ji,jj,jk,jn,Kmm) ) * umask(ji,jj,jk)
                  ztv(ji,jj) = ( pt(ji  ,jj+1,jk,jn,Kmm) - pt(ji,jj,jk,jn,Kmm) ) * vmask(ji,jj,jk)
               END_2D
               !
               DO_2D( 1, 0, 1, 0 )    ! Horizontal advective fluxes
                  zC2t_u = pt(ji,jj,jk,jn,Kmm) + pt(ji+1,jj  ,jk,jn,Kmm)   ! 2 x C2 interpolation of T at u- & v-points (x2)
                  zC2t_v = pt(ji,jj,jk,jn,Kmm) + pt(ji  ,jj+1,jk,jn,Kmm)
                  !                                                        ! C4 interpolation of T at u- & v-points (x2)
                  zC4t_u =  zC2t_u + r1_6 * ( ztu(ji-1,jj  ) - ztu(ji+1,jj  ) )
                  zC4t_v =  zC2t_v + r1_6 * ( ztv(ji  ,jj-1) - ztv(ji  ,jj+1) )
                  !                                                        ! C4 minus upstream advective fluxes
                  ztFu(ji,jj,jk) =  0.5_wp * pU(ji,jj,jk) * zC4t_u - ztFu(ji,jj,jk)
                  ztFv(ji,jj,jk) =  0.5_wp * pV(ji,jj,jk) * zC4t_v - ztFv(ji,jj,jk)
               END_2D
            END DO
            !
#if ! defined key_PSYCLONE_2p5p0
            DEALLOCATE( ztu, ztv )  
#endif
            !
         END SELECT
         !
         ! loop boundaries for ztFw (avoid a communication if .not.Aimp2)
         IF( ll_zAimp2 ) THEN ; kbnd = 0
         ELSE                 ; kbnd = 1 ; ENDIF
         !
         SELECT CASE( kn_fct_v )    !* vertical anti-diffusive fluxes (w-masked interior values)
         !
         CASE(  2  )                   !- 2nd order centered
            !
            DO_3D( kbnd, kbnd, kbnd, kbnd, 2, jpkm1 )
               ztFw(ji,jj,jk) =  (  pW(ji,jj,jk) * 0.5_wp * ( pt(ji,jj,jk,jn,Kmm) + pt(ji,jj,jk-1,jn,Kmm) )   &
                  &              - ztFw(ji,jj,jk)  ) * wmask(ji,jj,jk)
            END_3D
            !
         CASE(  4  )                   !- 4th order COMPACT
            !
#if ! defined key_PSYCLONE_2p5p0
            ALLOCATE( ztw(T2D(1),jpk) )
#endif
            !
            CALL interp_4th_cpt( pt(:,:,:,jn,Kmm) , ztw )   ! zwt = COMPACT interpolation of T at w-point
            DO_3D( kbnd, kbnd, kbnd, kbnd, 2, jpkm1 )
               ztFw(ji,jj,jk) = ( pW(ji,jj,jk) * ztw(ji,jj,jk) - ztFw(ji,jj,jk) ) * wmask(ji,jj,jk)
            END_3D
            !
#if ! defined key_PSYCLONE_2p5p0
            DEALLOCATE( ztw )
#endif
            !
         END SELECT

         IF( lk_linssh ) THEN    ! top ocean value: high order = upstream  ==>>  ztFw=0
            !                    !            note: for non linear ssh, ztFw at surface is alreay set to 0
            DO_2D( kbnd, kbnd, kbnd, kbnd )
               ztFw(ji,jj,1) = 0._wp   ! only ocean surface as interior ztFw values have been w-masked
            END_2D
         ENDIF
         !
         IF ( ll_zAimp2 ) THEN
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )    !* trend and after field with monotonic scheme
               !                                                ! total intermediate advective trends
               ztra = - (  ( ztFu(ji,jj,jk) - ztFu(ji-1,jj  ,jk  ) )   &   ! add () NP halo
                  &      + ( ztFv(ji,jj,jk) - ztFv(ji  ,jj-1,jk  ) )   &
                  &      + ( ztFw(ji,jj,jk) - ztFw(ji  ,jj  ,jk+1) ) ) * r1_e1e2t(ji,jj)
               ztmp(ji,jj,jk) = zta_up1(ji,jj,jk) + pdt * ztra / e3t(ji,jj,jk,Kaa) * tmask(ji,jj,jk)
            END_3D
            !
            CALL tridia_solver( zwdia, zwsup, zwinf, ztmp, ztmp, 0, kbnd=0 )
            !
            DO_3D( 0, 0, 0, 0, 2, jpkm1 )       ! Interior value ( multiplied by wmask)
               ztFw(ji,jj,jk) =  ztFw(ji,jj,jk) + ( MAX( wi(ji,jj,jk) , 0._wp ) * ztmp(ji,jj,jk) + &
                  &                                 MIN( wi(ji,jj,jk) , 0._wp ) * ztmp(ji,jj,jk-1) ) * e1e2t(ji,jj) * wmask(ji,jj,jk)
            END_3D
         END IF
         !
         !        !==  monotonicity algorithm  ==!
         !
         IF( ln_zad_Aimp .AND. kn_fct_imp == 2 ) THEN
            CALL lbc_lnk( 'traadv_fct', ztFu, 'U', -1.0_wp, ztFv, 'V', -1.0_wp, ztFw, 'T', 1.0_wp, zta_up1, 'T', 1.0_wp )
         ELSE
            CALL lbc_lnk( 'traadv_fct', ztFu, 'U', -1.0_wp, ztFv, 'V', -1.0_wp, zta_up1, 'T', 1.0_wp )
         ENDIF
         !
         ! -- Flux limiter
         ! ---------------
         CALL nonosc( Kaa, pt(:,:,:,jn,Kbb), ztFu, ztFv, ztFw, zta_up1, pdt )
         !
         !        !==  final trend with corrected fluxes  ==!
         !
         ! -- Final trend with corrected fluxes
         ! ------------------------------------
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            ztra = - (  ( ztFu(ji,jj,jk) - ztFu(ji-1,jj  ,jk  ) )   &   ! add () for NP reproducibility
               &      + ( ztFv(ji,jj,jk) - ztFv(ji  ,jj-1,jk  ) )   &
               &      + ( ztFw(ji,jj,jk) - ztFw(ji  ,jj  ,jk+1) ) ) * r1_e1e2t(ji,jj)
            !
            pt(ji,jj,jk,jn,Krhs) = pt(ji,jj,jk,jn,Krhs) + ztra * tmask(ji,jj,jk)                     !!clem
            ! update upstream (for implicit)
            zta_up1(ji,jj,jk) = zta_up1(ji,jj,jk) + pdt * ztra / e3t(ji,jj,jk,Kaa) * tmask(ji,jj,jk)
         END_3D
         !
         IF ( ll_zAimp2 ) THEN
            !
            ztmp(:,:,1) = 0._wp ; ztmp(:,:,jpk) = 0._wp
            DO_3D( 0, 0, 0, 0, 2, jpkm1 )      ! Interior value ( multiplied by wmask)
               ztmp(ji,jj,jk) =  - ( MAX( wi(ji,jj,jk) , 0._wp ) * zta_up1(ji,jj,jk) + &
                  &                  MIN( wi(ji,jj,jk) , 0._wp ) * zta_up1(ji,jj,jk-1) ) * wmask(ji,jj,jk)
               !
               ztFw(ji,jj,jk) = ztFw(ji,jj,jk) + ztmp(ji,jj,jk) * e1e2t(ji,jj)  ! update vertical fluxes
            END_3D
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               pt(ji,jj,jk,jn,Krhs) = pt(ji,jj,jk,jn,Krhs) - ( ztmp(ji,jj,jk) - ztmp(ji,jj,jk+1) ) 
            END_3D
         END IF
         !
         IF( l_trd .OR. l_hst ) THEN   ! trend diagnostics // heat/salt transport
            DO_3D( 1, 0, 1, 0, 1, jpkm1 )
               ztrdx(ji,jj,jk) = ztrdx(ji,jj,jk) + ztFu(ji,jj,jk)  ! <<< add anti-diffusive fluxes
               ztrdy(ji,jj,jk) = ztrdy(ji,jj,jk) + ztFv(ji,jj,jk)  !     to upstream fluxes
            END_3D                                                 !
            DO_3D( 0, 0, 0, 0, 2, jpkm1 )                          !
               ztrdz(ji,jj,jk) = ztrdz(ji,jj,jk) + ztFw(ji,jj,jk)  !
            END_3D                                                 !
            !
            IF( l_trd ) THEN              ! trend diagnostics 
               CALL trd_tra( kt, Kmm, Krhs, cdtype, jn, jptra_xad, ztrdx, pU, pt(:,:,:,jn,Kmm) )
               CALL trd_tra( kt, Kmm, Krhs, cdtype, jn, jptra_yad, ztrdy, pV, pt(:,:,:,jn,Kmm) )
               CALL trd_tra( kt, Kmm, Krhs, cdtype, jn, jptra_zad, ztrdz, pW, pt(:,:,:,jn,Kmm) )
            ENDIF
            !                             ! heat/salt transport
            IF( l_hst )   CALL dia_ar5_hst( jn, 'adv', ztrdx(:,:,:), ztrdy(:,:,:) )
            !
         ENDIF
         IF( l_ptr ) THEN              ! "Poleward" transports
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               zptry(ji,jj,jk) = zptry(ji,jj,jk) + ztFv(ji,jj,jk)  ! <<< add anti-diffusive fluxes
            END_3D
            CALL dia_ptr_hst( jn, 'adv', zptry(:,:,:) )
         ENDIF
         !
      END DO                     !==  end loop over the tracers  ==!
      !
      IF( ll_zAimp2 )          DEALLOCATE( ztmp, zwdia, zwinf, zwsup )
      IF( l_trd .OR. l_hst )   DEALLOCATE( ztrdx, ztrdy, ztrdz )
      IF( l_ptr )              DEALLOCATE( zptry )
      !
   END SUBROUTINE tra_adv_fct


   SUBROUTINE fct_up1_1stp( Kbb, Kmm, Kaa, pDt, pt_b, pU, pV, pW, ptFu, ptFv, ptFw, pt_up1, pt_rhs )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE upstream  ***
      !!
      !! **  Purpose :   compute the upstream fluxes and upstream guess of tracer
      !! 
!! gm / smg: this routine is not used in RK3 and so it should be removed. 
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::   Kbb, Kmm, Kaa    ! ocean time level indices
      REAL(wp)                            , INTENT(in   ) ::   pDt              ! tracer time-step
      REAL(wp), DIMENSION(:,:,:)          , INTENT(in   ) ::   pt_b             ! before tracer
      REAL(wp), DIMENSION(T2D(nn_hls),jpk), INTENT(in   ) ::   pU, pV, pW       ! volume transport components [m^3/s]
      REAL(wp), DIMENSION(T2D(nn_hls),jpk), INTENT(  out) ::   ptFu, ptFv, ptFw ! upstream fluxes
      REAL(wp), DIMENSION(T2D(nn_hls),jpk), INTENT(  out) ::   pt_up1           ! tracer or upstream guess of tracer
      REAL(wp), DIMENSION(:,:,:)          , INTENT(inout) ::   pt_rhs           ! RHS tendency
      !
      INTEGER  ::   ji, jj, jk    ! dummy loop indices
      INTEGER  ::   ik            ! local integer
      REAL(wp) ::   ztra   ! local scalar
      !!----------------------------------------------------------------------
!!$      IF( .NOT. ll_upsxy ) THEN         !** no alternate directions **!
!!$      ELSE                              !** alternate directions **!
!!$      ENDIF
      !
      ! -- upstream fluxes

      ! horizontal
      DO_3D( 1, 0, 1, 0, 1, jpkm1 )
         ptFu(ji,jj,jk) = ( MAX( pU(ji,jj,jk) , 0._wp ) * pt_b(ji,jj,jk) + MIN( pU(ji,jj,jk) , 0._wp ) * pt_b(ji+1,jj  ,jk) )
         ptFv(ji,jj,jk) = ( MAX( pV(ji,jj,jk) , 0._wp ) * pt_b(ji,jj,jk) + MIN( pV(ji,jj,jk) , 0._wp ) * pt_b(ji  ,jj+1,jk) )
      END_3D

      ! vertical
      DO_3D( 1, 1, 1, 1, 2, jpkm1 )
         ptFw(ji,jj,jk) = ( MAX( pW(ji,jj,jk) , 0._wp ) * pt_b(ji,jj,jk) + MIN( pW(ji,jj,jk) , 0._wp ) * pt_b(ji,jj,jk-1)   ) * wmask(ji,jj,jk) !!clem
         !                                                                                                                    ! should not be necessary to mask by wmask because W is already masked
      END_3D
      IF( lk_linssh ) THEN               ! top ocean value (only in linear free surface as ptFw has been w-masked)
         IF( ln_isfcav ) THEN                        ! top of the ice-shelf cavities and at the ocean surface
            DO_2D( 1, 1, 1, 1 )
               ik = mikt(ji,jj)
               ptFw(ji,jj,ik) = pW(ji,jj,ik) * pt_b(ji,jj,ik)
            END_2D
         ELSE                                        ! no cavities: only at the ocean surface
            DO_2D( 1, 1, 1, 1 )
               ptFw(ji,jj,1) = pW(ji,jj,1) * pt_b(ji,jj,1)
            END_2D
         ENDIF
      ELSE
         DO_2D( 1, 1, 1, 1 )
            ptFw(ji,jj,1) = 0._wp
         END_2D
      ENDIF

      ! -- after tracer with upstream scheme

      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         ztra = - (  ( ptFu(ji,jj,jk) - ptFu(ji-1,jj  ,jk  ) )   &
            &      + ( ptFv(ji,jj,jk) - ptFv(ji  ,jj-1,jk  ) )   &
            &      + ( ptFw(ji,jj,jk) - ptFw(ji  ,jj  ,jk+1) ) ) * r1_e1e2t(ji,jj)

         pt_rhs(ji,jj,jk) = pt_rhs(ji,jj,jk) + ztra / e3t(ji,jj,jk,Kmm) * tmask(ji,jj,jk) !!clem

         IF ( ll_zAimp1 ) THEN
            ztra = ztra - ( wi(ji,jj,jk) - wi(ji,jj,jk+1) ) * pt_b(ji,jj,jk)
         END IF

         pt_up1(ji,jj,jk) = ( e3t(ji,jj,jk,Kbb) * pt_b(ji,jj,jk) + pDt * ztra ) / e3t(ji,jj,jk,Kaa) * tmask(ji,jj,jk) !!clem
      END_3D

      IF ( ll_zAimp2 ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            zwdia(ji,jj,jk) = 1._wp + pDt * ( MAX( wi(ji,jj,jk  ) , 0._wp ) - MIN( wi(ji,jj,jk+1) , 0._wp ) ) / e3t(ji,jj,jk,Kaa)
            zwinf(ji,jj,jk) =         pDt *   MIN( wi(ji,jj,jk  ) , 0._wp )                                   / e3t(ji,jj,jk,Kaa)
            zwsup(ji,jj,jk) =       - pDt *   MAX( wi(ji,jj,jk+1) , 0._wp )                                   / e3t(ji,jj,jk,Kaa)
         END_3D
         !
         CALL tridia_solver( zwdia, zwsup, zwinf, pt_up1, pt_up1, 0, kbnd=0 )
         !
         ztmp(:,:,1) = 0._wp ; ztmp(:,:,jpk) = 0._wp
         DO_3D( 0, 0, 0, 0, 2, jpkm1 )       ! Interior value ( multiplied by wmask)
            ztmp(ji,jj,jk) =  ( MAX( wi(ji,jj,jk) , 0._wp ) * pt_up1(ji,jj,jk) + &
               &                MIN( wi(ji,jj,jk) , 0._wp ) * pt_up1(ji,jj,jk-1) ) * wmask(ji,jj,jk)
            !
            !!clem: not sure we have to update the vertical flux?
            ptFw(ji,jj,jk) = ptFw(ji,jj,jk) + ztmp(ji,jj,jk) * e1e2t(ji,jj)  ! update vertical fluxes
         END_3D
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            pt_rhs(ji,jj,jk) = pt_rhs(ji,jj,jk) - ( ztmp(ji,jj,jk) - ztmp(ji,jj,jk+1) ) / e3t(ji,jj,jk,Kmm)
         END_3D
      ENDIF
      !
   END SUBROUTINE fct_up1_1stp

   
   SUBROUTINE fct_up1_2stp( Kbb, Kmm, Kaa, pDt, pt_b, pU, pV, pW, ptFu, ptFv, ptFw, pt_up1, pt_rhs )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE upstream  ***
      !!
      !! **  Purpose :   compute the upstream fluxes and upstream guess of tracer
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::   Kbb, Kmm, Kaa    ! ocean time level indices
      REAL(wp)                            , INTENT(in   ) ::   pDt              ! tracer time-step
      REAL(wp), DIMENSION(:,:,:)          , INTENT(in   ) ::   pt_b             ! before tracer
      REAL(wp), DIMENSION(T2D(nn_hls),jpk), INTENT(in   ) ::   pU, pV, pW       ! ocean volume flux components [m^3/s]
      REAL(wp), DIMENSION(T2D(nn_hls),jpk), INTENT(  out) ::   ptFu, ptFv, ptFw ! upstream fluxes
      REAL(wp), DIMENSION(T2D(nn_hls),jpk), INTENT(  out) ::   pt_up1           ! tracer or upstream guess of tracer
      REAL(wp), DIMENSION(:,:,:)          , INTENT(inout) ::   pt_rhs           ! RHS tendency
      !
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      INTEGER  ::   ik              ! local integer
      REAL(wp) ::   zimp, zDt, ztra ! local scalar
      !!----------------------------------------------------------------------
!!$      IF( .NOT. ll_upsxy ) THEN         !** no alternate directions **!
!!$      ELSE                              !** alternate directions **!
!!$      ENDIF     
      !

      ! *** 1st step
      zDt = 0.5_wp * pDt
      
      ! -- upstream fluxes

      ! horizontal
      DO_3D( 2, 1, 2, 1, 1, jpkm1 )
         ptFu(ji,jj,jk) = MAX( pU(ji,jj,jk) , 0._wp ) * pt_b(ji,jj,jk) + MIN( pU(ji,jj,jk) , 0._wp ) * pt_b(ji+1,jj  ,jk)
         ptFv(ji,jj,jk) = MAX( pV(ji,jj,jk) , 0._wp ) * pt_b(ji,jj,jk) + MIN( pV(ji,jj,jk) , 0._wp ) * pt_b(ji  ,jj+1,jk)
      END_3D
      
      ! vertical
      DO_3D( 1, 1, 1, 1, 2, jpkm1 )
         ptFw(ji,jj,jk) = MAX( pW(ji,jj,jk) , 0._wp ) * pt_b(ji,jj,jk) + MIN( pW(ji,jj,jk) , 0._wp ) * pt_b(ji,jj,jk-1) * wmask(ji,jj,jk) !!clem
         !                                                                                                               ! should not be necessary to mask by wmask because W is already masked
      END_3D
      IF( lk_linssh ) THEN               ! top ocean value (only in linear free surface as ptFw has been w-masked)
         IF( ln_isfcav ) THEN                        ! top of the ice-shelf cavities and at the ocean surface
            DO_2D( 1, 1, 1, 1 )
               ik = mikt(ji,jj)
               ptFw(ji,jj,ik) = pW(ji,jj,ik) * pt_b(ji,jj,ik)
            END_2D
         ELSE                                        ! no cavities: only at the ocean surface
            DO_2D( 1, 1, 1, 1 )
               ptFw(ji,jj,1) = pW(ji,jj,1) * pt_b(ji,jj,1)
            END_2D
         ENDIF
      ELSE
         DO_2D( 1, 1, 1, 1 )
            ptFw(ji,jj,1) = 0._wp
         END_2D
      ENDIF

      ! -- middle time step tracer with upstream scheme
      
      DO_3D( 1, 1, 1, 1, 1, jpkm1 )
         ztra = - (  ( ptFu(ji,jj,jk) - ptFu(ji-1,jj  ,jk  ) )   &
            &      + ( ptFv(ji,jj,jk) - ptFv(ji  ,jj-1,jk  ) )   &
            &      + ( ptFw(ji,jj,jk) - ptFw(ji  ,jj  ,jk+1) ) ) * r1_e1e2t(ji,jj)

         IF( ll_zAimp1 )   ztra = ztra - ( wi(ji,jj,jk) - wi(ji,jj,jk+1) ) * pt_b(ji,jj,jk)
         !
         ! guess upstream at Kmm
         pt_up1(ji,jj,jk) = ( e3t(ji,jj,jk,Kbb) * pt_b(ji,jj,jk) + zDt * ztra ) / e3t(ji,jj,jk,Kmm) * tmask(ji,jj,jk) !!clem
      END_3D


      IF ( ll_zAimp2 ) THEN
         DO_3D( 1, 1, 1, 1, 1, jpkm1 )
            zwdia(ji,jj,jk) = 1._wp + zDt * ( MAX( wi(ji,jj,jk  ) , 0._wp ) - MIN( wi(ji,jj,jk+1) , 0._wp ) ) / e3t(ji,jj,jk,Kmm)
            zwinf(ji,jj,jk) =         zDt *   MIN( wi(ji,jj,jk  ) , 0._wp )                                   / e3t(ji,jj,jk,Kmm)
            zwsup(ji,jj,jk) =       - zDt *   MAX( wi(ji,jj,jk+1) , 0._wp )                                   / e3t(ji,jj,jk,Kmm)
         END_3D
         !
         CALL tridia_solver( zwdia, zwsup, zwinf, pt_up1, pt_up1, 0, kbnd=1 )
         !
         ztmp(:,:,1) = 0._wp ; ztmp(:,:,jpk) = 0._wp
         DO_3D( 1, 1, 1, 1, 2, jpkm1 )       ! Interior value ( multiplied by wmask)
            ztmp(ji,jj,jk) =  ( MAX( wi(ji,jj,jk) , 0._wp ) * pt_up1(ji,jj,jk) + &
               &                MIN( wi(ji,jj,jk) , 0._wp ) * pt_up1(ji,jj,jk-1) ) * wmask(ji,jj,jk)
            !
            !!clem: not sure we have to update the vertical flux. I would say no...
            ptFw(ji,jj,jk) = ptFw(ji,jj,jk) + ztmp(ji,jj,jk) * e1e2t(ji,jj)  ! update vertical fluxes
         END_3D
      ENDIF

      
      ! *** 2nd step
      zDt = pDt
      
      ! -- upstream fluxes

      ! horizontal
      DO_3D( 1, 0, 1, 0, 1, jpkm1 )
         ptFu(ji,jj,jk) = 0.5_wp * ( ptFu(ji,jj,jk) + MAX( pU(ji,jj,jk) , 0._wp ) * pt_up1(ji  ,jj,jk) &
            &                                       + MIN( pU(ji,jj,jk) , 0._wp ) * pt_up1(ji+1,jj,jk) )
         ptFv(ji,jj,jk) = 0.5_wp * ( ptFv(ji,jj,jk) + MAX( pV(ji,jj,jk) , 0._wp ) * pt_up1(ji,jj  ,jk) &
            &                                       + MIN( pV(ji,jj,jk) , 0._wp ) * pt_up1(ji,jj+1,jk) )
      END_3D
      
      ! vertical
      DO_3D( 1, 1, 1, 1, 2, jpkm1 )
         ptFw(ji,jj,jk) = 0.5_wp * ( ptFw(ji,jj,jk) + ( MAX( pW(ji,jj,jk) , 0._wp ) * pt_up1(ji,jj,jk  ) &
            &                                       +   MIN( pW(ji,jj,jk) , 0._wp ) * pt_up1(ji,jj,jk-1) ) * wmask(ji,jj,jk) ) ! wmask for cavities and linssh
         !                                                                                                                       
      END_3D
      IF( lk_linssh ) THEN               ! top ocean value (only in linear free surface as ptFw has been w-masked)
         IF( ln_isfcav ) THEN                        ! top of the ice-shelf cavities and at the ocean surface
            DO_2D( 1, 1, 1, 1 )
               ik = mikt(ji,jj)
               IF( ik == 1 ) THEN
                  ptFw(ji,jj,ik) = 0.5_wp * ( ptFw(ji,jj,ik) + pW(ji,jj,ik) * pt_up1(ji,jj,ik) )
               ELSE
                  ptFw(ji,jj,ik) = ( ptFw(ji,jj,ik) + 0.5_wp * pW(ji,jj,ik) * pt_up1(ji,jj,ik) )
               ENDIF
            END_2D
         ELSE                                        ! no cavities: only at the ocean surface
            DO_2D( 1, 1, 1, 1 )
               ptFw(ji,jj,1) = 0.5_wp * ( ptFw(ji,jj,1) + pW(ji,jj,1) * pt_up1(ji,jj,1) )
            END_2D
         ENDIF
      ENDIF

      ! -- after tracer with upstream scheme
      
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         ztra = - (  ( ptFu(ji,jj,jk) - ptFu(ji-1,jj  ,jk  ) )   &
            &      + ( ptFv(ji,jj,jk) - ptFv(ji  ,jj-1,jk  ) )   &
            &      + ( ptFw(ji,jj,jk) - ptFw(ji  ,jj  ,jk+1) ) ) * r1_e1e2t(ji,jj)

         pt_rhs(ji,jj,jk) = pt_rhs(ji,jj,jk) + ztra                     * tmask(ji,jj,jk) 
         
         IF( ll_zAimp1 )   ztra = ztra - ( wi(ji,jj,jk) - wi(ji,jj,jk+1) ) * pt_b(ji,jj,jk)
         !
         pt_up1(ji,jj,jk) = ( e3t(ji,jj,jk,Kbb) * pt_b(ji,jj,jk) + zDt * ztra ) / e3t(ji,jj,jk,Kaa) * tmask(ji,jj,jk) !!clem
      END_3D

      IF ( ll_zAimp2 ) THEN
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            zwdia(ji,jj,jk) = 1._wp + 0.5_wp * zDt * ( MAX( wi(ji,jj,jk  ) , 0._wp ) - MIN( wi(ji,jj,jk+1) , 0._wp ) ) / e3t(ji,jj,jk,Kaa)
            zwinf(ji,jj,jk) =         0.5_wp * zDt *   MIN( wi(ji,jj,jk  ) , 0._wp )                                   / e3t(ji,jj,jk,Kaa)
            zwsup(ji,jj,jk) =       - 0.5_wp * zDt *   MAX( wi(ji,jj,jk+1) , 0._wp )                                   / e3t(ji,jj,jk,Kaa)
         END_3D
         !
         CALL tridia_solver( zwdia, zwsup, zwinf, pt_up1, pt_up1, 0, kbnd=0 )
         !
         ztmp(:,:,1) = 0._wp ; ztmp(:,:,jpk) = 0._wp
         DO_3D( 0, 0, 0, 0, 2, jpkm1 )       ! Interior value ( multiplied by wmask)
            ztmp(ji,jj,jk) =  0.5_wp * ( MAX( wi(ji,jj,jk) , 0._wp ) * pt_up1(ji,jj,jk) + &
               &                         MIN( wi(ji,jj,jk) , 0._wp ) * pt_up1(ji,jj,jk-1) ) * wmask(ji,jj,jk)
            !
            !!clem: not sure we have to update the vertical flux
            ptFw(ji,jj,jk) = ptFw(ji,jj,jk) + ztmp(ji,jj,jk) * e1e2t(ji,jj)  ! update vertical fluxes
         END_3D
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            pt_rhs(ji,jj,jk) = pt_rhs(ji,jj,jk) - ( ztmp(ji,jj,jk) - ztmp(ji,jj,jk+1) )
         END_3D
         ! 
         ! udpdate tridiag matrix with full time step (not half) => clem: is it used?
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
            zwdia(ji,jj,jk) =  1._wp + zDt * ( MAX( wi(ji,jj,jk  ) , 0._wp ) - MIN( wi(ji,jj,jk+1) , 0._wp ) ) / e3t(ji,jj,jk,Kaa)
            zwinf(ji,jj,jk) =          zDt *   MIN( wi(ji,jj,jk  ) , 0._wp )                                   / e3t(ji,jj,jk,Kaa)
            zwsup(ji,jj,jk) =        - zDt *   MAX( wi(ji,jj,jk+1) , 0._wp )                                   / e3t(ji,jj,jk,Kaa)
         END_3D
         !
      ENDIF
      !
   END SUBROUTINE fct_up1_2stp

   
   SUBROUTINE nonosc_org( Kaa, pbef, paa, pbb, pcc, paft, pdt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nonosc  ***
      !!
      !! **  Purpose :   compute monotonic tracer fluxes from the upstream
      !!       scheme and the before field by a nonoscillatory algorithm
      !!
      !! **  Method  : 
      !!       warning : pbef and paft must be masked, but the boundaries
      !!       conditions on the fluxes are not necessary zalezak (1979)
      !!       drange (1995) multi-dimensional forward-in-time and upstream-
      !!       in-space based differencing for fluid
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::   Kaa             ! time level index
      REAL(wp)                            , INTENT(in   ) ::   pdt             ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk)    , INTENT(in   ) ::   pbef            ! before field
      REAL(wp), DIMENSION(T2D(nn_hls),jpk), INTENT(in   ) ::   paft            ! after field
      REAL(wp), DIMENSION(T2D(nn_hls),jpk), INTENT(inout) ::   paa, pbb, pcc   ! monotonic fluxes in the 3 directions
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ikm1         ! local integer
      REAL(wp) ::   zpos, zneg, zbt, zbig                 ! local scalars
      REAL(wp) ::   zup, zdo                              !   -      -
      REAL(wp), DIMENSION(T2D(2),jpk) :: zbetup, zbetdo, zbup, zbdo
      !!----------------------------------------------------------------------
      !
      zbig = HUGE(1._wp)
      zbetup(:,:,jpk) = zbig   ;   zbetdo(:,:,jpk) = zbig

      ! Search local extrema
      ! --------------------
      ! max/min of pbef & paft with large negative/positive value (-/+zbig) inside land
      DO_3D( 2, 2, 2, 2, 1, jpk )
         IF( tmask(ji,jj,jk) == 1._wp ) THEN
            zbup(ji,jj,jk) = MAX( pbef(ji,jj,jk), paft(ji,jj,jk) )
            zbdo(ji,jj,jk) = MIN( pbef(ji,jj,jk), paft(ji,jj,jk) )
         ELSE
            zbup(ji,jj,jk) = -zbig
            zbdo(ji,jj,jk) =  zbig
         ENDIF
      END_3D

      DO jk = 1, jpkm1
         ikm1 = MAX(jk-1,1)
         DO_2D( 1, 1, 1, 1 )

            ! search maximum in neighbourhood
            zup = MAX(  zbup(ji  ,jj  ,jk  ),   &
               &        zbup(ji-1,jj  ,jk  ), zbup(ji+1,jj  ,jk  ),   &
               &        zbup(ji  ,jj-1,jk  ), zbup(ji  ,jj+1,jk  ),   &
               &        zbup(ji  ,jj  ,ikm1), zbup(ji  ,jj  ,jk+1)  )

            ! search minimum in neighbourhood
            zdo = MIN(  zbdo(ji  ,jj  ,jk  ),   &
               &        zbdo(ji-1,jj  ,jk  ), zbdo(ji+1,jj  ,jk  ),   &
               &        zbdo(ji  ,jj-1,jk  ), zbdo(ji  ,jj+1,jk  ),   &
               &        zbdo(ji  ,jj  ,ikm1), zbdo(ji  ,jj  ,jk+1)  )

            ! positive part of the flux
            zpos = MAX( 0._wp, paa(ji-1,jj  ,jk  ) ) - MIN( 0._wp, paa(ji  ,jj  ,jk  ) )   &
               & + MAX( 0._wp, pbb(ji  ,jj-1,jk  ) ) - MIN( 0._wp, pbb(ji  ,jj  ,jk  ) )   &
               & + MAX( 0._wp, pcc(ji  ,jj  ,jk+1) ) - MIN( 0._wp, pcc(ji  ,jj  ,jk  ) )

            ! negative part of the flux
            zneg = MAX( 0._wp, paa(ji  ,jj  ,jk  ) ) - MIN( 0._wp, paa(ji-1,jj  ,jk  ) )   &
               & + MAX( 0._wp, pbb(ji  ,jj  ,jk  ) ) - MIN( 0._wp, pbb(ji  ,jj-1,jk  ) )   &
               & + MAX( 0._wp, pcc(ji  ,jj  ,jk  ) ) - MIN( 0._wp, pcc(ji  ,jj  ,jk+1) )

            ! up & down beta terms
            zbt = e1e2t(ji,jj) * e3t(ji,jj,jk,Kaa) / pdt
            IF( zup /= -zbig .AND. zpos /= 0._wp ) THEN   ;   zbetup(ji,jj,jk) = ( zup - paft(ji,jj,jk) ) / zpos * zbt
            ELSE                                          ;   zbetup(ji,jj,jk) = zbig
            ENDIF
            IF( zdo /=  zbig .AND. zneg /= 0._wp ) THEN   ;   zbetdo(ji,jj,jk) = ( paft(ji,jj,jk) - zdo ) / zneg * zbt
            ELSE                                          ;   zbetdo(ji,jj,jk) = zbig
            ENDIF
         END_2D
      END DO

      ! 3. monotonic flux in the i, j and k direction (paa, pbb and pcc)
      ! ----------------------------------------------------------------
      DO_3D( 1, 0, 1, 0, 1, jpkm1 )
         IF( paa(ji,jj,jk) > 0._wp ) THEN   ;   paa(ji,jj,jk) = paa(ji,jj,jk) * MIN( 1._wp, zbetdo(ji,jj,jk), zbetup(ji+1,jj,jk) )
         ELSE                               ;   paa(ji,jj,jk) = paa(ji,jj,jk) * MIN( 1._wp, zbetup(ji,jj,jk), zbetdo(ji+1,jj,jk) )
         ENDIF
         IF( pbb(ji,jj,jk) > 0._wp ) THEN   ;   pbb(ji,jj,jk) = pbb(ji,jj,jk) * MIN( 1._wp, zbetdo(ji,jj,jk), zbetup(ji,jj+1,jk) )
         ELSE                               ;   pbb(ji,jj,jk) = pbb(ji,jj,jk) * MIN( 1._wp, zbetup(ji,jj,jk), zbetdo(ji,jj+1,jk) )
         ENDIF
      END_3D
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         IF( pcc(ji,jj,jk+1) > 0._wp ) THEN ;   pcc(ji,jj,jk+1) = pcc(ji,jj,jk+1) * MIN( 1._wp, zbetdo(ji,jj,jk+1), zbetup(ji,jj,jk) )
         ELSE                               ;   pcc(ji,jj,jk+1) = pcc(ji,jj,jk+1) * MIN( 1._wp, zbetup(ji,jj,jk+1), zbetdo(ji,jj,jk) )
         ENDIF
      END_3D
      !
   END SUBROUTINE nonosc_org

   SUBROUTINE nonosc( Kaa, pbef, paa, pbb, pcc, paft, pdt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nonosc  ***
      !!
      !!              **  memory optimization version  **
      !!                   **  nn_hls = 2 only  **
      !!
      !! **  Purpose :   compute monotonic tracer fluxes from the upstream
      !!       scheme and the before field by a nonoscillatory algorithm
      !!
      !! **  Method  :   ... ???
      !!       warning : pbef and paft must be masked, but the boundaries
      !!       conditions on the fluxes are not necessary zalezak (1979)
      !!       drange (1995) multi-dimensional forward-in-time and upstream-
      !!       in-space based differencing for fluid
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::   Kaa             ! time level index
      REAL(wp)                            , INTENT(in   ) ::   pdt             ! tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk)    , INTENT(in   ) ::   pbef            ! masked before field
      REAL(wp), DIMENSION(T2D(nn_hls),jpk), INTENT(in   ) ::   paft            ! masked after  field
      REAL(wp), DIMENSION(T2D(nn_hls),jpk), INTENT(inout) ::   paa, pbb, pcc   ! monotonic fluxes in the 3 directions
      !
      INTEGER  ::   ji, jj, jk, jl        ! dummy loop indices
      INTEGER  ::   ikm1, ik, ikp1, iik   ! local integer
      !
      REAL(wp) ::   zpos, zneg, zbig, zsmall, zup, zdo, zbt, zland, z1_Dt   ! local scalars
      REAL(wp) ::   zau, zbu, zcu, zav, zbv, zcv, za, zb, zc, zcoef         !   -      -
      !
      REAL(wp), DIMENSION(T2D(2),3) ::   zbup, zbdo
      REAL(wp), DIMENSION(T2D(1),3) ::   zbetup, zbetdo
      !!----------------------------------------------------------------------
      !
      zbig   = HUGE(1._wp)
      !!zsmall = TINY(1._wp)
      
!!gm NB   pbef = ts(Kbb) |
!!        paft = zwi     |===>>>> both are masked before CALL nonosc 
!!                                can be removed from zbup, zbdd computation

! algorithm :
!
! first level :
!   ===>>>  zb.      computed at level 1 and 2
!   ===>>>  zbet.    computed at level 1 only  and use zbup/do at level 1 and 2 only, ikm1 = ik = 1
!   ===>>>  paa,pbb  update at level 1   using zbet. at 1 only
!           pcc      update at level 1   using zbet. at 1 only, zbet.(ikm1) = zbet.(1).
!
! interior levels jk :
!   ===>>>  zb.      computed at level ikp1       (level ikm1 ik are known )
!   ===>>>  zbet.    computed at level ikk only   using zbup/do at ikm1,ik,ikp1
!   ===>>>  paa,pbb  update at level jk   using zbetup/do at ikk only
!           pcc      update at level jk   using zbetup/do at ikk and ikkm1

      z1_Dt = 1._wp / pdt

      DO jk = 1, jpkm1
         !
         !        !==  set slices and compute zbup/zbdo  ==!   (zbup/do = max/min of pbef and paft with a large +/- value inland)
         !
         IF( jk ==  1 ) THEN          ! set indices - start with (3-1-2)
            !
            ikm1 = 3    ;   ik  = 1     ;   ikp1 = 2
            !
            DO_2D( 2, 2, 2, 2 )
               ! 1st option
!!$               zland = zbig * ( 1._wp - tmask(ji,jj,1) )
!!$               zbup(ji,jj,ik) = MAX( pbef(ji,jj,1) - zland , paft(ji,jj,1) - zland  )
!!$               zbdo(ji,jj,ik) = MIN( pbef(ji,jj,1) + zland , paft(ji,jj,1) + zland  )
!!$
!!$               zland = zbig * ( 1._wp - tmask(ji,jj,2) )
!!$               zbup(ji,jj,ikp1) = MAX( pbef(ji,jj,2) - zland , paft(ji,jj,2) - zland  )
!!$               zbdo(ji,jj,ikp1) = MIN( pbef(ji,jj,2) + zland , paft(ji,jj,2) + zland  )

               ! 2nd option
               zbup(ji,jj,ik)   = MERGE( MAX( pbef(ji,jj,1), paft(ji,jj,1) ), -zbig, tmask(ji,jj,1) == 1._wp )
               zbdo(ji,jj,ik)   = MERGE( MIN( pbef(ji,jj,1), paft(ji,jj,1) ),  zbig, tmask(ji,jj,1) == 1._wp )

               zbup(ji,jj,ikp1) = MERGE( MAX( pbef(ji,jj,2), paft(ji,jj,2) ), -zbig, tmask(ji,jj,2) == 1._wp )
               zbdo(ji,jj,ikp1) = MERGE( MIN( pbef(ji,jj,2), paft(ji,jj,2) ),  zbig, tmask(ji,jj,2) == 1._wp )

               zbup(ji,jj,ikm1) = zbup(ji,jj,ik)
               zbdo(ji,jj,ikm1) = zbdo(ji,jj,ik)
               !
            END_2D                    
            !
         ELSEIF( jk == jpkm1 ) THEN   ! swap indices
            !
            ikm1 = ik   ;   ik = ikp1 ! ikp1 = jpkm1 already ==> keep this value
            !
         ELSE                         ! swap indices
            !
            iik = ikm1
            ikm1 = ik   ;   ik = ikp1   ;   ikp1 = iik
            !
            DO_2D( 2, 2, 2, 2 )
               ! 1st option
!!$               zland = zbig * ( 1._wp - tmask(ji,jj,jk+1) )
!!$               zbup(ji,jj,ikp1) = MAX( pbef(ji,jj,jk+1) - zland , paft(ji,jj,jk+1) - zland  )
!!$               zbdo(ji,jj,ikp1) = MIN( pbef(ji,jj,jk+1) + zland , paft(ji,jj,jk+1) + zland  )
               ! 2nd option
               zbup(ji,jj,ikp1) = MERGE( MAX( pbef(ji,jj,jk+1), paft(ji,jj,jk+1) ), -zbig, tmask(ji,jj,jk+1) == 1._wp )
               zbdo(ji,jj,ikp1) = MERGE( MIN( pbef(ji,jj,jk+1), paft(ji,jj,jk+1) ),  zbig, tmask(ji,jj,jk+1) == 1._wp )
            END_2D                    
            !
         ENDIF
         !
         !        !==  compute the beta term  ==!   (zbetup/do)
         !
         DO_2D( 1, 1, 1, 1 )
            !                                            ! search maximum in neighbourhood
            zup = MAX(  zbup(ji  ,jj  ,ik  ),                         &
               &        zbup(ji-1,jj  ,ik  ), zbup(ji+1,jj  ,ik  ),   &
               &        zbup(ji  ,jj-1,ik  ), zbup(ji  ,jj+1,ik  ),   &
               &        zbup(ji  ,jj  ,ikm1), zbup(ji  ,jj  ,ikp1)  )
            !                                            ! search minimum in neighbourhood
            zdo = MIN(  zbdo(ji  ,jj  ,ik  ),                         &
               &        zbdo(ji-1,jj  ,ik  ), zbdo(ji+1,jj  ,ik  ),   &
               &        zbdo(ji  ,jj-1,ik  ), zbdo(ji  ,jj+1,ik  ),   &
               &        zbdo(ji  ,jj  ,ikm1), zbdo(ji  ,jj  ,ikp1)  )
            !                                            ! positive part of the flux
            zpos = MAX( 0., paa(ji-1,jj  ,jk  ) ) - MIN( 0., paa(ji  ,jj  ,jk  ) )   &
               & + MAX( 0., pbb(ji  ,jj-1,jk  ) ) - MIN( 0., pbb(ji  ,jj  ,jk  ) )   &
               & + MAX( 0., pcc(ji  ,jj  ,jk+1) ) - MIN( 0., pcc(ji  ,jj  ,jk  ) )
            !                                            ! negative part of the flux
            zneg = MAX( 0., paa(ji  ,jj  ,jk  ) ) - MIN( 0., paa(ji-1,jj  ,jk  ) )   &
               & + MAX( 0., pbb(ji  ,jj  ,jk  ) ) - MIN( 0., pbb(ji  ,jj-1,jk  ) )   &
               & + MAX( 0., pcc(ji  ,jj  ,jk  ) ) - MIN( 0., pcc(ji  ,jj  ,jk+1) )
            !                                            ! up & down beta terms
            zbt = e1e2t(ji,jj) * e3t(ji,jj,jk,Kaa) * z1_Dt
            !
            IF( zup /= -zbig .AND. zpos /= 0._wp ) THEN   ;   zbetup(ji,jj,ik) = ( zup - paft(ji,jj,jk) ) / zpos * zbt
            ELSE                                          ;   zbetup(ji,jj,ik) = zbig
            ENDIF
            IF( zdo /=  zbig .AND. zneg /= 0._wp ) THEN   ;   zbetdo(ji,jj,ik) = ( paft(ji,jj,jk) - zdo ) / zneg * zbt
            ELSE                                          ;   zbetdo(ji,jj,ik) = zbig
            ENDIF
            ! better coding but some compilers may not like it (because of division by 0)
!!$            zbetup(ji,jj,ik) = MERGE( ( zup            - paft(ji,jj,jk) ) / zpos * zbt, zbig, zup /= -zbig .AND. zpos /= 0._wp )
!!$            zbetdo(ji,jj,ik) = MERGE( ( paft(ji,jj,jk) - zdo            ) / zneg * zbt, zbig, zdo /=  zbig .AND. zneg /= 0._wp )
            !
         END_2D
         !
         !
         !        !==  monotonic flux  ==!
         !
         DO_2D( 1, 0, 1, 0 )                ! update horizontal flux (paa,pbb)
            ! 1st option
!!$            zau = MIN( 1._wp, zbetdo(ji+1,jj,ik), zbetup(ji,jj,ik) )    ! i-direction
!!$            zbu = MIN( 1._wp, zbetup(ji+1,jj,ik), zbetdo(ji,jj,ik) )
!!$            zcu =         ( 0.5_wp  + SIGN( 0.5_wp , paa(ji,jj,jk) ) )
!!$            zcoef = ( zcu * zau + ( 1._wp - zcu ) * zbu )
!!$            
!!$            paa(ji,jj,jk) = paa(ji,jj,jk) * zcoef
!!$
!!$            zav = MIN( 1._wp, zbetdo(ji,jj+1,ik), zbetup(ji,jj,ik) )   ! j-direction
!!$            zbv = MIN( 1._wp, zbetup(ji,jj+1,ik), zbetdo(ji,jj,ik) )
!!$            zcv =         ( 0.5_wp  + SIGN( 0.5_wp , pbb(ji,jj,jk) ) )
!!$            zcoef = ( zcv * zav + ( 1._wp - zcv ) * zbv )
!!$
!!$            pbb(ji,jj,jk) = pbb(ji,jj,jk) * zcoef
            !
            ! 2nd option
            zcoef = MERGE( MIN( 1._wp, zbetdo(ji,jj,ik), zbetup(ji+1,jj,ik) ), &
               &           MIN( 1._wp, zbetup(ji,jj,ik), zbetdo(ji+1,jj,ik) ), &
               &           paa(ji,jj,jk) > 0._wp )
            paa(ji,jj,jk) = paa(ji,jj,jk) * zcoef

            zcoef = MERGE( MIN( 1._wp, zbetdo(ji,jj,ik), zbetup(ji,jj+1,ik) ), &
               &           MIN( 1._wp, zbetup(ji,jj,ik), zbetdo(ji,jj+1,ik) ), &
               &           pbb(ji,jj,jk) > 0._wp )
            pbb(ji,jj,jk) = pbb(ji,jj,jk) * zcoef           
            !
         END_2D

         IF( jk /= 1 ) THEN                   ! update vertical flux (pcc) (except at the surface where pcc is 0)
            DO_2D( 0, 0, 0, 0 )
               ! 1st option
!!$               za = MIN( 1._wp, zbetdo(ji,jj,ikm1), zbetup(ji,jj,ik) )
!!$               zb = MIN( 1._wp, zbetup(ji,jj,ikm1), zbetdo(ji,jj,ik) )
!!$               zc =         ( 0.5_wp  + SIGN( 0.5_wp , pcc(ji,jj,jk) ) )
!!$               zcoef = ( zc * za + ( 1._wp - zc ) * zb )
!!$               !
!!$               pcc(ji,jj,jk) = pcc(ji,jj,jk) * zcoef
               !
               ! 2nd option
               zcoef = MERGE( MIN( 1._wp, zbetdo(ji,jj,ik), zbetup(ji,jj,ikm1) ), &
                  &           MIN( 1._wp, zbetup(ji,jj,ik), zbetdo(ji,jj,ikm1) ), &
                  &           pcc(ji,jj,jk) > 0._wp )
               pcc(ji,jj,jk) = pcc(ji,jj,jk) * zcoef

            END_2D
         ENDIF
         !
      END DO    ! jk-loop   
      !
   END SUBROUTINE nonosc


   
   SUBROUTINE interp_4th_cpt_org( pt_in, pt_out )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interp_4th_cpt_org  ***
      !!
      !! **  Purpose :   Compute the interpolation of tracer at w-point
      !!
      !! **  Method  :   4th order compact interpolation
      !!----------------------------------------------------------------------
      REAL(wp),DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pt_in    ! now tracer fields
      REAL(wp),DIMENSION(T2D(0) ,jpk), INTENT(  out) ::   pt_out   ! now tracer field interpolated at w-pts
      !
      INTEGER :: ji, jj, jk   ! dummy loop integers
      REAL(wp),DIMENSION(T2D(0),jpk) :: zwd, zwi, zws, zwrm, zwt
      !!----------------------------------------------------------------------

      DO_3D( 0, 0, 0, 0, 3, jpkm1 )       !==  build the three diagonal matrix  ==!
         zwd (ji,jj,jk) = 4._wp
         zwi (ji,jj,jk) = 1._wp
         zws (ji,jj,jk) = 1._wp
         zwrm(ji,jj,jk) = 3._wp * ( pt_in(ji,jj,jk-1) + pt_in(ji,jj,jk) )
         !
         IF( tmask(ji,jj,jk+1) == 0._wp) THEN   ! Switch to second order centered at bottom
            zwd (ji,jj,jk) = 1._wp
            zwi (ji,jj,jk) = 0._wp
            zws (ji,jj,jk) = 0._wp
            zwrm(ji,jj,jk) = 0.5 * ( pt_in(ji,jj,jk-1) + pt_in(ji,jj,jk) )
         ENDIF
      END_3D
      !
      jk = 2                                    ! Switch to second order centered at top
      DO_2D( 0, 0, 0, 0 )
         zwd (ji,jj,jk) = 1._wp
         zwi (ji,jj,jk) = 0._wp
         zws (ji,jj,jk) = 0._wp
         zwrm(ji,jj,jk) = 0.5 * ( pt_in(ji,jj,jk-1) + pt_in(ji,jj,jk) )
      END_2D
      !
      !                       !==  tridiagonal solve  ==!
      DO_2D( 0, 0, 0, 0 )           ! first recurrence
         zwt(ji,jj,2) = zwd(ji,jj,2)
      END_2D
      DO_3D( 0, 0, 0, 0, 3, jpkm1 )
         zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) /zwt(ji,jj,jk-1)
      END_3D
      !
      DO_2D( 0, 0, 0, 0 )           ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
         pt_out(ji,jj,2) = zwrm(ji,jj,2)
      END_2D
      DO_3D( 0, 0, 0, 0, 3, jpkm1 )
         pt_out(ji,jj,jk) = zwrm(ji,jj,jk) - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *pt_out(ji,jj,jk-1)
      END_3D

      DO_2D( 0, 0, 0, 0 )           ! third recurrence: Xk = (Zk - Sk Xk+1 ) / Tk
         pt_out(ji,jj,jpkm1) = pt_out(ji,jj,jpkm1) / zwt(ji,jj,jpkm1)
      END_2D
      DO_3DS( 0, 0, 0, 0, jpk-2, 2, -1 )
         pt_out(ji,jj,jk) = ( pt_out(ji,jj,jk) - zws(ji,jj,jk) * pt_out(ji,jj,jk+1) ) / zwt(ji,jj,jk)
      END_3D
      !
   END SUBROUTINE interp_4th_cpt_org


   SUBROUTINE interp_4th_cpt( pt_in, pt_out )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interp_4th_cpt  ***
      !!
      !! **  Purpose :   Compute the interpolation of tracer at w-point
      !!
      !! **  Method  :   4th order compact interpolation
      !!----------------------------------------------------------------------
      REAL(wp),DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pt_in    ! field at t-point
      REAL(wp),DIMENSION(T2D(1) ,jpk), INTENT(  out) ::   pt_out   ! field interpolated at w-point
      !
      INTEGER ::   ji, jj, jk   ! dummy loop integers
      INTEGER ::   ikt, ikb     ! local integers
      REAL(wp),DIMENSION(T2D(1),jpk) :: zwd, zwi, zws, zwrm, zwt
      !!----------------------------------------------------------------------
      !
      !                      !==  build the three diagonal matrix & the RHS  ==!
      !
      DO_3D( 1, 1, 1, 1, 3, jpkm1 )    ! interior (from jk=3 to jpk-1)
         zwd (ji,jj,jk) = 3._wp * wmask(ji,jj,jk) + 1._wp                 !       diagonal
         zwi (ji,jj,jk) =         wmask(ji,jj,jk)                         ! lower diagonal
         zws (ji,jj,jk) =         wmask(ji,jj,jk)                         ! upper diagonal
         zwrm(ji,jj,jk) = 3._wp * wmask(ji,jj,jk)                     &   ! RHS
            &           *       ( pt_in(ji,jj,jk) + pt_in(ji,jj,jk-1) )
      END_3D
      !
!!gm
!      SELECT CASE( kbc )               !* boundary condition
!      CASE( np_NH   )   ! Neumann homogeneous at top & bottom
!      CASE( np_CEN2 )   ! 2nd order centered  at top & bottom
!      END SELECT
!!gm
      !
      IF ( ln_isfcav ) THEN            ! set level two values which may not be set in ISF case
         DO_2D( 1, 1, 1, 1 )
            zwd (ji,jj,2) = 1._wp
            zwi (ji,jj,2) = 0._wp
            zws (ji,jj,2) = 0._wp
            zwrm(ji,jj,2) = 0._wp
         END_2D
      END IF
      !
      DO_2D( 1, 1, 1, 1 )              ! 2nd order centered at top & bottom
         ikt = mikt(ji,jj) + 1            ! w-point below the 1st  wet point
         ikb = MAX(mbkt(ji,jj), 2)        !     -   above the last wet point
         !
         zwd (ji,jj,ikt) = 1._wp          ! top
         zwi (ji,jj,ikt) = 0._wp
         zws (ji,jj,ikt) = 0._wp
         zwrm(ji,jj,ikt) = 0.5_wp * ( pt_in(ji,jj,ikt-1) + pt_in(ji,jj,ikt) )
         !
         zwd (ji,jj,ikb) = 1._wp          ! bottom
         zwi (ji,jj,ikb) = 0._wp
         zws (ji,jj,ikb) = 0._wp
         zwrm(ji,jj,ikb) = 0.5_wp * ( pt_in(ji,jj,ikb-1) + pt_in(ji,jj,ikb) )
      END_2D
      !
      !                       !==  tridiagonal solver  ==!
      !
      DO_2D( 1, 1, 1, 1 )           !* 1st recurrence:   Tk = Dk - Ik Sk-1 / Tk-1
         zwt(ji,jj,2) = zwd(ji,jj,2)
      END_2D
      DO_3D( 1, 1, 1, 1, 3, jpkm1 )
         zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) /zwt(ji,jj,jk-1)
      END_3D
      !
      DO_2D( 1, 1, 1, 1 )           !* 2nd recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
         pt_out(ji,jj,2) = zwrm(ji,jj,2)
      END_2D
      DO_3D( 1, 1, 1, 1, 3, jpkm1 )
         pt_out(ji,jj,jk) = zwrm(ji,jj,jk) - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *pt_out(ji,jj,jk-1)
      END_3D

      DO_2D( 1, 1, 1, 1 )           !* 3d recurrence:    Xk = (Zk - Sk Xk+1 ) / Tk
         pt_out(ji,jj,jpkm1) = pt_out(ji,jj,jpkm1) / zwt(ji,jj,jpkm1)
      END_2D
      DO_3DS( 1, 1, 1, 1, jpk-2, 2, -1 )
         pt_out(ji,jj,jk) = ( pt_out(ji,jj,jk) - zws(ji,jj,jk) * pt_out(ji,jj,jk+1) ) / zwt(ji,jj,jk)
      END_3D
      !
   END SUBROUTINE interp_4th_cpt


   SUBROUTINE tridia_solver( pD, pU, pL, pRHS, pt_out, klev, kbnd )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tridia_solver  ***
      !!
      !! **  Purpose :   solve a symmetric 3diagonal system
      !!
      !! **  Method  :   solve M.t_out = RHS(t)  where M is a tri diagonal matrix ( jpk*jpk )
      !!
      !!             ( D_1 U_1  0   0   0  )( t_1 )   ( RHS_1 )
      !!             ( L_2 D_2 U_2  0   0  )( t_2 )   ( RHS_2 )
      !!             (  0  L_3 D_3 U_3  0  )( t_3 ) = ( RHS_3 )
      !!             (        ...          )( ... )   ( ...  )
      !!             (  0   0   0  L_k D_k )( t_k )   ( RHS_k )
      !!
      !!        M is decomposed in the product of an upper and lower triangular matrix.
      !!        The tri-diagonals matrix is given as input 3D arrays:   pD, pU, pL
      !!        (i.e. the Diagonal, the Upper diagonal, and the Lower diagonal).
      !!        The solution is pta.
      !!        The 3d array zwt is used as a work space array.
      !!----------------------------------------------------------------------
      REAL(wp),DIMENSION(T2D(nn_hls),jpk)     , INTENT(in   ) ::   pD, pU, PL    ! 3-diagonal matrix
      REAL(wp),DIMENSION(T2D(nn_hls),jpk)     , INTENT(in   ) ::   pRHS          ! Right-Hand-Side
      REAL(wp),DIMENSION(T2D(nn_hls),jpk)     , INTENT(  out) ::   pt_out        !!gm field at level=F(klev)
      INTEGER                                 , INTENT(in   ) ::   klev          ! =1 pt_out at w-level
      !                                                                          ! =0 pt at t-level
      INTEGER                       , OPTIONAL, INTENT(in   ) ::   kbnd          ! Loop bounds
      INTEGER ::   ji, jj, jk   ! dummy loop integers
      INTEGER ::   kstart, ibnd ! local indices
      REAL(wp),DIMENSION(T2D(1),jpk) ::   zwt   ! 3D work array
      !!----------------------------------------------------------------------
      ibnd = 1
      IF( PRESENT(kbnd) ) ibnd = kbnd
      !
      kstart =  1  + klev
      !
      DO_2D( ibnd, ibnd, ibnd, ibnd )                       !* 1st recurrence:   Tk = Dk - Ik Sk-1 / Tk-1
         zwt(ji,jj,kstart) = pD(ji,jj,kstart)
      END_2D
      DO_3D( ibnd, ibnd, ibnd, ibnd, kstart+1, jpkm1 )
         zwt(ji,jj,jk) = pD(ji,jj,jk) - pL(ji,jj,jk) * pU(ji,jj,jk-1) /zwt(ji,jj,jk-1)
      END_3D
      !
      DO_2D( ibnd, ibnd, ibnd, ibnd )                       !* 2nd recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
         pt_out(ji,jj,kstart) = pRHS(ji,jj,kstart)
      END_2D
      DO_3D( ibnd, ibnd, ibnd, ibnd, kstart+1, jpkm1 )
         pt_out(ji,jj,jk) = pRHS(ji,jj,jk) - pL(ji,jj,jk) / zwt(ji,jj,jk-1) *pt_out(ji,jj,jk-1)
      END_3D

      DO_2D( ibnd, ibnd, ibnd, ibnd )                       !* 3d recurrence:    Xk = (Zk - Sk Xk+1 ) / Tk
         pt_out(ji,jj,jpkm1) = pt_out(ji,jj,jpkm1) / zwt(ji,jj,jpkm1)
      END_2D
      DO_3DS( ibnd, ibnd, ibnd, ibnd, jpk-2, kstart, -1 )
         pt_out(ji,jj,jk) = ( pt_out(ji,jj,jk) - pU(ji,jj,jk) * pt_out(ji,jj,jk+1) ) / zwt(ji,jj,jk)
      END_3D
      !
   END SUBROUTINE tridia_solver

   !!======================================================================
END MODULE traadv_fct
