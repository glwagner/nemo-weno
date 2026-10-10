MODULE traadv
   !!==============================================================================
   !!                       ***  MODULE  traadv  ***
   !! Ocean active tracers:  advection trend
   !!==============================================================================
   !! History :  2.0  !  2005-11  (G. Madec)  Original code
   !!            3.3  !  2010-09  (C. Ethe, G. Madec)  merge TRC-TRA + switch from velocity to transport
   !!            3.6  !  2011-06  (G. Madec)  Addition of Mixed Layer Eddy parameterisation
   !!            3.7  !  2014-05  (G. Madec)  Add 2nd/4th order cases for CEN and FCT schemes
   !!             -   !  2014-12  (G. Madec) suppression of cross land advection option
   !!            3.6  !  2015-06  (E. Clementi) Addition of Stokes drift in case of wave coupling
   !!            4.5  !  2021-04  (G. Madec, S. Techene) add advective velocities as optional arguments
   !!            4.5  !  2022-06  (S. Techene, G, Madec) refactorization to reduce local memory usage
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv_trp   : compute transport for tracer advection
   !!   tra_adv       : compute ocean tracer advection trend
   !!   tra_adv_init  : control the different options of advection scheme
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE dom_oce        ! ocean space and time domain
   ! TEMP: [tiling] This change not necessary after all lbc_lnks removed in the nn_hls = 2 case in tra_adv_fct
   USE domtile
   USE domutl, ONLY : lbnd_ij
   USE sbcwave        ! wave module
   USE sbc_oce        ! surface boundary condition: ocean
   USE traadv_cen     ! centered scheme            (tra_adv_cen  routine)
   USE traadv_fct     ! FCT      scheme            (tra_adv_fct  routine)
   USE traadv_mus     ! MUSCL    scheme            (tra_adv_mus  routine)
   USE traadv_ubs     ! UBS      scheme            (tra_adv_ubs  routine)
   USE traadv_qck     ! QUICKEST scheme            (tra_adv_qck  routine)
   USE traadv_weno    ! WENO     scheme            (tra_adv_weno routine)
   USE tramle         ! Mixed Layer Eddy transport (tra_mle_trp  routine)
   USE ldftra         ! Eddy Induced transport     (ldf_eiv_trp  routine)
   USE ldfslp         ! Lateral diffusion: slopes of neutral surfaces
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! trends manager: tracers
   USE diaptr         ! Poleward heat transport
   USE eosbn2
   USE divhor         ! horizontal divergence
   USE sshwzv         ! vertical velocity and ssh
   USE zdf_oce , ONLY : ln_zad_Aimp
   USE dynadv  , ONLY : ln_dynadv_vec
   USE zdf_oce , ONLY : ln_zad_Aimp
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O module
   USE prtctl         ! Print control
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv        ! called by stprk3_stg.F90
   PUBLIC   tra_adv_init   ! called by nemogcm.F90
   PUBLIC   tra_adv_trp    ! called by stprk3_stg.F90

   !                            !!* Namelist namtra_adv *
   LOGICAL ::   ln_traadv_OFF    ! no advection on T and S
   LOGICAL ::   ln_traadv_cen    ! centered scheme flag
   INTEGER ::      nn_cen_h, nn_cen_v   ! =2/4 : horizontal and vertical choices of the order of CEN scheme
   LOGICAL ::   ln_traadv_fct    ! FCT scheme flag
   INTEGER ::      nn_fct_h, nn_fct_v   ! =2/4 : horizontal and vertical choices of the order of FCT scheme
   INTEGER ::      nn_fct_imp           ! =1/2 : optimized or accurate treatment of implicit
   LOGICAL ::   ln_traadv_mus    ! MUSCL scheme flag
   LOGICAL ::      ln_mus_ups           ! use upstream scheme in vivcinity of river mouths
   LOGICAL ::   ln_traadv_ubs    ! UBS scheme flag
   INTEGER ::      nn_ubs_v             ! =2/4 : vertical choice of the order of UBS scheme
   LOGICAL ::   ln_traadv_qck    ! QUICKEST scheme flag
   LOGICAL ::   ln_traadv_weno   ! WENO scheme flag
   INTEGER ::      nn_weno_h, nn_weno_v ! =3 : horizontal and vertical choices of the order of WENO scheme

   INTEGER ::   nadv             ! choice of the type of advection scheme
   !                             ! associated indices:
   INTEGER, PARAMETER ::   np_NO_adv  = 0   ! no T-S advection
   INTEGER, PARAMETER ::   np_CEN     = 1   ! 2nd/4th order centered scheme
   INTEGER, PARAMETER ::   np_FCT     = 2   ! 2nd/4th order Flux Corrected Transport scheme
   INTEGER, PARAMETER ::   np_MUS     = 3   ! MUSCL scheme
   INTEGER, PARAMETER ::   np_UBS     = 4   ! 3rd order Upstream Biased Scheme
   INTEGER, PARAMETER ::   np_QCK     = 5   ! QUICK scheme
   INTEGER, PARAMETER ::   np_WEN     = 6   ! WENO scheme
   !
   REAL(wp), PUBLIC            ::   r_stb_thres_tra_v  ! starting Courant number threshold for adaptive implicit vertical advection
   REAL(wp), PUBLIC            ::   r_stb_cstra_tra_v  ! stability constraint for tracer advection (vertical)
   REAL(wp), PUBLIC            ::   r_stb_cstra_tra_h  ! stability constraint for tracer advection (horizontal)

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!                                                     RK3 time-stepping
   !!----------------------------------------------------------------------

   SUBROUTINE tra_adv_trp( kt, kstg, kit000, Kbb, Kmm, Kaa, Krhs, pFu, pFv, pFw )
      !!
      INTEGER                         , INTENT(in   ) ::   kt                  ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   kstg, kit000        ! RK3 stage and init index
      INTEGER                         , INTENT(in   ) ::   Kbb, Kmm, Kaa, Krhs ! time level indices
      REAL(wp), DIMENSION(:,:,:)      , INTENT(inout) ::   pFu, pFv            ! advective transport
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pFw                 !
      !!
      CALL tra_adv_trp_t( kt, kstg, kit000, Kbb, Kmm, Kaa, Krhs, pFu, pFv, lbnd_ij(pFu), pFw )
   END SUBROUTINE tra_adv_trp


   SUBROUTINE tra_adv_trp_t( kt, kstg, kit000, Kbb, Kmm, Kaa, Krhs, pFu, pFv, ktFuv, pFw )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_trp  ***
      !!
      !! ** Purpose :   compute transport for tracer advection 
      !!
      !! ** Method  : - Update pFu/v/w with stoke drift or eiv or mle
      !!----------------------------------------------------------------------
      INTEGER,  DIMENSION(2)                      , INTENT(in   ) ::   ktFuv
      INTEGER                                     , INTENT(in   ) ::   kt                  ! ocean time-step index
      INTEGER                                     , INTENT(in   ) ::   kstg, kit000        ! RK3 stage and init index
      INTEGER                                     , INTENT(in   ) ::   Kbb, Kmm, Kaa, Krhs ! time level indices
      REAL(wp), DIMENSION(AB2D(ktFuv),JPK)        , INTENT(inout) ::   pFu, pFv            ! advective volume transport  [m^3/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk)            , INTENT(inout) ::   pFw                 ! advective volume transport  [m^3/s]
      !
      INTEGER ::   ji, jj, jk        ! dummy loop index
      REAL(wp)::   z_stfp, z_2stfp   ! local scalar
      LOGICAL ::   ll_Fw             ! local logical
#if ! defined key_PSYCLONE_2p5p0
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE :: zFu_cor, zFv_cor   ! 2D workspace
#else
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zFu_cor, zFv_cor   ! 2D workspace
#endif
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_adv_trp')
      !
      ll_Fw = .false.
      !
      IF( ln_dynadv_vec )   ll_Fw = .true.                          ! VIF: activate vertical advective transport computation
      !
      IF( ln_shuman .AND. kstg == 3 ) THEN      ! shuman averaging (stage 3 only)
         ll_Fw = .true.
         !
         !   shuman time filter parameter (1/6)
         z_stfp  = 0.1666666_wp
         z_2stfp = 2._wp * z_stfp 
         !
#if ! defined key_PSYCLONE_2p5p0
         ALLOCATE( zFu_cor(T2D(nn_hls)), zFv_cor(T2D(nn_hls)) )
#endif
         !
         DO_2D( nn_hls, nn_hls-1, nn_hls, nn_hls-1 )
#if defined key_linssh
            zFu_cor(ji,jj) = ( z_2stfp * un_adv(ji,jj) * r1_hu_0(ji,jj)   &
               &             -  z_stfp * (  uu_b(ji,jj,Kaa) + uu_b(ji,jj,Kbb)  ) ) * e2u(ji,jj)
            zFv_cor(ji,jj) = ( z_2stfp * vn_adv(ji,jj) * r1_hv_0(ji,jj)   &
               &             -  z_stfp * (  vv_b(ji,jj,Kaa) + vv_b(ji,jj,Kbb)  ) ) * e1v(ji,jj)
#else
            zFu_cor(ji,jj) = ( z_2stfp * un_adv(ji,jj) * r1_hu_0(ji,jj)   &
               &             -  z_stfp * (  ( 1._wp + r3u(ji,jj,Kaa) ) * uu_b(ji,jj,Kaa)   &
               &                          + ( 1._wp + r3u(ji,jj,Kbb) ) * uu_b(ji,jj,Kbb)  ) &
               &             ) * e2u(ji,jj)
            zFv_cor(ji,jj) = ( z_2stfp * vn_adv(ji,jj) * r1_hv_0(ji,jj)   &
               &             -  z_stfp * (  ( 1._wp + r3v(ji,jj,Kaa) ) * vv_b(ji,jj,Kaa)   &
               &                          + ( 1._wp + r3v(ji,jj,Kbb) ) * vv_b(ji,jj,Kbb)  )&
               &             ) * e1v(ji,jj)
#endif
         END_2D
         !
         DO_3D( nn_hls, nn_hls-1, nn_hls, nn_hls-1, 1, jpkm1 )
            pFu(ji,jj,jk) = (            (1._wp - z_2stfp) * pFu(ji,jj,jk)                                    &
               &          +  z_stfp * (                  e3u(ji,jj,jk,Kbb)*uu(ji,jj,jk,Kbb)                   &
               &                                       + e3u(ji,jj,jk,Kaa)*uu(ji,jj,jk,Kaa)  ) * e2u(ji,jj)   &
               &          + e3u_0(ji,jj,jk) * zFu_cor(ji,jj)                                                  &
               &            ) * umask(ji,jj,jk)
            !
            pFv(ji,jj,jk) = (            (1._wp - z_2stfp) * pFv(ji,jj,jk)                                    &
               &          +  z_stfp * (                  e3v(ji,jj,jk,Kbb)*vv(ji,jj,jk,Kbb)                   &
               &                                       + e3v(ji,jj,jk,Kaa)*vv(ji,jj,jk,Kaa)  ) * e1v(ji,jj)   &
               &          + e3v_0(ji,jj,jk) * zFv_cor(ji,jj)                                                  &
               &            ) * vmask(ji,jj,jk)
         END_3D
         !
         ! Correct fluxes so that vertical integral matches barotropic mode
         ! estimate
         !
#if ! defined key_PSYCLONE_2p5p0
         DEALLOCATE( zFu_cor, zFv_cor )
#endif
         !
      ENDIF
      !
      IF( ln_sdw )  THEN
         ll_Fw = .true.
         !
         DO_3D( nn_hls, nn_hls-1, nn_hls, nn_hls-1, 1, jpkm1 )                   ! At all stages : Add the Stokes Drift
            pFu(ji,jj,jk) = pFu(ji,jj,jk) + e2u(ji,jj) * e3u(ji,jj,jk,Kmm) * usd(ji,jj,jk)
            pFv(ji,jj,jk) = pFv(ji,jj,jk) + e1v(ji,jj) * e3v(ji,jj,jk,Kmm) * vsd(ji,jj,jk)
         END_3D
      ENDIF
      !
      DO_2D( nn_hls, nn_hls-1, nn_hls, nn_hls-1 )
         pFu(ji,jj,jpk) = 0._wp                                                   ! no transport trough the bottom 
         pFv(ji,jj,jpk) = 0._wp
         pFw(ji,jj,jpk) = 0._wp
      END_2D
      !
      IF( kstg == 3 ) THEN                                   ! At stage 3 only
         IF( ln_ldfeiv .AND. .NOT. ln_traldf_triad ) THEN                     ! Add the eiv transport
            ll_Fw = .true.
            CALL ldf_eiv_trp( kt, kit000, pFu, pFv, pFw, Kmm, Krhs )
         ENDIF
         !
         IF( ln_mle    )   THEN                                              ! Add the mle transport
            ll_Fw = .true.
            CALL eos( ts, Kmm, rhd, rhop )                                    ! now in potential density for tra_mle computation
            CALL tra_mle_trp( kt, pFu, pFv, pFw, Kmm )
         ENDIF
      ENDIF
      !
      IF( ll_Fw ) THEN
         !
         CALL wzv( kt, Kbb, Kmm, Kaa, pFu, pFv, ww, np_transport )
         !
      ELSE 
         IF( ln_zad_Aimp .AND. kstg == 3 ) THEN 
            ! Retrieve total vertical velocity for subsequent ww/wi partitioning
            DO_3D( 1, 1, 1, 1, 2, jpkm1 )
               ww(ji,jj,jk) = ww(ji,jj,jk) + wi(ji,jj,jk)
            END_3D
         ENDIF
      ENDIF
      !                                              ! Partition ww/wwi at stage 3 only
      IF( ln_zad_Aimp .AND. kstg == 3 ) CALL wAimp( kt, Kmm, pFu, pFv, ww, wi, np_transport, &
                                             & r_stb_thres_tra_v, r_stb_cstra_tra_v, r_stb_cstra_tra_h, ld_diag=.TRUE. )
      !
      IF ( ll_Fw .OR. (ln_zad_Aimp.AND.kstg == 3) ) THEN
         DO_3D( 1, 1, 1, 1, 1, jpkm1 )
            pFw(ji,jj,jk) = e1e2t(ji,jj) * ww(ji,jj,jk)
         END_3D
      ENDIF
      !
      !
      IF( ln_timing )   CALL timing_stop( 'tra_adv_trp' )
      !
   END SUBROUTINE tra_adv_trp_t

   
   SUBROUTINE tra_adv( kt, Kbb, Kmm, Kaa, pts, Krhs, pau, pav, paw, kstg )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv  ***
      !!
      !! ** Purpose :   compute the ocean tracer advection trend.
      !!
      !! ** Method  : - Update ts(Krhs) with the advective trend following nadv
      !!----------------------------------------------------------------------
      INTEGER                                  , INTENT(in   ) ::   kt                  ! ocean time-step index
      INTEGER                                  , INTENT(in   ) ::   Kbb, Kmm, Kaa, Krhs ! time level indices
      INTEGER                                  , INTENT(in   ) ::   kstg                ! optional stage indicator
      REAL(wp), DIMENSION(jpi,jpj,jpk)         , INTENT(in   ) ::   pau, pav, paw       ! advective transport  [m^3/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) ::   pts                 ! active tracers and RHS of tracer equation
      !
      INTEGER ::   ji, jj, jk   ! dummy loop index
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdt, ztrds
      ! TEMP: [tiling] This change not necessary after all lbc_lnks removed in the nn_hls = 2 case in tra_adv_fct
      LOGICAL ::   lskip
      LOGICAL ::   ll_dofct
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_adv')
      !
      lskip    = .FALSE.
      ll_dofct = .TRUE.

      ! FCT at last stage
      IF( kstg/=3 ) ll_dofct = .FALSE.

      ! TEMP: [tiling] These changes not necessary after all lbc_lnks removed in the nn_hls = 2 case in tra_adv_fct
      IF( l_istiled .AND. nadv == np_FCT .AND. ll_dofct )  THEN
         IF( ntile == 1 ) THEN
            CALL dom_tile_stop( ldhold=.TRUE. )
         ELSE
            lskip = .TRUE.
         ENDIF
      ENDIF
      !
      IF( .NOT. lskip ) THEN
         !                                         !==  effective advective transport  ==!
         !
         IF( l_trdtra )   THEN                    !* Save ta and sa trends
            ALLOCATE( ztrdt(T2D(0),jpk), ztrds(T2D(0),jpk) )
            ztrdt(:,:,:) = pts(T2D(0),:,jp_tem,Krhs)
            ztrds(:,:,:) = pts(T2D(0),:,jp_sal,Krhs)
         ENDIF
         !
         !
         SELECT CASE ( nadv )                      !==  compute advection trend and add it to general trend  ==!
         !
         CASE ( np_CEN )                                 ! Centered scheme : 2nd / 4th order
            CALL tra_adv_cen    ( kt, nit000, 'TRA',      pau, pav, paw,      Kmm,      pts, jpts, Krhs, nn_cen_h, nn_cen_v )
         CASE ( np_FCT )                                 ! FCT scheme      : 2nd / 4th order
            IF ( ll_dofct ) THEN
               CALL tra_adv_fct ( kt, nit000, 'TRA', rDt, pau, pav, paw, Kbb, Kmm, Kaa, pts, jpts, Krhs, nn_fct_h, nn_fct_v, &
                  &                                                                                                nn_fct_imp )
            ELSE
               CALL tra_adv_cen ( kt, nit000, 'TRA',      pau, pav, paw,      Kmm,      pts, jpts, Krhs, nn_fct_h, nn_fct_v )
            ENDIF
         CASE ( np_MUS )                                 ! MUSCL
            CALL tra_adv_mus    ( kt, nit000, 'TRA', rDt, pau, pav, paw, Kbb, Kmm,      pts, jpts, Krhs, ln_mus_ups         )
         CASE ( np_UBS )                                 ! UBS
            CALL tra_adv_ubs    ( kt, nit000, 'TRA', rDt, pau, pav, paw, Kbb, Kmm,      pts, jpts, Krhs, nn_ubs_v           )
         CASE ( np_QCK )                                 ! QUICKEST
            CALL tra_adv_qck    ( kt, nit000, 'TRA', rDt, pau, pav, paw, Kbb, Kmm,      pts, jpts, Krhs                     )
         CASE ( np_WEN )                                 ! WENO
            CALL tra_adv_weno   ( kt, nit000, 'TRA',      pau, pav, paw,      Kmm,      pts, jpts, Krhs, nn_weno_h, nn_weno_v )
         !
         END SELECT
         !
         CALL iom_put( "uocetr_eff", pau )                                        ! output effective transport [m^3/s]
         CALL iom_put( "vocetr_eff", pav )
         CALL iom_put( "wocetr_eff", paw )
         !
!!gm ???
         IF( l_diaptr ) CALL dia_ptr( kt, Kmm, pav(:,:,:) )                          ! diagnose the effective MSF
!!gm ???
         !
         !
         IF( l_trdtra )   THEN                      ! save the advective trends for further diagnostics
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
               ztrdt(ji,jj,jk) = pts(ji,jj,jk,jp_tem,Krhs) - ztrdt(ji,jj,jk)
               ztrds(ji,jj,jk) = pts(ji,jj,jk,jp_sal,Krhs) - ztrds(ji,jj,jk)
            END_3D
            CALL trd_tra( kt, Kmm, Krhs, 'TRA', jp_tem, jptra_adv, ztrdt )
            CALL trd_tra( kt, Kmm, Krhs, 'TRA', jp_sal, jptra_adv, ztrds )
            DEALLOCATE( ztrdt, ztrds )
         ENDIF

         ! TEMP: [tiling] This change not necessary after all lbc_lnks removed in the nn_hls = 2 case in tra_adv_fct
         IF( ln_tile .AND. .NOT. l_istiled ) CALL dom_tile_start( ldhold=.TRUE. )
         !
      ENDIF
      !                                              ! print mean trends (used for debugging)
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=pts(:,:,:,jp_tem,Krhs), clinfo1=' adv  - Ta: ', mask1=tmask, &
         &                                  tab3d_2=pts(:,:,:,jp_sal,Krhs), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( ln_timing )   CALL timing_stop( 'tra_adv' )
      !
   END SUBROUTINE tra_adv


   SUBROUTINE tra_adv_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_init  ***
      !!
      !! ** Purpose :   Control the consistency between namelist options for
      !!              tracer advection schemes and set nadv
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ios   ! Local integers
      !
      NAMELIST/namtra_adv/ ln_traadv_OFF,                                    &   ! No advection
         &                 ln_traadv_cen , nn_cen_h, nn_cen_v,               &   ! CEN
         &                 ln_traadv_fct , nn_fct_h, nn_fct_v, nn_fct_imp,   &   ! FCT
         &                 ln_traadv_mus , ln_mus_ups,                       &   ! MUSCL
         &                 ln_traadv_ubs ,           nn_ubs_v,               &   ! UBS
         &                 ln_traadv_qck ,                                   &   ! QCK
         &                 ln_traadv_weno, nn_weno_h, nn_weno_v                  ! WENO
      !!----------------------------------------------------------------------
      !
      !                                !==  Namelist  ==!
      READ_NML_REF(numnam,namtra_adv)
      READ_NML_CFG(numnam,namtra_adv)
      IF(lwm) WRITE( numond, namtra_adv )
      !
      IF( ln_mus_ups .AND. .NOT. ln_rnf ) THEN
         CALL ctl_warn( 'tra_adv_init: with no runoff at all, ln_mus_ups is set to .FALSE. ')
         ln_mus_ups = .FALSE.   ! must be set for all processes, so out of the "IF(lwp)"
      ENDIF
      !
      IF(lwp) THEN                           ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_adv_init : choice/control of the tracer advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_adv : chose a advection scheme for tracers'
         WRITE(numout,*) '      No advection on T & S                     ln_traadv_OFF = ', ln_traadv_OFF
         WRITE(numout,*) '      centered scheme                           ln_traadv_cen = ', ln_traadv_cen
         WRITE(numout,*) '            horizontal 2nd/4th order               nn_cen_h   = ', nn_cen_h
         WRITE(numout,*) '            vertical   2nd/4th order               nn_cen_v   = ', nn_cen_v
         WRITE(numout,*) '      Flux Corrected Transport scheme           ln_traadv_fct = ', ln_traadv_fct
         WRITE(numout,*) '            horizontal 2nd/4th order               nn_fct_h   = ', nn_fct_h
         WRITE(numout,*) '            vertical   2nd/4th order               nn_fct_v   = ', nn_fct_v
         WRITE(numout,*) '            implicit   optimized(1)/accurate(2)    nn_fct_imp = ', nn_fct_imp
         WRITE(numout,*) '      MUSCL scheme                              ln_traadv_mus = ', ln_traadv_mus
         WRITE(numout,*) '            + upstream scheme near river mouths    ln_mus_ups = ', ln_mus_ups
         WRITE(numout,*) '      UBS scheme                                ln_traadv_ubs = ', ln_traadv_ubs
         WRITE(numout,*) '            vertical   2nd/4th order               nn_ubs_v   = ', nn_ubs_v
         WRITE(numout,*) '      QUICKEST scheme                           ln_traadv_qck = ', ln_traadv_qck
         WRITE(numout,*) '      WENO scheme                              ln_traadv_weno = ', ln_traadv_weno
         WRITE(numout,*) '            horizontal order                       nn_weno_h  = ', nn_weno_h
         WRITE(numout,*) '            vertical   order                       nn_weno_v  = ', nn_weno_v
      ENDIF
      !
      !                                !==  Parameter control & set nadv ==!
      ioptio = 0
      IF( ln_traadv_OFF ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_NO_adv   ;   ENDIF
      IF( ln_traadv_cen ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_CEN      ;   ENDIF
      IF( ln_traadv_fct ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_FCT      ;   ENDIF
      IF( ln_traadv_mus ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_MUS      ;   ENDIF
      IF( ln_traadv_ubs ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_UBS      ;   ENDIF
      IF( ln_traadv_qck ) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_QCK      ;   ENDIF
      IF( ln_traadv_weno) THEN   ;   ioptio = ioptio + 1   ;   nadv = np_WEN      ;   ENDIF
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'tra_adv_init: Choose ONE advection option in namelist namtra_adv' )
      !
      IF( ln_traadv_cen .AND. ( nn_cen_h /= 2 .AND. nn_cen_h /= 4 )   &          ! Centered
                        .AND. ( nn_cen_v /= 2 .AND. nn_cen_v /= 4 )   ) THEN
        CALL ctl_stop( 'tra_adv_init: CEN scheme, choose 2nd or 4th order' )
      ENDIF
      IF( ln_traadv_fct .AND. ( nn_fct_h /= 2 .AND. nn_fct_h /= 4 )   &          ! FCT
         &              .AND. ( nn_fct_v /= 2 .AND. nn_fct_v /= 4 ) .AND. ( nn_fct_imp /= 1 .AND. nn_fct_imp /= 2 ) ) THEN
        CALL ctl_stop( 'tra_adv_init: FCT scheme, choose 2nd or 4th order, and optimized or accurate treatment of implicit' )
      ENDIF
      ! TEMP: [tiling] This change not necessary after all lbc_lnks removed in the nn_hls = 2 case in tra_adv_fct
      IF( ln_traadv_fct .AND. ln_tile ) THEN
         CALL ctl_warn( 'tra_adv_init: FCT scheme does not yet work with tiling' )
      ENDIF
      IF( ln_traadv_ubs .AND. ( nn_ubs_v /= 2 .AND. nn_ubs_v /= 4 )   ) THEN     ! UBS
         CALL ctl_stop( 'tra_adv_init: UBS scheme, choose 2nd or 4th order' )
      ENDIF
      IF( ln_traadv_ubs .AND. nn_ubs_v == 4 ) THEN
         CALL ctl_warn( 'tra_adv_init: UBS scheme, only 2nd FCT scheme available on the vertical. It will be used' )
      ENDIF
      IF( ln_traadv_ubs .AND. ln_zad_Aimp )   CALL ctl_stop( 'tra_adv_init: UBS scheme, vertical implicit is not coded' )
      IF( ln_traadv_weno .AND. ( nn_weno_h /= 3 .OR. nn_weno_v /= 3 ) ) THEN     ! WENO
         CALL ctl_stop( 'tra_adv_init: WENO scheme, only 3rd order is available (higher orders need nn_hls > 2)' )
      ENDIF
      IF( ln_isfcav ) THEN                                                       ! ice-shelf cavities
         IF(  ln_traadv_cen .AND. nn_cen_v == 4    .OR.   &                            ! NO 4th order with ISF
            & ln_traadv_fct .AND. nn_fct_v == 4   )   CALL ctl_stop( 'tra_adv_init: 4th order COMPACT scheme not allowed with ISF' )
      ENDIF
      !
      IF( ln_traadv_cen .OR. ln_traadv_fct) THEN
         IF ( nn_fct_h==2 ) THEN ; r_stb_cstra_tra_h = pp_stb_cstra_tra_ce2 ; ELSE ; r_stb_cstra_tra_h = pp_stb_cstra_tra_ce4 ; ENDIF
         IF ( nn_fct_v==2 ) THEN ; r_stb_cstra_tra_v = pp_stb_cstra_tra_ce2 ; ELSE ; r_stb_cstra_tra_v = pp_stb_cstra_tra_co4 ; ENDIF
         IF ( nn_fct_v==2 ) THEN ; r_stb_thres_tra_v = pp_stb_thres_tra_ce2 ; ELSE ; r_stb_thres_tra_v = pp_stb_thres_tra_co4 ; ENDIF
         !
      ELSEIF ( ln_traadv_mus ) THEN
         r_stb_cstra_tra_h = pp_stb_cstra_tra_mus
         r_stb_cstra_tra_v = pp_stb_cstra_tra_mus
         r_stb_thres_tra_v = pp_stb_thres_tra_mus
         !
      ELSEIF ( ln_traadv_ubs ) THEN
         r_stb_cstra_tra_h = pp_stb_cstra_tra_ubs
         IF ( nn_fct_v==2 ) THEN ; r_stb_cstra_tra_v = pp_stb_cstra_tra_ce2 ; ELSE ; r_stb_cstra_tra_v = pp_stb_cstra_tra_co4 ; ENDIF
         IF ( nn_fct_v==2 ) THEN ; r_stb_thres_tra_v = pp_stb_thres_tra_ce2 ; ELSE ; r_stb_thres_tra_v = pp_stb_thres_tra_co4 ; ENDIF
         !
      ELSEIF ( ln_traadv_qck ) THEN
         r_stb_cstra_tra_h = pp_stb_cstra_tra_qck
         r_stb_cstra_tra_v = pp_stb_cstra_tra_qck
         r_stb_thres_tra_v = pp_stb_thres_tra_qck
         !
      ELSEIF ( ln_traadv_weno ) THEN
         r_stb_cstra_tra_h = pp_stb_cstra_tra_weno
         r_stb_cstra_tra_v = pp_stb_cstra_tra_weno
         r_stb_thres_tra_v = pp_stb_thres_tra_weno
         !
      ENDIF
      !
      !                                !==  Print the choice  ==!
      IF(lwp) THEN
         WRITE(numout,*)
         SELECT CASE ( nadv )
         CASE( np_NO_adv  )   ;   WRITE(numout,*) '   ==>>>   NO T-S advection'
         CASE( np_CEN     )   ;   WRITE(numout,*) '   ==>>>   CEN      scheme is used. Horizontal order  : ', nn_cen_h, &
            &                                                                        ' Vertical   order  : ', nn_cen_v
         CASE( np_FCT     )   ;   WRITE(numout,*) '   ==>>>   FCT      scheme is used. Horizontal order  : ', nn_fct_h, &
            &                                                                        ' Vertical   order  : ', nn_fct_v, &
            &                                                                        ' Implicit treatment: ', nn_fct_imp
         CASE( np_MUS     )   ;   WRITE(numout,*) '   ==>>>   MUSCL    scheme is used'
         CASE( np_UBS     )   ;   WRITE(numout,*) '   ==>>>   UBS      scheme is used'
         CASE( np_QCK     )   ;   WRITE(numout,*) '   ==>>>   QUICKEST scheme is used'
         CASE( np_WEN     )   ;   WRITE(numout,*) '   ==>>>   WENO     scheme is used. Horizontal order  : ', nn_weno_h, &
            &                                                                        ' Vertical   order  : ', nn_weno_v
         END SELECT
      ENDIF
      !
      CALL tra_mle_init            !== initialisation of the Mixed Layer Eddy parametrisation (MLE)  ==!
      !
   END SUBROUTINE tra_adv_init

  !!======================================================================
END MODULE traadv
