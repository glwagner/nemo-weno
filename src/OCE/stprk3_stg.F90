MODULE stprk3_stg
   !!======================================================================
   !!                       ***  MODULE stp_stg  ***
   !! Time-stepping   : manager of the ocean, tracer and ice time stepping
   !!                   using a 3rd order Runge-Kutta  with fixed or quasi-eulerian coordinate
   !!======================================================================
   !! History :  4.2  !  2021-01  (S. Techene, G. Madec, N. Ducousso, F. Lemarie)  Original code
   !!   NEMO     5.0  !  2024-01  (S. Techene, G. Madec)  restructuration for Shuman averaging  
   !!             -   !  2024-01  (S. Techene, G. Madec)  optimized 1st stage 3D RHS due to 2D RHS optimization 
   !!----------------------------------------------------------------------
#if defined key_qco   ||   defined key_linssh
   !!----------------------------------------------------------------------
   !!   'key_qco'                        Quasi-Eulerian vertical coordinate
   !!                          OR
   !!   'key_linssh                       Fixed in time vertical coordinate
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!    stp_stg  : NEMO 3rd order Runge-Kutta stage with qco or linssh
   !!----------------------------------------------------------------------
   USE step_oce       ! time stepping used modules
   USE domqco         ! quasi-eulerian coordinate      (dom_qco_r3c routine)
   USE bdydyn         ! ocean open boundary conditions (define bdy_dyn)
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE tramle         ! ML eddy induced transport (tra_adv_mle  routine)
!!smg
   USE trdtra,   ONLY : trd_tra_tot ! total T,S thickness weighted tendency 
!!smg
# if defined key_top
   USE trcstp_rk3
# endif
# if defined key_agrif
   USE agrif_oce_interp
# endif
   !
   USE prtctl         ! print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp_stg   ! called by nemogcm.F90

   INTEGER,  PUBLIC, PARAMETER ::   np_LIN = 0   ! linear interpolation of barotropic mode at each stage: ORCA1 OK
   INTEGER,  PUBLIC, PARAMETER ::   np_IMP = 1   ! implicitation of barotropic mode at each stage: ORCA2 OK
   INTEGER,  PUBLIC, PARAMETER ::   np_HYB = 2   ! hybrid version linear interpolation of ssh implicitation barotropic velocities
!!smg   
   LOGICAL  ::   lo_trd_tot = .FALSE.   ! flag for total e3t*tracer trend 
!!smg
   INTEGER  ::   n_baro_upd =  np_HYB   ! Forward in time update of barotropic mode
                                        ! at stages 1 & 2

   REAL(wp) ::   r1_2 = 1._wp / 2._wp
   REAL(wp) ::   r1_3 = 1._wp / 3._wp
   REAL(wp) ::   r2_3 = 2._wp / 3._wp
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ssha               ! sea-surface height  at N+1
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ua_b, va_b         ! barotropic velocity at N+1
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   r3ta, r3ua, r3va   ! ssh/h_0 ratio at t,u,v-points at N+1
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   r3fb, r3fa         ! ssh/h_0 ratio at f-point at N and N+1

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE stp_stg( kstg, kstp, Kbb, Kmm, Krhs, Kaa )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp_stg  ***
      !!
      !! ** Purpose : - stage of RK3 time stepping of OCE and TOP
      !!
      !! ** Method  :   input: computed by stp_2D
      !!              ssh             sea surface height  at N, N+1
      !!              (uu_b,vv_b)     barotropic velocity at N, N+1
      !!              (un_adv,vn_adv) barotropic transport from N to N+1
      !!              (uu,vv)_Krhs    3D RHS at 1st stage in Vector Inv. Form only
      !!              ,
      !!              -1- set ssh(Naa) (Naa = N+1/3, N+1/2, or N+1)
      !!              -2- set the Kmm advective velocity (zFu,zFv,zFw) 
      !!                  and vertical velocity (ww,wi)
      !!              -3- Compute the after (Naa) velocity (uu,vv)_Kaa
      !!              -4- Compute the after (Naa) T-S
      !!              -5- 
      !!              -6- 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kstg                        ! RK3 stage
      INTEGER, INTENT(in) ::   kstp, Kbb, Kmm, Krhs, Kaa   ! ocean time-step and time-level indices
      !
      INTEGER  ::   ji, jj, jk, jn, jtile    ! dummy loop indices
      REAL(wp) ::   ze3Tb, ze3Tr, z1_e3t     ! local scalars
      REAL(wp) ::   ze3Sb, ze3Sr             !   -      -
#if ! defined key_PSYCLONE_2p5p0
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   zub, zvb        ! barotropic velocity correction at Kmm
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zFu, zFv, zFw   ! advective transport
#else
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zub, zvb        ! barotropic velocity correction at Kmm
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zFu, zFv, zFw   ! advective transport
#endif
      !! ---------------------------------------------------------------------

      IF( ln_timing )   CALL timing_start('stp_stg')
      !
      IF( kstp == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'stp_stg : Runge Kutta 3rd order at stage ', kstg
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         lo_trd_tot = iom_use('ttrd_tot') .OR. iom_use('strd_tot')
      ENDIF
      !
#if ! defined key_PSYCLONE_2p5p0
      ALLOCATE( zFu(jpi,jpj,jpk), zFv(jpi,jpj,jpk), zFw(jpi,jpj,jpk) )
#endif
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !>>>             Set ∆t and external mode fields at Kaa              <<<
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      SELECT CASE( kstg )
      !                    !---------------!
      CASE ( 1 )           !==  Stage 1  ==!   Kbb = Kmm = N  ;  Kaa = N+1/3
         !                 !---------------!
         !
         ALLOCATE( ssha(jpi,jpj) , ua_b(jpi,jpj) , va_b(jpi,jpj) )
         !
         rDt = r1_3 * rn_Dt            ! set time-step : rn_Dt/3
         r1_Dt = 1._wp / rDt
         !
         ssha(:,:) = ssh (:,:,Kaa)     ! save ssh, uu_b, vv_b at N+1  (computed in dynspg_ts)
         ua_b(:,:) = uu_b(:,:,Kaa)
         va_b(:,:) = vv_b(:,:,Kaa)

         SELECT CASE( n_baro_upd )
         CASE ( np_IMP )
            !                          ! Set ssh and (uu_b,vv_b) at N+1 value
            ssh (:,:,Kaa) = ssha(:,:)
            uu_b(:,:,Kaa) = ua_b(:,:)
            vv_b(:,:,Kaa) = va_b(:,:)
         CASE ( np_LIN )
            !                          ! interpolated ssh and (uu_b,vv_b) at Kaa (N+1/3)
            ssh (:,:,Kaa) = r2_3 * ssh (:,:,Kbb) + r1_3 * ssha(:,:)
            uu_b(:,:,Kaa) = r2_3 * uu_b(:,:,Kbb) + r1_3 * ua_b(:,:)
            vv_b(:,:,Kaa) = r2_3 * vv_b(:,:,Kbb) + r1_3 * va_b(:,:)
         CASE ( np_HYB ) 
            !                          ! hybrid : interpolated ssh at Kaa (N+1/3) and implict (uu_b,vv_b) at N+1 values
            ssh (:,:,Kaa) = r2_3 * ssh (:,:,Kbb) + r1_3 * ssha(:,:)
            uu_b(:,:,Kaa) = ua_b(:,:)
            vv_b(:,:,Kaa) = va_b(:,:)
         END SELECT
         !
         !
         !                     !==  ssh/h0 ratio at Kaa  ==!
         !
         IF( .NOT.lk_linssh ) THEN     ! "after" ssh/h_0 ratio at t,u,v-column computed at N+1 stored in r3.a
            !
            ALLOCATE( r3ta(jpi,jpj) , r3ua(jpi,jpj) , r3va(jpi,jpj) , r3fa(jpi,jpj) , r3fb(jpi,jpj) )
            !
            r3fb(:,:) = r3f(:,:)
            CALL dom_qco_r3c_RK3( ssha, r3ta, r3ua, r3va, r3fa )
            !
            CALL lbc_lnk( 'stp_stg', r3ua, 'U', 1._wp, r3va, 'V', 1._wp, r3fa, 'F', 1._wp )
            !                          !
            SELECT CASE( n_baro_upd )
            CASE ( np_IMP )
               r3t(:,:,Kaa) = r3ta(:,:)   ! at N+1
               r3u(:,:,Kaa) = r3ua(:,:)
               r3v(:,:,Kaa) = r3va(:,:)
            CASE ( np_LIN, np_HYB ) 
               r3t(:,:,Kaa) = r2_3 * r3t(:,:,Kbb) + r1_3 * r3ta(:,:)   ! at N+1/3 (Kaa)
               r3u(:,:,Kaa) = r2_3 * r3u(:,:,Kbb) + r1_3 * r3ua(:,:)
               r3v(:,:,Kaa) = r2_3 * r3v(:,:,Kbb) + r1_3 * r3va(:,:)
            END SELECT
            ! r3f already properly set up                           ! at N     (Kmm)
         ENDIF
         !
         !                 !---------------!
      CASE ( 2 )           !==  Stage 2  ==!   Kbb = N   ;   Kmm = N+1/3   ;   Kaa = N+1/2
         !                 !---------------!
         !
         rDt = r1_2 * rn_Dt            ! set time-step : rn_Dt/2
         r1_Dt = 1._wp / rDt
         !
         SELECT CASE( n_baro_upd )
         CASE ( np_IMP )
            !                          ! set ssh and (uu_b,vv_b) at N+1
            ssh (:,:,Kaa) = ssha(:,:)
            uu_b(:,:,Kaa) = ua_b(:,:)
            vv_b(:,:,Kaa) = va_b(:,:)
            !
            IF( .NOT.lk_linssh ) THEN
               r3t(:,:,Kaa) = r3ta(:,:)    ! at N+1 
               r3u(:,:,Kaa) = r3ua(:,:)
               r3v(:,:,Kaa) = r3va(:,:)
               r3f(:,:)     = r3fa(:,:) 
            ENDIF
         CASE ( np_LIN )
            !                          ! set ssh and (uu_b,vv_b) at N+1/2  (Kaa)
            ssh (:,:,Kaa) = r1_2 * ( ssh (:,:,Kbb) + ssha(:,:) )
            uu_b(:,:,Kaa) = r1_2 * ( uu_b(:,:,Kbb) + ua_b(:,:) )
            vv_b(:,:,Kaa) = r1_2 * ( vv_b(:,:,Kbb) + va_b(:,:) )
            IF( .NOT.lk_linssh ) THEN
               r3t(:,:,Kaa) = r1_2 * ( r3t(:,:,Kbb) + r3ta(:,:) )   ! at N+1/2 (Kaa)
               r3u(:,:,Kaa) = r1_2 * ( r3u(:,:,Kbb) + r3ua(:,:) )
               r3v(:,:,Kaa) = r1_2 * ( r3v(:,:,Kbb) + r3va(:,:) )
               r3f(:,:)     = r2_3 * r3fb(:,:) + r1_3 * r3fa(:,:)   ! at N+1/3 (Kmm)
            ENDIF
         CASE ( np_HYB )
            !                          ! set ssh and (uu_b,vv_b) at N+1/2  (Kaa)
            ssh (:,:,Kaa) = r1_2 * ( ssh (:,:,Kbb) + ssha(:,:) )
            uu_b(:,:,Kaa) = ua_b(:,:)
            vv_b(:,:,Kaa) = va_b(:,:)
            IF( .NOT.lk_linssh ) THEN
               r3t(:,:,Kaa) = r1_2 * ( r3t(:,:,Kbb) + r3ta(:,:) )   ! at N+1/2 (Kaa)
               r3u(:,:,Kaa) = r1_2 * ( r3u(:,:,Kbb) + r3ua(:,:) )
               r3v(:,:,Kaa) = r1_2 * ( r3v(:,:,Kbb) + r3va(:,:) )
               r3f(:,:)     = r2_3 * r3fb(:,:) + r1_3 * r3fa(:,:)   ! at N+1/3 (Kmm)
            ENDIF
         END SELECT
         !
         !                 !---------------!
      CASE ( 3 )           !==  Stage 3  ==!   Kbb = N   ;   Kmm = N+1/2   ;   Kaa = N+1
         !                 !---------------!
         !
         rDt = rn_Dt                   ! set time-step : rn_Dt
         r1_Dt = 1._wp / rDt
         !
         ssh (:,:,Kaa) = ssha(:,:)     ! recover ssh and (uu_b,vv_b) at N + 1
         uu_b(:,:,Kaa) = ua_b(:,:)
         vv_b(:,:,Kaa) = va_b(:,:)
         !
         DEALLOCATE( ssha , ua_b , va_b )
         !
         IF( .NOT.lk_linssh ) THEN
            r3t(:,:,Kaa) = r3ta(:,:)                          ! at N+1   (Kaa)
            r3u(:,:,Kaa) = r3ua(:,:)
            r3v(:,:,Kaa) = r3va(:,:)
            r3f(:,:    ) = r1_2 * ( r3fb(:,:) + r3fa(:,:) )   ! at N+1/2 (Kmm)
            DEALLOCATE( r3ta, r3ua, r3va, r3fb )              ! deallocate all r3. except r3fa which will be
            !                                                 ! saved in r3f at the end of the time integration and then deallocated
            !
         ENDIF
         !
      END SELECT

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !>>>            Dynamic : RHS computation + time-stepping            <<<
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      IF( ln_tile ) CALL dom_tile_start         ! [tiling] DYN tiling loop (2)
      DO jtile = 1, nijtile
         IF( ln_tile ) CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )
         !
         !           !=============================================!
         !           !==   Set fields used in advection at Kmm   ==!
         !           !=============================================!
         !
#if ! defined key_PSYCLONE_2p5p0
         ALLOCATE( zub(T2D(nn_hls)), zvb(T2D(nn_hls)) )
#endif
         !              !- horizontal transport (zFu,zFv) -!   VF for tracers only ; FF for momentum & tracers
         !
         SELECT CASE( n_baro_upd )
         CASE ( np_IMP )
            DO_2D( nn_hls, nn_hls-1, nn_hls, nn_hls-1 )
               zub(ji,jj) = r1_Dt * rn_Dt * un_adv(ji,jj)*r1_hu(ji,jj,Kmm) - uu_b(ji,jj,Kmm)    ! barotropic velocity correction
               zvb(ji,jj) = r1_Dt * rn_Dt * vn_adv(ji,jj)*r1_hv(ji,jj,Kmm) - vv_b(ji,jj,Kmm)
            END_2D
         CASE ( np_LIN, np_HYB )
            DO_2D( nn_hls, nn_hls-1, nn_hls, nn_hls-1 )
               zub(ji,jj) = un_adv(ji,jj)*r1_hu(ji,jj,Kmm) - uu_b(ji,jj,Kmm)    ! barotropic velocity correction
               zvb(ji,jj) = vn_adv(ji,jj)*r1_hv(ji,jj,Kmm) - vv_b(ji,jj,Kmm)
            END_2D
         END SELECT
         !
         DO_3D( nn_hls, nn_hls-1, nn_hls, nn_hls-1, 1, jpkm1 )               ! advective transport
            zFu(ji,jj,jk) = e2u(ji,jj)*e3u(ji,jj,jk,Kmm) * ( uu(ji,jj,jk,Kmm) + zub(ji,jj)*umask(ji,jj,jk) )
            zFv(ji,jj,jk) = e1v(ji,jj)*e3v(ji,jj,jk,Kmm) * ( vv(ji,jj,jk,Kmm) + zvb(ji,jj)*vmask(ji,jj,jk) )
         END_3D

#if ! defined key_PSYCLONE_2p5p0
         DEALLOCATE( zub, zvb )
#endif
         !
         !              !- vertical velocity and transport (ww,wi,zFw) -!
         !
         !                                        !* adaptive-implicit vertical advection is  zero  at stage 1 and stage 2
         IF( ln_zad_Aimp .AND. kstg == 1 )   wi(T2D(1),:) = 0._wp     ! and is update at stage 3 only
         !
         IF( ln_dynadv_vec ) THEN                 !* Vector invariant Form : no use of ww at stage 1 as 3D RHS computed in stp_2D 
            !                                                                zFw used in tracers only and computed in tra_adv_trp
            IF( kstg /= 1 ) THEN
               !                                      ! ww cross-level velocity at Kmm consistent with (uu,vv)
               CALL wzv( kstp, Kbb, Kmm, Kaa, uu(:,:,:,Kmm), vv(:,:,:,Kmm), ww, np_velocity )
               !                                      ! ww / wi Partition at stage 3 only
               IF( ln_zad_Aimp .AND. kstg == 3 )   CALL wAimp( kstp, Kmm, uu(:,:,:,Kmm), vv(:,:,:,Kmm), ww, wi, np_velocity, & 
                                                        &     r_stb_thres_dyn_v, r_stb_cstra_dyn_v, r_stb_cstra_dyn_h )
               !                                      ! ww needed over 2 halo points (no tiling with ln_dynadv_weno)
               IF( ln_dynadv_weno )   CALL lbc_lnk( 'stp_RK3_stg', ww, 'T', 1._wp )
               !
            ENDIF
         ELSE                                     !* Flux Form : set ww, wi at 3rd stage and zFw
            !                                         ! ww cross-level velocity at Kmm consistent with (zFu,zFv)
            CALL wzv( kstp, Kbb, Kmm, Kaa, zFu, zFv, ww, np_transport )
            !                                         ! ww / wi Partition at stage 3 only
            IF( ln_zad_Aimp .AND. kstg == 3 )   CALL wAimp( kstp, Kmm, zFu, zFv, ww, wi, np_transport, & 
                                                        &     r_stb_thres_dyn_v, r_stb_cstra_dyn_v, r_stb_cstra_dyn_h )
            DO_3D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpkm1 )
               zFw(ji,jj,jk) = e1e2t(ji,jj) * ww(ji,jj,jk)
            END_3D
            !
         ENDIF
         !           !============================================!
         !           !==   RHS : part computed at each stages   ==!
         !           !============================================!
         !
         SELECT CASE( kstg )
         !                 !---------------!   Flux Form        : add the missing ADV to the 1st stage 3D RHS
         CASE ( 1 )        !==  Stage 1  ==!
            !              !---------------!   Vector Inv. Form : 1st stage 3D RHS already entirely computed in stp_2D
            !
            !                             !* Flux Form advection only   ==> 3D RHS
            IF( .NOT.ln_dynadv_vec )   CALL dyn_adv( kstp, Kmm, Kmm, uu, vv, Krhs, zFu, zFv, zFw )
            !
         !                 !-------------------!   Flux Form        : HPG + VOR (COR+MET) + ADV 
         CASE ( 2 , 3 )    !==  Stage 2 & 3  ==!
            !              !-------------------!   Vector Inv. Form : HPG + VOR (COR+RVO) +ADV (KEG +ZAD)
            !
            !                             !*  hydrostatic pressure gradient (HPG))  *!   always called FIRST
            CALL    eos    (        ts, Kmm, rhd, rhop )     ! Kmm in situ density anomaly for hpg computation
            !                                                ! and potential density for BGC 
            CALL    dyn_hpg( kstp,      Kmm, uu, vv, Krhs )  ! Hydrostratic Pressure Gradient (HPG)
            !
            !                             !* Coriolis / vorticity       ==> 3D RHS
            CALL    dyn_vor( kstp,      Kmm, uu, vv, Krhs )
            !
            !                             !* advection (VIF or FF)      ==> 3D RHS
            IF( ln_dynadv_vec ) THEN                        ! VIF: only velocities used for momentum advection
               CALL dyn_adv( kstp, Kmm, Kmm, uu, vv, Krhs)
            ELSE                                            ! FF : advective transports used for momentum advection
               CALL dyn_adv( kstp, Kmm, Kmm, uu, vv, Krhs, zFu, zFv, zFw )
            ENDIF
            !
         END SELECT
         !
!!gm ===>>>>>> Probably useless since uu_b(Kaa) will be imposed at the end of stage 1 and 2
!                   but may be necessary in stage 3 due to implicite in dynzdf.
!                   except if my idea for the matrice construction is OK !
!         !                                         ! grad_h of ps          ==> RHS
!         DO_3D( 0, 0, 0, 0, 1, jpkm1 )
!            uu(ji,jj,jk,Krhs) = uu(ji,jj,jk,Krhs) - grav * ( ssh(ji+1,jj  ,Kmm) - ssh(ji,jj,Kmm) )
!            vv(ji,jj,jk,Krhs) = vv(ji,jj,jk,Krhs) - grav * ( ssh(ji  ,jj+1,Kmm) - ssh(ji,jj,Kmm) )
!         END_3D
!!gm

!!mjb assimilation increments are applied at all RK3 stages because tests applying them only at stage 3
!!mjb gave strangely poor results. The additional computational expense is minor.  
!!mjb tra_sbc_RK3 and ssh_asm_div should almost certainly be applied on all RK3 stages 
      IF(  lk_asminc .AND. ln_asmiau ) THEN               ! apply assimilation increment
         IF( ln_dyninc )   CALL dyn_asm_inc( kstp, Kbb, Kmm, uu, vv, Krhs )   ! dynamics   ==> RHS
      ENDIF

         !          !=================================================================!
         !          !==   stage 1 & 2 : time-stepping                               ==!
         !          !==   stage 3     : time-stepping with all remaining RHS trends ==!
         !          !=================================================================!
         !
         SELECT CASE( kstg )
         !                 !-------------------!
         CASE ( 1 , 2 )    !==  Stage 1 & 2  ==!   time stepping
            !              !-------------------!
            !
            IF( ln_dynadv_vec .OR. lk_linssh ) THEN   !* applied on velocity
               DO_3D( 0, 0, 0, 0, 1, jpkm1 )
                  uu(ji,jj,jk,Kaa) = ( uu(ji,jj,jk,Kbb) + rDt * uu(ji,jj,jk,Krhs) ) * umask(ji,jj,jk)
                  vv(ji,jj,jk,Kaa) = ( vv(ji,jj,jk,Kbb) + rDt * vv(ji,jj,jk,Krhs) ) * vmask(ji,jj,jk)
               END_3D
            ELSE                                      !* applied on thickness weighted velocity
               DO_3D( 0, 0, 0, 0, 1, jpkm1 )
#  if defined key_qco
                  uu(ji,jj,jk,Kaa) = (         ( 1._wp + r3u(ji,jj,Kbb) ) * uu(ji,jj,jk,Kbb )  &
                     &                 + rDt * ( 1._wp + r3u(ji,jj,Kmm) ) * uu(ji,jj,jk,Krhs)  )   &
                     &             /           ( 1._wp + r3u(ji,jj,Kaa) ) * umask(ji,jj,jk)
                  vv(ji,jj,jk,Kaa) = (         ( 1._wp + r3v(ji,jj,Kbb) ) * vv(ji,jj,jk,Kbb )  &
                     &                 + rDt * ( 1._wp + r3v(ji,jj,Kmm) ) * vv(ji,jj,jk,Krhs)  )   &
                     &             /           ( 1._wp + r3v(ji,jj,Kaa) ) * vmask(ji,jj,jk)
#  else
                  uu(ji,jj,jk,Kaa) = (         e3u(ji,jj,jk,Kbb) * uu(ji,jj,jk,Kbb )  &
                     &                 + rDt * e3u(ji,jj,jk,Kmm) * uu(ji,jj,jk,Krhs)  )   &
                     &             /           e3u(ji,jj,jk,Kaa) * umask(ji,jj,jk)
                  vv(ji,jj,jk,Kaa) = (         e3v(ji,jj,jk,Kbb) * vv(ji,jj,jk,Kbb )  &
                     &                 + rDt * e3v(ji,jj,jk,Kmm) * vv(ji,jj,jk,Krhs)  )   &
                     &             /           e3v(ji,jj,jk,Kaa) * vmask(ji,jj,jk)
#endif
               END_3D
            ENDIF
            !
            IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=uu(:,:,:,Kaa), clinfo1='stp stg   - Ua: ', mask1=umask,   &
               &                                  tab3d_2=vv(:,:,:,Kaa), clinfo2=           ' Va: ', mask2=vmask, clinfo3='dyn' )
            !
            !
            !              !---------------!
         CASE ( 3 )        !==  Stage 3  ==!   add left over RHS terms + time stepping
            !              !---------------!
            !
            !                  !*  complete the momentum RHS *!   except ZDF (implicit)
            !
                               CALL dyn_ldf( kstp, Kbb, Kmm, uu, vv, Krhs )     ! lateral mixing    
            !
            IF( ln_zdfosm  )   CALL dyn_osm( kstp,      Kmm, uu, vv, Krhs )     ! OSMOSIS non-local velocity fluxes
            !
            IF( ln_bdy     )   CALL bdy_dyn3d_dmp( kstp, Kbb, uu, vv, Krhs )    ! bdy damping trends
            !
            IF( ln_dyndmp .AND. ln_c1d )  CALL dyn_dmp( kstp, Kbb, Kmm, uu(:,:,:,Krhs), vv(:,:,:,Krhs), Nrhs )   ! internal damping trends- momentum
            !

         END SELECT

# if defined key_agrif
      END DO
      IF( ln_tile ) CALL dom_tile_stop

      IF(.NOT. Agrif_Root() ) THEN                              ! AGRIF:   sponge ==> momentum and tracer RHS
         IF( kstg == 3 )   CALL Agrif_Sponge_dyn
      ENDIF

      ! [TILING OVERLAP] wi is updated in tra_adv_trp on halo points and read in dyn_zdf
      IF( ln_tile .AND. kstg == 3 .AND. ln_zad_Aimp .AND. .NOT. ln_shuman ) CALL dom_tile_copyin( 'wi', wi )
      ! [TILING OVERLAP] zF[uv] are updated in place on halo points in tra_adv_trp
      IF( ln_tile .AND. .NOT. ln_shuman ) CALL dom_tile_copyin( 'zFu', zFu, 'zFv', zFv )
      IF( ln_tile ) CALL dom_tile_start                               ! [tiling] DYN tiling loop (2, continued)
      DO jtile = 1, nijtile
         IF( ln_tile ) CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )
         IF( ln_tile .AND. kstg == 3 .AND. ln_zad_Aimp .AND. .NOT. ln_shuman ) CALL dom_tile_copyin( 'wi', wi )
# endif
         !                     !*  DYN time integration + ZDF  *!   ∆t = rDt
         !
         IF( kstg == 3 )   CALL dyn_zdf( kstp, Kbb, Kmm, Krhs, uu, vv, Kaa  )  ! vertical diffusion and time integration

         !
         !                 !==  All stages: correct the barotropic component ==!   at Kaa = N+1/3, N+1/2 or N+1
         !                                                   
#if ! defined key_PSYCLONE_2p5p0
         ALLOCATE( zub(T2D(0)), zvb(T2D(0)) )
#endif
         !
         DO_2D( 0, 0, 0, 0 )             ! barotropic velocity correction
            zub(ji,jj) = uu_b(ji,jj,Kaa) - SUM( e3u_0(ji,jj,:)*uu(ji,jj,:,Kaa) ) * r1_hu_0(ji,jj)
            zvb(ji,jj) = vv_b(ji,jj,Kaa) - SUM( e3v_0(ji,jj,:)*vv(ji,jj,:,Kaa) ) * r1_hv_0(ji,jj)
         END_2D
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )   ! corrected horizontal velocity
            uu(ji,jj,jk,Kaa) = uu(ji,jj,jk,Kaa) + zub(ji,jj)*umask(ji,jj,jk)
            vv(ji,jj,jk,Kaa) = vv(ji,jj,jk,Kaa) + zvb(ji,jj)*vmask(ji,jj,jk)
         END_3D
         !
#if ! defined key_PSYCLONE_2p5p0
         DEALLOCATE( zub, zvb )
#endif

         !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         !>>>            Tracers : RHS computation + time-stepping            <<<
         !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         !
         IF( .NOT.ln_shuman ) THEN             ! No Shuman averaging- tra_adv_trp can be in tiling loop
            !                                  !      clem: check why we cannot use this statement: IF( .NOT.ln_shuman .OR. ( ln_shuman .AND. kstg /= 3 ) )
            !
#if defined key_agrif
            IF( ln_tile ) CALL dom_tile_copyin( 'zFu', zFu, 'zFv', zFv )
#endif
            !                        !* Update or Compute (VIF) advective transport
            CALL tra_adv_trp( kstp, kstg, nit000, Kbb, Kmm, Kaa, Krhs, zFu, zFv, zFw )
            !
!!gm ==>>> Question: it probably be better to activate BBL at each stages (==> it have to be tested)
            !                                            ! BBL coefficients required for both passive- and active-tracer transport within
            !                                            ! the BBL (stage 3 only, requires uu, vv, gdept at Kmm)
            IF( kstg == 3 .AND. ln_trabbl )   CALL bbl( kstp, nit000, Kbb, Kmm )
            !
#if defined key_agrif
            IF( ln_tile ) CALL dom_tile_copyout( 'zFu', zFu, 'zFv', zFv )
            IF( ln_tile .AND. kstg == 3 .AND. ln_zad_Aimp ) CALL dom_tile_copyout( 'wi', wi )
#endif
         ENDIF
         !
      END DO
      !
      IF( ln_tile ) CALL dom_tile_stop
#if defined key_agrif
      IF( ln_tile .AND. .NOT. ln_shuman ) CALL dom_tile_copyout( 'zFu', zFu, 'zFv', zFv )
      IF( ln_tile .AND. kstg == 3 .AND. ln_zad_Aimp .AND. .NOT. ln_shuman ) CALL dom_tile_copyout( 'wi', wi )
#endif
      !
      IF( ln_shuman ) THEN                     ! Shuman averaging- tra_adv_trp not in tiling loop due to lbc_lnk
         !                                     !   clem: check why we cannot use this statement: IF( ln_shuman .AND. kstg == 3 )
         CALL lbc_lnk( 'stp_stg', uu(:,:,:,Kaa), 'U', -1._wp, vv(:,:,:,Kaa), 'V', -1._wp, ldfull=.TRUE. )
         !
         IF( ln_bdy ) THEN                               !* BDY open boundaries
            IF( ln_dynspg_exp )   CALL bdy_dyn( kstp, Kbb, uu, vv, Kaa )
            IF( ln_dynspg_ts  )   CALL bdy_dyn( kstp, Kbb, uu, vv, Kaa, dyn3d_only=.true. )
         ENDIF
         !
         !                           !* Update or Compute (VIF) advective transport
         CALL tra_adv_trp( kstp, kstg, nit000, Kbb, Kmm, Kaa, Krhs, zFu, zFv, zFw )
         !
         !                                            ! BBL coefficients required for both passive- and active-tracer transport within
         !                                            ! the BBL (stage 3 only, requires uu, vv, gdept at Kmm)
         IF( kstg == 3 .AND. ln_trabbl )   CALL bbl( kstp, nit000, Kbb, Kmm )
         !
      ENDIF
      !
# if defined key_top
      !                       !==  Passive Tracer  ==!
      IF( ln_top )   CALL trc_stp( kstg, kstp, Kbb, Kmm, Krhs, Kaa, zFu, zFv, zFw )
      !
# endif
      !
      !                       !==  T-S Tracers  ==!
      !
!!gm ===>>> Modify tra_adv_...  routines so that Krhs to zero useless N.B. issue with passive tracer (trc_sbc and _sms first !)
      DO jn = 1, jpts
         ts(:,:,:,jn,Krhs) = 0._wp                                                ! set tracer trends to zero (:,:,:) needed otherwise it does not work (?)
      END DO
      !
      IF( ln_tile )   CALL dom_tile_start         ! [tiling] TRA tiling loop (1)
      DO jtile = 1, nijtile
         IF( ln_tile )   CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )
         !
         CALL tra_adv    ( kstp, Kbb, Kmm, Kaa, ts, Krhs, zFu, zFv, zFw, kstg )   ! horizontal & vertical advection
         !
         CALL tra_sbc( kstp, Kbb, Kmm,      ts, Krhs,                kstg )   ! surface boundary condition

!!gm ===>>> Question: I'm surprised not to see kstg in argument of tra_isf : potential BUG.to be checked..
         IF( ln_isf )   CALL tra_isf( kstp, Kmm, ts, Krhs )                       ! ice shelf heat flux
      
      END DO
      IF( ln_tile ) CALL dom_tile_stop
      !
!!mjb assimilation increments are applied at all RK3 stages because tests applying them only at stage 3
!!mjb gave strangely poor results. The additional computational expense is minor.  
      IF(  lk_asminc .AND. ln_asmiau ) THEN               ! apply assimilation increment
         IF( ln_trainc )   CALL tra_asm_inc( kstp, Kbb, Kmm, ts    , Krhs )   ! tracers    ==> RHS
      ENDIF

      !   !=============================================================!
      !   !==   stage 1 & 2 : time-stepping                           ==!
      !   !==   stage 3     : time-stepping with left over RHS trends ==!
      !   !=============================================================!
      !
      SELECT CASE( kstg )
      !                 !-------------------!
      CASE ( 1 , 2 )    !==  Stage 1 & 2  ==!    time stepping
         !              !-------------------!
         !
         DO jn = 1, jpts
            DO_3D( 0, 0, 0, 0, 1, jpkm1 )
#  if defined key_linssh
               !              ! linear ssh : applied on tracer
               ts(ji,jj,jk,jn,Kaa) = ( ts(ji,jj,jk,jn,Kbb ) + rDt * ts(ji,jj,jk,jn,Krhs) &
                                                                  / e3t(ji,jj,jk,Krhs)   * tmask(ji,jj,jk) )   !!smg
!!smg          ts(ji,jj,jk,jn,Kaa) = ( ts(ji,jj,jk,jn,Kbb ) + rDt * ts(ji,jj,jk,jn,Krhs) * tmask(ji,jj,jk) )
#  elif defined key_qco  
               !              ! qco : thickness weighted time-stepping using (1+r3.) only
!!smg               ts(ji,jj,jk,jn,Kaa) = (        ( 1._wp + r3t(ji,jj,Kbb) )*ts(ji,jj,jk,jn,Kbb )                      &
!!smg                  &                   + rDt * ( 1._wp + r3t(ji,jj,Kmm) )*ts(ji,jj,jk,jn,Krhs)*tmask(ji,jj,jk)  )   &
!!smg                  &                /          ( 1._wp + r3t(ji,jj,Kaa) )
               ts(ji,jj,jk,jn,Kaa) = (        ( 1._wp + r3t(ji,jj,Kbb) )*ts(ji,jj,jk,jn,Kbb )                                      &  
                  &                   + rDt                             *ts(ji,jj,jk,jn,Krhs)/e3t_0(ji,jj,jk)*tmask(ji,jj,jk)  )   &  !!smg 
                  &                /          ( 1._wp + r3t(ji,jj,Kaa) )
#  else   
               !              ! 4D e3t : thickness weighted time-stepping using e3t
!!smg               ts(ji,jj,jk,jn,Kaa) = (         e3t(ji,jj,jk,Kbb)*ts(ji,jj,jk,jn,Kbb )                      &
!!smg                  &                    + rDt * e3t(ji,jj,jk,Kmm)*ts(ji,jj,jk,jn,Krhs)*tmask(ji,jj,jk)  )   &
!!smg                  &                /           e3t(ji,jj,jk,Kaa)
               ts(ji,jj,jk,jn,Kaa) = (         e3t(ji,jj,jk,Kbb)*ts(ji,jj,jk,jn,Kbb )                      &
                  &                    + rDt *                   ts(ji,jj,jk,jn,Krhs)*tmask(ji,jj,jk)  )   &  !!smg
                  &                /           e3t(ji,jj,jk,Kaa)
#  endif
            END_3D
         END DO
         !
         IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=ts(:,:,:,jp_tem,Kaa), clinfo1='stp stg   - Ta: ', mask1=tmask,   &
            &                                  tab3d_2=ts(:,:,:,jp_sal,Kaa), clinfo2=           ' Sa: ', mask2=tmask, clinfo3='tra' )
         !
         !              !---------------!
      CASE ( 3 )        !==  Stage 3  ==!   add all left over RHS terms
         !              !---------------!
         !
         !
         IF( ln_bdy     )   CALL bdy_tra_dmp   ( kstp, Kbb, ts    , Krhs )   ! bdy damping trends     ==> RHS
         !
# if defined key_agrif
         IF(.NOT. Agrif_Root() ) THEN                        ! AGRIF:   sponge ==> momentum and tracer RHS
            CALL Agrif_Sponge_tra
         ENDIF
# endif
         IF( ln_tile )   CALL dom_tile_start         ! [tiling] TRA tiling loop (2)
         DO jtile = 1, nijtile
            IF( ln_tile )   CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )
            !
            !           !==  complete the tracers RHS  ==!   except ZDF (implicit)
            !
            IF( ln_traqsr  )   CALL tra_qsr( kstp,      Kmm, ts, Krhs )  ! penetrative solar radiation qsr
                               CALL tra_ldf( kstp, Kbb, Kmm, ts, Krhs )  ! lateral mixing
            IF( ln_trabbc  )   CALL tra_bbc( kstp,      Kmm, ts, Krhs )  ! bottom heat flux
            IF( ln_trabbl  )   CALL tra_bbl( kstp, Kbb, Kmm, ts, Krhs )  ! advective (and/or diffusive) bottom boundary layer scheme
            IF( ln_tradmp  )   CALL tra_dmp( kstp, Kbb, Kmm, ts, Krhs )  ! internal damping trends

            IF( ln_zdfmfc  )   CALL tra_mfc( kstp, Kbb,      ts, Krhs )  ! Mass Flux Convection
            IF( ln_zdfosm  ) THEN
                               CALL tra_osm( kstp,      Kmm, ts, Krhs )  ! OSMOSIS non-local tracer fluxes ==> RHS
            ENDIF
            !
            !           !== TRA time integration + ZDF  ==!   
            !
                               CALL tra_zdf( kstp, Kbb, Kmm, Krhs, ts    , Kaa  )  ! vertical mixing and after tracer fields
            IF( ln_zdfnpc  )   CALL tra_npc( kstp,      Kmm, Krhs, ts    , Kaa  )  ! update after fields by non-penetrative convection
            !
         END DO
         IF( ln_tile ) CALL dom_tile_stop

         !!smg 
         IF( lo_trd_tot )   CALL trd_tra_tot( Kbb, Kaa )
         
         IF( ln_zdfosm .AND. lrst_oce ) CALL osm_rst( kstp, Kmm, 'WRITE' )   ! Write OSMOSIS fields and ww to restart file
         !
         IF( .NOT.lk_linssh ) THEN
            r3f(:,:) = r3fa(:,:)                                         ! save r3fa in r3f before deallocation
            DEALLOCATE( r3fa )                                           ! (r3f = r3f(Kbb) of the next time step)
         ENDIF
         !
      END SELECT

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Set boundary conditions
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
!!gm ===>>>  I think the boundary cond. should be done just after the time-stepping of DYN and TRA, not here !
!!           this will probably remove the Shuman case since the is a already a lbc_lnk in shuman casea
# if defined key_agrif
            CALL Agrif_tra( kstp, kstg )   !* AGRIF zoom boundaries
            CALL Agrif_dyn( kstp, kstg )
# endif
      !                                              !* local domain boundaries
      IF( ln_shuman ) THEN
         CALL lbc_lnk( 'stp_RK3_stg', ts(:,:,:,jp_tem,Kaa), 'T',  1._wp, ts(:,:,:,jp_sal,Kaa), 'T',  1._wp )
      ELSE
         CALL lbc_lnk( 'stp_RK3_stg', uu(:,:,:,       Kaa), 'U', -1._wp, vv(:,:,:       ,Kaa), 'V', -1._wp   &
            &                       , ts(:,:,:,jp_tem,Kaa), 'T',  1._wp, ts(:,:,:,jp_sal,Kaa), 'T',  1._wp )
      ENDIF
      !                                              !* BDY open boundaries
      IF( ln_bdy )   THEN
                                  CALL bdy_tra( kstp, Kbb, ts,     Kaa )
         IF( .NOT.ln_shuman ) THEN
            IF( ln_dynspg_exp )   CALL bdy_dyn( kstp, Kbb, uu, vv, Kaa )
            IF( ln_dynspg_ts  )   CALL bdy_dyn( kstp, Kbb, uu, vv, Kaa, dyn3d_only=.TRUE. )
         ENDIF
      ENDIF
      !
#if ! defined key_PSYCLONE_2p5p0
      DEALLOCATE( zFu, zFv, zFw )
#endif
      !
      IF( ln_timing )   CALL timing_stop('stp_stg')
      !
   END SUBROUTINE stp_stg

#else
   !!----------------------------------------------------------------------
   !!   default option             EMPTY MODULE           qco not activated
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE stprk3_stg
