MODULE stp2d
   !!======================================================================
   !!                       ***  MODULE stp2d  ***
   !! Time-stepping   : manager of the ocean, tracer and ice time stepping
   !!                   using a 3rd order Rung Kuta  with fixed or quasi-eulerian coordinate
   !!======================================================================
   !! History :  5.0  !  2021-01  (S. Techene, G. Madec, N. Ducousso, F. Lemarie)  Original code
   !!   NEMO     5.0  !  2024-01  (S. Techene, G. Madec)  remove duplicated calcul of 3D RHS in stp_2D and 1st stage of RK3
   !!----------------------------------------------------------------------
#if defined key_qco   ||   defined key_linssh
   !!----------------------------------------------------------------------
   !!   'key_qco'                        Quasi-Eulerian vertical coordinate
   !!                          OR
   !!   'key_linssh                       Fixed in time vertical coordinate
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!    stp_2D       : RK3 case
   !!----------------------------------------------------------------------
   USE step_oce       ! time stepping used modules
   USE domqco         ! quasi-eulerian coordinate     (dom_qco_r3c  routine)
   USE dynspg_ts      ! 2D mode integration
   USE dynadv_cen2    ! centred flux form advection   (dyn_adv_cen2 routine)
   USE dynadv_up3     ! UP3 flux form advection       (dyn_adv_up3  routine)
   USE dynkeg         ! kinetic energy gradient       (dyn_keg      routine)
   USE dynzad         ! vertical advection            (dyn_zad      routine)
   USE sbc_ice , ONLY : snwice_mass, snwice_mass_b
   USE sbcapr         ! surface boundary condition: atmospheric pressure
   USE sbcwave,  ONLY : bhd_wave
#if defined key_agrif
   USE agrif_oce_interp
   USE agrif_oce_sponge
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp_2D    ! called by stprk3.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE stp_2D( kt, Kbb, Kmm, Kaa, Krhs )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp_2D  ***
      !!
      !! ** Purpose : - Compute sea-surface height and barotropic velocity at Kaa
      !!                in single 1st RK3.
      !!
      !! ** Method  : -1- Compute the 3D to 2D Right Hand Side
      !!                 * Momentum (uu,vv)_rhs is 3D :
      !!                    -- HPG Hor. pressure gradient
      !!                    -- LDF Lateral mixing 
      !!                    -- VOR (COR+RVO) (Vector Inv. Form) (ViF)
      !!                           (COR+MET) (Flux Form) (FF)
      !!                    -- KEG + ZAD     (ViF) 
      !!                       ADV            (FF)
      !!                    In ViF it is also the 1st stage RHS (not recomputed in stprk3_stg)
      !!                    In FF  only ADV has to be recomputed at the 1st stage 
      !!
      !!                 * (Ue_rhs, Ve_rhs) 3D to 2D dynamics : vertical averaging of 3D RHS
      !!                    Note that in FF case, ADV return the vertical averaging and does not 
      !!                    update the 3D RHS, so that ADV can be added to the 3D RHS at 1st stage of RK3 
      !!
      !!                 * External forcings added to 2D RHS:
      !!                    -- baroclinic drag
      !!                    -- wind 
      !!                    -- atmospheric pressure
      !!                    -- snow+ice load
      !!                    -- surface wave load
      !!                 * ssh (sshe_rhs) :
      !!                      Net column average freshwater flux
      !!
      !!              -2- Solve the external mode Eqs. using sub-time stepping
      !!                  by a call to dyn_spg_ts (will be renamed stp_2D_solver)
      !!
      !! ** action  :   ssh         at Kaa : N+1 sea surface height 
      !!                (uu_b,vv_b) at Kaa : N+1 barotropic velocity 
      !!                (un_adv,vn_adv): barotropic transport from N to N+1  [m2/s]
      !!                (uu,vv)_rhs : 1st stage RHS in ViF, except ADV in FF
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, Kbb, Kmm, Kaa, Krhs   ! ocean time-step and time-level indices
      !
      INTEGER  ::   ji, jj, jk, jtile                    ! dummy loop indices
      REAL(wp) ::   zg_2, zintp, zgrho0r, zld, zztmp     ! local scalars
#if ! defined key_PSYCLONE_2p5p0
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zpice   ! 2D workspace
#else
      REAL(wp), DIMENSION(T2D(1)) ::   zpice             ! 2D workspace
#endif
      !! ---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('stp_2D')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'stp_2D : barotropic field in single first '
         IF(lwp) WRITE(numout,*) '~~~~~~'
      ENDIF
      !
      IF( lk_linssh ) THEN    !==  Compute ww(:,:,1)  ==!   (needed for momentum advection)
!!gm  only in Flux Form, Vector Form  dzU_z=0 assumed to be zero 
!!gm  ww(k=1) = div_h(uu_b) ==> modif dans dynadv                        <<<=== TO BE DONE
      ENDIF

      ALLOCATE( sshe_rhs(jpi,jpj) , Ue_rhs(A2D(0)) , Ve_rhs(A2D(0)) , CdU_u(jpi,jpj) , CdU_v(jpi,jpj) )
      IF( nn_bt_frc > 1 )   ALLOCATE( Ue_xtp(A2D(0)) , Ve_xtp(A2D(0)) )   ! advective part of the 2D RHS (extrapolated in dynspg_ts)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
      !                      RHS of barotropic momentum  Eq.                  !
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!

      !                 !======================================!   Flux Form        : HPG (+ LDF) + COR + MET + ADV 
      !                 !==  Dynamics 2D RHS from 3D trends  ==!
      !                 !======================================!   Vector Inv. Form : HPG (+ LDF) + VOR + RVO + KEG + ZAD
      !
      IF( ln_tile ) CALL dom_tile_start         ! [tiling] DYN tiling loop (1)
         DO jtile = 1, nijtile
         IF( ln_tile ) CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )

         !                             !*  hydrostatic pressure gradient (HPG))  *!   always called FIRST
         CALL eos    ( ts, Kbb, rhd )                          ! in situ density anomaly at Kbb
         CALL dyn_hpg( kt, Kbb     , uu, vv, Krhs )            ! horizontal gradient of Hydrostatic pressure
         !
         !                             !*  lateral viscosity (LDF)  *!
         CALL dyn_ldf( kt, Kbb, Kbb, uu, vv, Krhs )
#if defined key_agrif
      END DO
      IF( ln_tile ) CALL dom_tile_stop

      IF(.NOT. Agrif_Root() ) THEN     !*  AGRIF: sponge  *!
         CALL Agrif_Sponge_dyn
      ENDIF

      IF( ln_tile ) CALL dom_tile_start         ! [tiling] DYN tiling loop (1, continued)
      DO jtile = 1, nijtile
         IF( ln_tile ) CALL dom_tile( ntsi, ntsj, ntei, ntej, ktile = jtile )
#endif
         !                             !*  COR + MET  *!   Flux Form        : Coriolis + Metric Term  
         !                             !*     VOR     *!   Vector Inv. Form : Coriolis + relative Vorticity
         IF( nn_bt_frc > 1 .AND. n_dynadv /= np_LIN_dyn ) THEN   ! planetary vorticity (+ Stokes-Coriolis) first, then save
            CALL dyn_vor( kt,      Kbb, uu, vv, Krhs, ncor )      ! the depth mean of the 2D RHS before its advective part
            DO_2D( 0, 0, 0, 0 )
               Ue_xtp(ji,jj) = - SUM( e3u_0(ji,jj,1:jpkm1)*uu(ji,jj,1:jpkm1,Krhs)*umask(ji,jj,1:jpkm1) ) * r1_hu_0(ji,jj)
               Ve_xtp(ji,jj) = - SUM( e3v_0(ji,jj,1:jpkm1)*vv(ji,jj,1:jpkm1,Krhs)*vmask(ji,jj,1:jpkm1) ) * r1_hv_0(ji,jj)
            END_2D
            CALL dyn_vor( kt,      Kbb, uu, vv, Krhs, nrvm )      ! relative vorticity (VIF) or metric term (FF)
         ELSE
            CALL dyn_vor( kt,      Kbb, uu, vv, Krhs )
            IF( nn_bt_frc > 1 ) THEN                              ! linear dynamics: no advective forcing
               DO_2D( 0, 0, 0, 0 )
                  Ue_xtp(ji,jj) = - SUM( e3u_0(ji,jj,1:jpkm1)*uu(ji,jj,1:jpkm1,Krhs)*umask(ji,jj,1:jpkm1) ) * r1_hu_0(ji,jj)
                  Ve_xtp(ji,jj) = - SUM( e3v_0(ji,jj,1:jpkm1)*vv(ji,jj,1:jpkm1,Krhs)*vmask(ji,jj,1:jpkm1) ) * r1_hv_0(ji,jj)
               END_2D
            ENDIF
         ENDIF
         !
         !                             !*  compute advection  *!   (only KEG + ZAD in Vector Inv. Form)
         IF( .NOT.lk_linssh ) THEN
            DO_2D( 1, nn_hls, 1, nn_hls )      ! loop bounds limited by ssh definition in ssh_nxt
               r3t(ji,jj,Kaa) =  ssh(ji,jj,Kaa) * r1_ht_0(ji,jj)               ! "after" ssh/h_0 ratio guess at t-column at Kaa (n+1)
            END_2D
         ENDIF
         !
         CALL wzv    ( kt, Kbb, Kbb, Kaa , uu(:,:,:,Kbb), vv(:,:,:,Kbb), ww, np_velocity )  ! ww guess at Kbb (n)
         !
         !                             !*  KEG + ZAD  *!   Vector Inv. Form : KE gradient + vertical advection
         !                             !*     ADV     *!   Flux Form        : flux form advection
         SELECT CASE( n_dynadv )
         !
         CASE( np_VEC_c2  )               != 2nd order Vector Form =!                       ==>> 3D RHS
            !
            CALL dyn_keg( kt, nn_dynkeg, Kbb, uu, vv, Krhs )      !- horizontal Gradient of KE
            !
            CALL dyn_zad( kt, Kbb, uu, vv, Krhs )                 !- vertical advection
            !
         CASE( np_FLX_c2  )               !=  2n order Flux Form =!   (CEN2)                ==>> 2D RHS only
            !
            CALL dyn_adv_cen2( kt     , Kbb, uu, vv, Krhs, pUe=Ue_rhs, pVe=Ve_rhs )
            !
         CASE( np_FLX_up3 )               != 3rd order Flux Form =!   (UP3)                 ==>> 2D RHS only
            CALL dyn_adv_up3 ( kt, Kbb, Kbb, uu, vv, Krhs, pUe=Ue_rhs, pVe=Ve_rhs )
            !
         END SELECT
         !
         !                             !*  vertical averaging  *!
         SELECT CASE( n_dynadv )
         CASE( np_VEC_c2, np_LIN_dyn )       ! Vector Inv. Form   ==>> averaged 3D RHS only
            DO_2D( 0, 0, 0, 0 )
               Ue_rhs(ji,jj) = SUM( e3u_0(ji,jj,1:jpkm1)*uu(ji,jj,1:jpkm1,Krhs)*umask(ji,jj,1:jpkm1) ) * r1_hu_0(ji,jj)
               Ve_rhs(ji,jj) = SUM( e3v_0(ji,jj,1:jpkm1)*vv(ji,jj,1:jpkm1,Krhs)*vmask(ji,jj,1:jpkm1) ) * r1_hv_0(ji,jj)
            END_2D
         CASE ( np_FLX_c2, np_FLX_up3 )      ! Flux Form          ==>> cumulated ADV 2D RHS with 3D RHS
            DO_2D( 0, 0, 0, 0 )
               Ue_rhs(ji,jj) = Ue_rhs(ji,jj) + SUM( e3u_0(ji,jj,1:jpkm1)*uu(ji,jj,1:jpkm1,Krhs)*umask(ji,jj,1:jpkm1) ) * r1_hu_0(ji,jj)
               Ve_rhs(ji,jj) = Ve_rhs(ji,jj) + SUM( e3v_0(ji,jj,1:jpkm1)*vv(ji,jj,1:jpkm1,Krhs)*vmask(ji,jj,1:jpkm1) ) * r1_hv_0(ji,jj)
            END_2D
         END SELECT
         !
         IF( nn_bt_frc > 1 ) THEN          !* advective part of the 2D RHS (extrapolated in time in dynspg_ts)
            DO_2D( 0, 0, 0, 0 )             !  VIF : RVO + KEG + ZAD   ;   FF : MET + ADV
               Ue_xtp(ji,jj) = Ue_xtp(ji,jj) + Ue_rhs(ji,jj)
               Ve_xtp(ji,jj) = Ve_xtp(ji,jj) + Ve_rhs(ji,jj)
            END_2D
         ENDIF

         !              !=====================================!
         !              !==  Dynamics: 2D momentum forcing  ==!
         !              !=====================================!
         !
         !
         !                             !* baroclinic drag forcing *!   (also provide the barotropic drag coeff.)
         CALL dyn_drg_init( Kbb, uu, vv, uu_b, vv_b, Ue_rhs, Ve_rhs, CdU_u, CdU_v )
         !
         !                             !* wind forcing *!
         DO_2D( 0, 0, 0, 0 )
            Ue_rhs(ji,jj) =  Ue_rhs(ji,jj) + r1_rho0 * utauU(ji,jj) * r1_hu(ji,jj,Kbb)
            Ve_rhs(ji,jj) =  Ve_rhs(ji,jj) + r1_rho0 * vtauV(ji,jj) * r1_hv(ji,jj,Kbb)
         END_2D
         !
         !                             !* atmospheric pressure forcing *!
         IF( ln_apr_dyn ) THEN
            DO_2D( 0, 0, 0, 0 )
               Ue_rhs(ji,jj) = Ue_rhs(ji,jj) + grav * (  ssh_ib (ji+1,jj  ) - ssh_ib (ji,jj) ) * r1_e1u(ji,jj)
               Ve_rhs(ji,jj) = Ve_rhs(ji,jj) + grav * (  ssh_ib (ji  ,jj+1) - ssh_ib (ji,jj) ) * r1_e2v(ji,jj)
            END_2D
         ENDIF
         !
         !                             !* snow+ice load *!   (embedded sea ice)
         IF( ln_ice_embd ) THEN
#if ! defined key_PSYCLONE_2p5p0
            ALLOCATE( zpice(T2D(1)) )
#endif
            zintp = REAL( MOD( kt-1, nn_fsbc ) ) / REAL( nn_fsbc )
            zgrho0r = - grav * r1_rho0
            DO_2D( 1, 1, 1, 1 )
               zpice(ji,jj) = ( zintp * snwice_mass(ji,jj) + (1._wp - zintp) * snwice_mass_b(ji,jj) ) * zgrho0r
            END_2D
            DO_2D( 0, 0, 0, 0 )
               Ue_rhs(ji,jj) = Ue_rhs(ji,jj) + ( zpice(ji+1,jj) - zpice(ji,jj) ) * r1_e1u(ji,jj)
               Ve_rhs(ji,jj) = Ve_rhs(ji,jj) + ( zpice(ji,jj+1) - zpice(ji,jj) ) * r1_e2v(ji,jj)
            END_2D
#if ! defined key_PSYCLONE_2p5p0
            DEALLOCATE( zpice )
#endif
         ENDIF
         !
         !                             !* surface wave load *!   (Bernoulli head)
         !
         IF( ln_bern_srfc ) THEN
            DO_2D( 0, 0, 0, 0 )
               Ue_rhs(ji,jj) = Ue_rhs(ji,jj) + ( bhd_wave(ji+1,jj) - bhd_wave(ji,jj) ) * r1_e1u(ji,jj)   !++ bhd_wave from wave model in m2/s2 [BHD parameters in WW3]
               Ve_rhs(ji,jj) = Ve_rhs(ji,jj) + ( bhd_wave(ji,jj+1) - bhd_wave(ji,jj) ) * r1_e2v(ji,jj)
            END_2D
         ENDIF
      END DO

      IF( ln_tile ) CALL dom_tile_stop

      !                 !=======================================!
      !                 !==   2D sea surface height forcing   ==!
      !                 !=======================================!
      !
      !                             !*  Net water flux forcing   (applied to a water column)
                     sshe_rhs(:,:) =                 emp(:,:)
      IF( ln_rnf )   sshe_rhs(:,:) = sshe_rhs(:,:) - rnf(:,:)
      IF( ln_isf )   sshe_rhs(:,:) = sshe_rhs(:,:) - fwfisf_cav(:,:) - fwfisf_par(:,:)
                     sshe_rhs(:,:) = r1_rho0 * sshe_rhs(:,:)
      !
      !                             !* Stokes drift divergence
      IF( ln_sdw )   sshe_rhs(:,:) = sshe_rhs(:,:) + div_sd(:,:)
      !
      !
      !                             !* ice sheet coupling
      IF( ln_isf .AND. ln_isfcpl ) THEN
         IF( ln_rstart .AND. kt == nit000 )   sshe_rhs(:,:) = sshe_rhs(:,:) + risfcpl_ssh(:,:)
         IF( ln_isfcpl_cons               )   sshe_rhs(:,:) = sshe_rhs(:,:) + risfcpl_cons_ssh(:,:)
      ENDIF
      !
#if defined key_asminc
      !                             !* Add the IAU weighted SSH increment
      IF( lk_asminc .AND. ln_sshinc .AND. ln_asmiau )   sshe_rhs(:,:) = sshe_rhs(:,:) - ssh_iau(:,:)
#endif
      !
#if defined key_agrif
      !                             !* AGRIF : fill boundary data arrays (on both )
      IF( .NOT.Agrif_Root() )   CALL agrif_dta_ts( kt )
#endif

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !             Compute ssh and (uu_b,vv_b)  at N+1  (Kaa)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      !    using a split-explicit time integration in forward mode 
      !    ( ABM3-AM4 time-integration Shchepetkin et al. OM2005) with temporal diffusion (Demange et al. JCP2019) )

      IF( ln_dynspg_ts )   &
         &   CALL dyn_spg_ts( kt, Kbb, Kbb, ssh, uu_b, vv_b, Kaa ) ! time-splitting
                   

      DEALLOCATE( sshe_rhs , Ue_rhs , Ve_rhs , CdU_u , CdU_v )
      IF( nn_bt_frc > 1 )   DEALLOCATE( Ue_xtp , Ve_xtp )
      !
      IF( ln_timing )   CALL timing_stop('stp_2D')
      !
   END SUBROUTINE stp_2D


#else
   !!----------------------------------------------------------------------
   !!   default option             EMPTY MODULE           qco not activated
   !!----------------------------------------------------------------------
#endif
   
   !!======================================================================
END MODULE stp2d
