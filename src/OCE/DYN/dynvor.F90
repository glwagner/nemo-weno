MODULE dynvor
   !!======================================================================
   !!                       ***  MODULE  dynvor  ***
   !! Ocean dynamics: Update the momentum trend with the relative and
   !!                 planetary vorticity trends
   !!======================================================================
   !! History :  OPA  ! 1989-12  (P. Andrich)  vor_ens: Original code
   !!            5.0  ! 1991-11  (G. Madec)  vor_ene, vor_mix: Original code
   !!            6.0  ! 1996-01  (G. Madec)  s-coord, suppress work arrays
   !!   NEMO     0.5  ! 2002-08  (G. Madec)  F90: Free form and module
   !!            1.0  ! 2004-02  (G. Madec)  vor_een: Original code
   !!             -   ! 2003-08  (G. Madec)  add vor_ctl
   !!             -   ! 2005-11  (G. Madec)  add dyn_vor (new step architecture)
   !!            2.0  ! 2006-11  (G. Madec)  flux form advection: add metric term
   !!            3.2  ! 2009-04  (R. Benshila)  vvl: correction of een scheme
   !!            3.3  ! 2010-10  (C. Ethe, G. Madec)  reorganisation of initialisation phase
   !!            3.7  ! 2014-04  (G. Madec)  trend simplification: suppress jpdyn_trd_dat vorticity
   !!             -   ! 2014-06  (G. Madec)  suppression of velocity curl from in-core memory
   !!             -   ! 2016-12  (G. Madec, E. Clementi) add Stokes-Coriolis trends 
   !!            4.0  ! 2017-07  (G. Madec)  linear dynamics + trends diag. with Stokes-Coriolis
   !!             -   ! 2018-03  (G. Madec)  add two new schemes (ln_dynvor_enT and ln_dynvor_eet)
   !!             -   ! 2018-04  (G. Madec)  add pre-computed gradient for metric term calculation
   !!            4.x  ! 2020-03  (G. Madec, A. Nasser)  make ln_dynvor_msk truly efficient on relative vorticity
   !!            4.2  ! 2020-12  (G. Madec, E. Clementi) add vortex force trends (ln_vortex_force=T)
   !!            4.5  ! 2022-06  (S. Techene, G, Madec) refactorization to reduce local memory usage
   !!            5.0  ! 2025-06  (J, Chanut) suppress mixed energy/enstrophy scheme (ln_dynvor_mix) 
   !!            5.x  ! 2026-10  (G. Wagner) WENO upwind relative vorticity (ln_dynvor_weno), following Oceananigans.jl
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_vor       : Update the momentum trend with the vorticity trend
   !!       vor_enT   : energy conserving scheme at T-pt  (ln_dynvor_enT=T)
   !!       vor_ene   : energy conserving scheme          (ln_dynvor_ene=T)
   !!       vor_ens   : enstrophy conserving scheme       (ln_dynvor_ens=T)
   !!       vor_een   : energy and enstrophy conserving   (ln_dynvor_een=T)
   !!       vor_weno  : WENO upwind relative vorticity    (ln_dynvor_weno=T)
   !!   dyn_vor_init  : set and control of the different vorticity option
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE dommsk         ! ocean mask
   USE dynadv         ! momentum advection
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   USE sbcwave        ! Surface Waves (add Stokes-Coriolis force)
   USE sbc_oce,  ONLY : ln_sdw, ln_vortex_force   ! use Stoke-Coriolis force
   !
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE
   
   INTERFACE dyn_vor
      MODULE PROCEDURE dyn_vor_3D, dyn_vor_RK3
   END INTERFACE

   PUBLIC   dyn_vor        ! routine called by stp2d.F90 and stprk3_stg.F90
   PUBLIC   dyn_vor_init   ! routine called by nemogcm.F90

   !                                   !!* Namelist namdyn_vor: vorticity term
   LOGICAL, PUBLIC ::   ln_dynvor_ens   !: enstrophy conserving scheme          (ENS)
   LOGICAL, PUBLIC ::   ln_dynvor_ene   !: f-point energy conserving scheme     (ENE)
   LOGICAL, PUBLIC ::   ln_dynvor_enT   !: t-point energy conserving scheme     (ENT)
   LOGICAL, PUBLIC ::   ln_dynvor_een   !: energy & enstrophy conserving scheme (EEN)
   LOGICAL, PUBLIC ::   ln_dynvor_msk   !: vorticity multiplied by fmask (=T) or not (=F) (all vorticity schemes)
   INTEGER, PUBLIC ::   nn_e3f_typ      !: e3f=masked averaging of e3t divided by 4 (=0) or by the sum of mask (=1)
   LOGICAL, PUBLIC ::   ln_dynvor_weno  !: relative vorticity reconstructed by WENO (Coriolis from the chosen scheme)
   INTEGER, PUBLIC ::      nn_dynvor_weno  !: order of the WENO reconstruction (=3/5)

   INTEGER, PUBLIC ::   nvor_scheme     !: choice of the type of advection scheme
   !                                       ! associated indices:
   INTEGER, PUBLIC, PARAMETER ::   np_ENS = 0   ! ENS scheme
   INTEGER, PUBLIC, PARAMETER ::   np_ENE = 1   ! ENE scheme
   INTEGER, PUBLIC, PARAMETER ::   np_ENT = 2   ! ENT scheme (t-point vorticity)
   INTEGER, PUBLIC, PARAMETER ::   np_EEN = 3   ! EEN scheme

   !                                    !: choice of calculated vorticity
   INTEGER, PUBLIC ::   ncor, nrvm, ntot   ! Coriolis, relative vorticity, total vorticity
   !                                       ! associated indices:
   INTEGER, PUBLIC, PARAMETER ::   np_COR = 1         ! Coriolis (planetary)
   INTEGER, PUBLIC, PARAMETER ::   np_RVO = 2         ! relative vorticity
   INTEGER, PUBLIC, PARAMETER ::   np_MET = 3         ! metric term
   INTEGER, PUBLIC, PARAMETER ::   np_CRV = 4         ! relative + planetary (total vorticity)
   INTEGER, PUBLIC, PARAMETER ::   np_CME = 5         ! Coriolis + metric term

   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   di_e2u_2        ! = di(e2u)/2          used in T-point metric term calculation
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   dj_e1v_2        ! = dj(e1v)/2           -        -      -       -
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   di_e2v_2e1e2f   ! = di(e2u)/(2*e1e2f)  used in F-point metric term calculation
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   dj_e1u_2e1e2f   ! = dj(e1v)/(2*e1e2f)   -        -      -       -
   !
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:), PUBLIC ::   e3f_0vor   ! e3f used in EEN, ENE and ENS cases (key_qco only)

   REAL(wp) ::   r1_4  = 0.250_wp         ! =1/4
   REAL(wp) ::   r1_8  = 0.125_wp         ! =1/8
   REAL(wp) ::   r1_12 = 1._wp / 12._wp   ! 1/12

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE dyn_vor_3D( kt, Kmm, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose :   compute the lateral ocean tracer physics.
      !!
      !! ** Action : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the now vorticity term trend
      !!             - save the trends in (ztrdu,ztrdv) in 2 parts (relative
      !!               and planetary vorticity trends) and send them to trd_dyn
      !!               for futher diagnostics (l_trddyn=T)
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in  ) ::   kt          ! ocean time-step index
      INTEGER                             , INTENT( in  ) ::   Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::   puu, pvv    ! ocean velocity field and RHS of momentum equation
      !
      INTEGER ::   ivor   ! vorticity computed by the chosen scheme
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  ztrdu, ztrdv
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_vor_3D')
      !
      IF( ln_dynvor_weno ) THEN   ;   ivor = ncor      ! Coriolis only, relative vorticity from vor_weno
      ELSE                        ;   ivor = ntot
      ENDIF
      !
      IF( l_trddyn ) THEN     !==  trend diagnostics case : split the added trend in two parts  ==!
         !
         ALLOCATE( ztrdu(T2D(0),jpk), ztrdv(T2D(0),jpk) )
         !
         ztrdu(:,:,:) = puu(T2D(0),:,Krhs)            !* planetary vorticity trend
         ztrdv(:,:,:) = pvv(T2D(0),:,Krhs)
         SELECT CASE( nvor_scheme )
         CASE( np_ENS )           ;   CALL vor_ens( kt, Kmm, ncor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! enstrophy conserving scheme
         CASE( np_ENE )           ;   CALL vor_ene( kt, Kmm, ncor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! energy conserving scheme
         CASE( np_ENT )           ;   CALL vor_enT( kt, Kmm, ncor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! energy conserving scheme (T-pts)
         CASE( np_EEN )           ;   CALL vor_een( kt, Kmm, ncor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! energy & enstrophy scheme
         END SELECT
         ztrdu(:,:,:) = puu(T2D(0),:,Krhs) - ztrdu(:,:,:)
         ztrdv(:,:,:) = pvv(T2D(0),:,Krhs) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_pvo, kt, Kmm )
         !
         IF( n_dynadv /= np_LIN_dyn ) THEN   !* relative vorticity or metric trend (only in non-linear case)
            ztrdu(:,:,:) = puu(T2D(0),:,Krhs)
            ztrdv(:,:,:) = pvv(T2D(0),:,Krhs)
            IF( ln_dynvor_weno ) THEN
                                     CALL vor_weno( kt, Kmm,       puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! WENO upwind scheme
            ELSE
               SELECT CASE( nvor_scheme )
               CASE( np_ENT )        ;   CALL vor_enT( kt, Kmm, nrvm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! energy conserving scheme (T-pts)
               CASE( np_ENE )        ;   CALL vor_ene( kt, Kmm, nrvm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! energy conserving scheme
               CASE( np_ENS )        ;   CALL vor_ens( kt, Kmm, nrvm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! enstrophy conserving scheme
               CASE( np_EEN )        ;   CALL vor_een( kt, Kmm, nrvm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! energy & enstrophy scheme
               END SELECT
            ENDIF
            ztrdu(:,:,:) = puu(T2D(0),:,Krhs) - ztrdu(:,:,:)
            ztrdv(:,:,:) = pvv(T2D(0),:,Krhs) - ztrdv(:,:,:)
            CALL trd_dyn( ztrdu, ztrdv, jpdyn_rvo, kt, Kmm )
         ENDIF
         !
         DEALLOCATE( ztrdu, ztrdv )
         !
      ELSE              !==  total vorticity trend added to the general trend  ==!
         !
         SELECT CASE ( nvor_scheme )      !==  vorticity trend added to the general trend  ==!
         CASE( np_ENT )                        !* energy conserving scheme  (T-pts)
                             CALL vor_enT( kt, Kmm, ivor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! total vorticity trend
            IF( ln_sdw .AND. .NOT. ln_vortex_force )  THEN
                             CALL vor_enT( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
            ELSE IF( ln_sdw .AND. ln_vortex_force )   THEN
                             CALL vor_enT( kt, Kmm, ntot, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend and vortex force
            ENDIF
         CASE( np_ENE )                        !* energy conserving scheme
                             CALL vor_ene( kt, Kmm, ivor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! total vorticity trend
            IF( ln_sdw .AND. .NOT. ln_vortex_force )  THEN
                             CALL vor_ene( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
            ELSE IF( ln_sdw .AND. ln_vortex_force )   THEN
                             CALL vor_ene( kt, Kmm, ntot, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend and vortex force
            ENDIF
         CASE( np_ENS )                        !* enstrophy conserving scheme
                             CALL vor_ens( kt, Kmm, ivor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! total vorticity trend

            IF( ln_sdw .AND. .NOT. ln_vortex_force )  THEN
                             CALL vor_ens( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! add the Stokes-Coriolis trend
            ELSE IF( ln_sdw .AND. ln_vortex_force )   THEN
                             CALL vor_ens( kt, Kmm, ntot, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! add the Stokes-Coriolis trend and vortex force
            ENDIF
         CASE( np_EEN )                        !* energy and enstrophy conserving scheme
                             CALL vor_een( kt, Kmm, ivor, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! total vorticity trend
               IF( ln_sdw .AND. .NOT. ln_vortex_force )  THEN
                             CALL vor_een( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
               ELSE IF( ln_sdw .AND. ln_vortex_force )   THEN
                             CALL vor_een( kt, Kmm, ntot, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend and vortex force
               ENDIF
         END SELECT
         !
         IF( ln_dynvor_weno )   CALL vor_weno( kt, Kmm, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! relative vorticity trend
         !
      ENDIF
      !
      !                       ! print sum trends (used for debugging)
      IF(sn_cfctl%l_prtctl) CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' vor  - Ua: ', mask1=umask,               &
         &                                tab3d_2=pvv(:,:,:,Krhs), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_vor_3D')
      !
   END SUBROUTINE dyn_vor_3D


   SUBROUTINE dyn_vor_RK3( kt, Kmm, puu, pvv, Krhs, knoco )
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose :   compute the lateral ocean tracer physics.
      !!
      !! ** Action : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the now vorticity term trend
      !!             - save the trends in (ztrdu,ztrdv) in 2 parts (relative
      !!               and planetary vorticity trends) and send them to trd_dyn
      !!               for futher diagnostics (l_trddyn=T)
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::   kt          ! ocean time-step index
      INTEGER                             , INTENT(in   ) ::   Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::   puu, pvv    ! ocean velocity field and RHS of momentum equation
      INTEGER                             , INTENT(in   ) ::   knoco       ! specified vorticity trend (= np_MET or np_RVO)
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_vor_RK3')
      !
      !              !==  total vorticity trend added to the general trend  ==!
      !!st WARNING 22/02 !!!!!!!! stoke drift or not stoke drift ? => bar to do later !!!
      !! stoke drift a garder pas vortex force a priori !!
      !! ATTENTION déja appelé avec Kbb !!
      
         !
         SELECT CASE ( nvor_scheme )      !==  vorticity trend added to the general trend  ==!
         CASE( np_ENT )                        !* energy conserving scheme  (T-pts)
                             CALL vor_enT( kt, Kmm, knoco, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! total vorticity trend
            IF( ln_sdw .AND. .NOT. ln_vortex_force )  THEN
                             CALL vor_enT( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
            ELSE IF( ln_sdw .AND. ln_vortex_force )   THEN
                             CALL vor_enT( kt, Kmm, ntot, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend and vortex force
            ENDIF
         CASE( np_ENE )                        !* energy conserving scheme
                             CALL vor_ene( kt, Kmm, knoco, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! total vorticity trend
            IF( ln_sdw .AND. .NOT. ln_vortex_force )  THEN
                             CALL vor_ene( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
            ELSE IF( ln_sdw .AND. ln_vortex_force )   THEN
                             CALL vor_ene( kt, Kmm, ntot, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend and vortex force
            ENDIF
         CASE( np_ENS )                        !* enstrophy conserving scheme
                             CALL vor_ens( kt, Kmm, knoco, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! total vorticity trend

            IF( ln_sdw .AND. .NOT. ln_vortex_force )  THEN
                             CALL vor_ens( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! add the Stokes-Coriolis trend
            ELSE IF( ln_sdw .AND. ln_vortex_force )   THEN
                             CALL vor_ens( kt, Kmm, ntot, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! add the Stokes-Coriolis trend and vortex force
            ENDIF
         CASE( np_EEN )                        !* energy and enstrophy conserving scheme
                             CALL vor_een( kt, Kmm, knoco, puu(:,:,:,Kmm) , pvv(:,:,:,Kmm) , puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! total vorticity trend
            IF( ln_sdw .AND. .NOT. ln_vortex_force )  THEN
                             CALL vor_een( kt, Kmm, ncor, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend
            ELSE IF( ln_sdw .AND. ln_vortex_force )   THEN
                             CALL vor_een( kt, Kmm, ntot, usd, vsd, puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )   ! add the Stokes-Coriolis trend and vortex force
            ENDIF
         END SELECT
      !
      !                       ! print sum trends (used for debugging)
      IF(sn_cfctl%l_prtctl) CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' vor  - Ua: ', mask1=umask,               &
         &                                tab3d_2=pvv(:,:,:,Krhs), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_vor_RK3')
      !
   END SUBROUTINE dyn_vor_RK3


   SUBROUTINE vor_enT( kt, Kmm, kvor, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE vor_enT  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!       and t-point evaluation of vorticity (planetary and relative).
      !!       conserves the horizontal kinetic energy.
      !!         The general trend of momentum is increased due to the vorticity
      !!       term which is given by:
      !!          voru = 1/bu  mj[ ( mi(mj(bf*rvor))+bt*f_t)/e3t  mj[vn] ]
      !!          vorv = 1/bv  mi[ ( mi(mj(bf*rvor))+bt*f_t)/e3f  mj[un] ]
      !!       where rvor is the relative vorticity at f-point
      !!
      !! ** Action : - Update (pu_rhs,pv_rhs) with the now vorticity term trend
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt               ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kmm              ! ocean time level index
      INTEGER                         , INTENT(in   ) ::   kvor             ! total, planetary, relative, or metric
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu, pv           ! now velocities
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs   ! total v-trend
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   zx1, zy1, zx2, zy2   ! local scalars
      REAL(wp), DIMENSION(T2D(1))           ::   zwt   ! 2D workspace
#if ! defined key_PSYCLONE_2p5p0
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zwz   ! 3D workspace, jpkm1 -> avoid lbc_lnk on jpk that is not defined
#else
      REAL(wp), DIMENSION(T2D(1))           ::   zwz   ! 3D workspace, jpkm1 -> avoid lbc_lnk on jpk that is not defined
#endif
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:vor_enT : vorticity term: t-point energy conserving scheme'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         ENDIF
      ENDIF
      !
#if ! defined key_PSYCLONE_2p5p0
      SELECT CASE( kvor )        ! allocate zwz if necessary
      CASE ( np_RVO , np_CRV )   ;   ALLOCATE( zwz(T2D(1)) )
      END SELECT
#endif
      !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         SELECT CASE( kvor )                 !== relative vorticity considered  ==!
            !
         CASE ( np_RVO , np_CRV )                  !* relative vorticity at f-point is used
            DO_2D( 1, 1, 1, 1 )
               zwz(ji,jj) = (  ( e2v(ji+1,jj) * pv(ji+1,jj,jk) - e2v(ji,jj) * pv(ji,jj,jk) )  &                   ! add () for
                  &          - ( e1u(ji,jj+1) * pu(ji,jj+1,jk) - e1u(ji,jj) * pu(ji,jj,jk) )  ) * r1_e1e2f(ji,jj) ! NP repro
            END_2D
            IF( ln_dynvor_msk ) THEN                     ! mask relative vorticity
               DO_2D( 1, 1, 1, 1 )
                  zwz(ji,jj) = zwz(ji,jj) * fmask(ji,jj,jk)
               END_2D
            ENDIF
            !
         END SELECT
         !
         SELECT CASE( kvor )                 !==  volume weighted vorticity considered  ==!
         !
         CASE ( np_COR )                           !* Coriolis (planetary vorticity)
            DO_2D( 0, 1, 0, 1 )
               zwt(ji,jj) = ff_t(ji,jj) * e1e2t(ji,jj)*e3t(ji,jj,jk,Kmm)
            END_2D
         CASE ( np_RVO )                           !* relative vorticity
            DO_2D( 0, 1, 0, 1 )
               zwt(ji,jj) = r1_4 * (   ( zwz(ji-1,jj  ) + zwz(ji,jj-1) )      &   ! need additional () for
                  &                  + ( zwz(ji-1,jj-1) + zwz(ji,jj  ) )  )   &   ! reproducibility around NP
                  &              * e1e2t(ji,jj)*e3t(ji,jj,jk,Kmm)
            END_2D
         CASE ( np_MET )                           !* metric term
            DO_2D( 0, 1, 0, 1 )
               zwt(ji,jj) = (   ( pv(ji,jj,jk) + pv(ji,jj-1,jk) ) * di_e2u_2(ji,jj)       &
                  &           - ( pu(ji,jj,jk) + pu(ji-1,jj,jk) ) * dj_e1v_2(ji,jj)   )   &
                  &       * e3t(ji,jj,jk,Kmm)
            END_2D
         CASE ( np_CRV )                           !* Coriolis + relative vorticity
            DO_2D( 0, 1, 0, 1 )
               zwt(ji,jj) = (  ff_t(ji,jj) + r1_4 * ( ( zwz(ji-1,jj  ) + zwz(ji,jj-1) )        &   ! need additional () for
                  &                                 + ( zwz(ji-1,jj-1) + zwz(ji,jj  ) ) )  )   &   ! reproducibility around NP
                  &       * e1e2t(ji,jj)*e3t(ji,jj,jk,Kmm)
            END_2D
         CASE ( np_CME )                           !* Coriolis + metric
            DO_2D( 0, 1, 0, 1 )
               zwt(ji,jj) = (  ff_t(ji,jj) * e1e2t(ji,jj)                               &
                    &        + ( pv(ji,jj,jk) + pv(ji,jj-1,jk) ) * di_e2u_2(ji,jj)      &
                    &        - ( pu(ji,jj,jk) + pu(ji-1,jj,jk) ) * dj_e1v_2(ji,jj)  )   &
                    &     * e3t(ji,jj,jk,Kmm)
            END_2D
         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','dyn_vor: wrong value for kvor')
         END SELECT
         !
         !                                   !==  compute and add the vorticity term trend  =!
         DO_2D( 0, 0, 0, 0 )
            pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + r1_4 * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)                    &
               &                                * (  zwt(ji+1,jj) * ( pv(ji+1,jj,jk) + pv(ji+1,jj-1,jk) )   &
               &                                   + zwt(ji  ,jj) * ( pv(ji  ,jj,jk) + pv(ji  ,jj-1,jk) )   )
               !
            pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) - r1_4 * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)                    &
               &                                * (  zwt(ji,jj+1) * ( pu(ji,jj+1,jk) + pu(ji-1,jj+1,jk) )   &
               &                                   + zwt(ji,jj  ) * ( pu(ji,jj  ,jk) + pu(ji-1,jj  ,jk) )   )
         END_2D
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      !
#if ! defined key_PSYCLONE_2p5p0
      SELECT CASE( kvor )        ! deallocate zwz if necessary
      CASE ( np_RVO , np_CRV )   ;   DEALLOCATE( zwz )
      END SELECT
#endif
      !
   END SUBROUTINE vor_enT


   SUBROUTINE vor_ene( kt, Kmm, kvor, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE vor_ene  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!       and the Sadourny (1975) flux form formulation : conserves the
      !!       horizontal kinetic energy.
      !!         The general trend of momentum is increased due to the vorticity
      !!       term which is given by:
      !!          voru = 1/e1u  mj-1[ (rvor+f)/e3f  mi(e1v*e3v pvv(:,:,:,Kmm)) ]
      !!          vorv = 1/e2v  mi-1[ (rvor+f)/e3f  mj(e2u*e3u puu(:,:,:,Kmm)) ]
      !!       where rvor is the relative vorticity
      !!
      !! ** Action : - Update (pu_rhs,pv_rhs) with the now vorticity term trend
      !!
      !! References : Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt          ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kmm              ! ocean time level index
      INTEGER                         , INTENT(in   ) ::   kvor        ! total, planetary, relative, or metric
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu, pv    ! now velocities
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs    ! total v-trend
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   zx1, zy1, zx2, zy2, ze3f, zmsk   ! local scalars
      REAL(wp), DIMENSION(T2D(1)) ::   zwx, zwy, zwz   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:vor_ene : vorticity term: energy conserving scheme'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         ENDIF
      ENDIF
      !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         SELECT CASE( kvor )                 !==  vorticity considered  ==!
         CASE ( np_COR )                           !* Coriolis (planetary vorticity)
            DO_2D( 1, 0, 1, 0 )
               zwz(ji,jj) = ff_f(ji,jj)
            END_2D
         CASE ( np_RVO )                           !* relative vorticity
            DO_2D( 1, 0, 1, 0 )
               zwz(ji,jj) = (  ( e2v(ji+1,jj  ) * pv(ji+1,jj  ,jk) - e2v(ji,jj) * pv(ji,jj,jk) )    &                 ! add () for
                  &          - ( e1u(ji  ,jj+1) * pu(ji  ,jj+1,jk) - e1u(ji,jj) * pu(ji,jj,jk) )  ) * r1_e1e2f(ji,jj) ! NP repro
            END_2D
            IF( ln_dynvor_msk ) THEN                     ! mask the relative vorticity
               DO_2D( 1, 0, 1, 0 )
                  zwz(ji,jj) = zwz(ji,jj) * fmask(ji,jj,jk)
               END_2D
            ENDIF
         CASE ( np_MET )                           !* metric term
            DO_2D( 1, 0, 1, 0 )
               zwz(ji,jj) = ( pv(ji+1,jj  ,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &       - ( pu(ji  ,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)
            END_2D
         CASE ( np_CRV )                           !* Coriolis + relative vorticity
            DO_2D( 1, 0, 1, 0 )
               zwz(ji,jj) = ff_f(ji,jj) + (  ( e2v(ji+1,jj) * pv(ji+1,jj,jk) - e2v(ji,jj) * pv(ji,jj,jk) ) & ! add () for NP repro
                  &                        - ( e1u(ji,jj+1) * pu(ji,jj+1,jk) - e1u(ji,jj) * pu(ji,jj,jk) )  ) * r1_e1e2f(ji,jj)
            END_2D
            IF( ln_dynvor_msk ) THEN                     ! mask the relative vorticity (NOT the Coriolis term)
               DO_2D( 1, 0, 1, 0 )
                  zwz(ji,jj) = ( zwz(ji,jj) - ff_f(ji,jj) ) * fmask(ji,jj,jk) + ff_f(ji,jj)
               END_2D
            ENDIF
         CASE ( np_CME )                           !* Coriolis + metric
            DO_2D( 1, 0, 1, 0 )
               zwz(ji,jj) = ff_f(ji,jj) + ( pv(ji+1,jj  ,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &                     - ( pu(ji  ,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)
            END_2D
         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','dyn_vor: wrong value for kvor'  )
         END SELECT
         !
#if defined key_qco   ||   defined key_linssh
         DO_2D( 1, 0, 1, 0 )                 !==  potential vorticity  ==!   (key_qco)
            zwz(ji,jj) = zwz(ji,jj) / e3f_vor(ji,jj,jk)
         END_2D
#else
         SELECT CASE( nn_e3f_typ  )           !==  potential vorticity  ==!
         CASE ( 0 )                                   ! original formulation  (masked averaging of e3t divided by 4)
            DO_2D( 1, 0, 1, 0 )
               ze3f = (  ( e3t(ji  ,jj+1,jk,Kmm)*tmask(ji  ,jj+1,jk)     &   ! need additional () for
                  &    +   e3t(ji+1,jj  ,jk,Kmm)*tmask(ji+1,jj  ,jk) )   &   ! reproducibility around NP
                  &    + ( e3t(ji  ,jj  ,jk,Kmm)*tmask(ji  ,jj  ,jk)     &
                  &    +   e3t(ji+1,jj+1,jk,Kmm)*tmask(ji+1,jj+1,jk) )  )
               IF( ze3f /= 0._wp ) THEN   ;   zwz(ji,jj) = zwz(ji,jj) * 4._wp / ze3f
               ELSE                       ;   zwz(ji,jj) = 0._wp
               ENDIF
            END_2D
         CASE ( 1 )                                   ! new formulation  (masked averaging of e3t divided by the sum of mask)
            DO_2D( 1, 0, 1, 0 )
               ze3f = (   ( e3t(ji  ,jj+1,jk,Kmm)*tmask(ji  ,jj+1,jk)     &   ! need additional () for
                  &     +   e3t(ji+1,jj  ,jk,Kmm)*tmask(ji+1,jj  ,jk) )   &   ! reproducibility around NP
                  &     + ( e3t(ji  ,jj  ,jk,Kmm)*tmask(ji  ,jj  ,jk)     &
                  &     +   e3t(ji+1,jj+1,jk,Kmm)*tmask(ji+1,jj+1,jk) )   )
               zmsk = (   (tmask(ji,jj+1,jk) + tmask(ji+1,jj  ,jk))   &
                  &     + (tmask(ji,jj  ,jk) + tmask(ji+1,jj+1,jk))   )
               IF( ze3f /= 0._wp ) THEN   ;   zwz(ji,jj) = zwz(ji,jj) * zmsk / ze3f
               ELSE                       ;   zwz(ji,jj) = 0._wp
               ENDIF
            END_2D
         END SELECT
#endif
         !                                   !==  horizontal fluxes  ==!
         DO_2D( 1, 1, 1, 1 )
            zwx(ji,jj) = e2u(ji,jj) * e3u(ji,jj,jk,Kmm) * pu(ji,jj,jk)
            zwy(ji,jj) = e1v(ji,jj) * e3v(ji,jj,jk,Kmm) * pv(ji,jj,jk)
         END_2D
         !
         !                                   !==  compute and add the vorticity term trend  =!
         DO_2D( 0, 0, 0, 0 )
            !zy1 = zwy(ji,jj-1) + zwy(ji+1,jj-1)
            !zy2 = zwy(ji,jj  ) + zwy(ji+1,jj  )
            !zx1 = zwx(ji-1,jj) + zwx(ji-1,jj+1)
            !zx2 = zwx(ji  ,jj) + zwx(ji  ,jj+1)
            pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + r1_4 * r1_e1u(ji,jj) * ( (zwz(ji  ,jj-1)*zwy(ji+1,jj-1) + zwz(ji,jj)*zwy(ji  ,jj  ) ) &
                             &                                            +(zwz(ji  ,jj-1)*zwy(ji  ,jj-1) + zwz(ji,jj)*zwy(ji+1,jj  ) ) )
            pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) - r1_4 * r1_e2v(ji,jj) * ( (zwz(ji-1,jj  )*zwx(ji-1,jj  ) + zwz(ji,jj)*zwx(ji  ,jj+1) ) &
                             &                                            +(zwz(ji-1,jj  )*zwx(ji-1,jj+1) + zwz(ji,jj)*zwx(ji  ,jj  ) ) )
         END_2D
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
   END SUBROUTINE vor_ene


   SUBROUTINE vor_ens( kt, Kmm, kvor, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE vor_ens  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      and the Sadourny (1975) flux FORM formulation : conserves the
      !!      potential enstrophy of a horizontally non-divergent flow. the
      !!      trend of the vorticity term is given by:
      !!          voru = 1/e1u  mj-1[ (rvor+f)/e3f ]  mj-1[ mi(e1v*e3v pvv(:,:,:,Kmm)) ]
      !!          vorv = 1/e2v  mi-1[ (rvor+f)/e3f ]  mi-1[ mj(e2u*e3u puu(:,:,:,Kmm)) ]
      !!      Add this trend to the general momentum trend:
      !!          (u(rhs),v(Krhs)) = (u(rhs),v(Krhs)) + ( voru , vorv )
      !!
      !! ** Action : - Update (pu_rhs,pv_rhs)) arrays with the now vorticity term trend
      !!
      !! References : Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt          ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kmm              ! ocean time level index
      INTEGER                         , INTENT(in   ) ::   kvor        ! total, planetary, relative, or metric
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu, pv    ! now velocities
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs    ! total v-trend
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zuav, zvau, ze3f, zmsk   ! local scalars
      REAL(wp), DIMENSION(T2D(1)) ::   zwx, zwy, zwz   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:vor_ens : vorticity term: enstrophy conserving scheme'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         ENDIF
      ENDIF
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         SELECT CASE( kvor )                 !==  vorticity considered  ==!
         CASE ( np_COR )                           !* Coriolis (planetary vorticity)
            DO_2D( 1, 0, 1, 0 )
               zwz(ji,jj) = ff_f(ji,jj)
            END_2D
         CASE ( np_RVO )                           !* relative vorticity
            DO_2D( 1, 0, 1, 0 )
               zwz(ji,jj) = (  ( e2v(ji+1,jj  ) * pv(ji+1,jj  ,jk) - e2v(ji,jj) * pv(ji,jj,jk) )    &   ! add () for NP repro
                  &          - ( e1u(ji  ,jj+1) * pu(ji  ,jj+1,jk) - e1u(ji,jj) * pu(ji,jj,jk) )  ) * r1_e1e2f(ji,jj)
            END_2D
            IF( ln_dynvor_msk ) THEN                     ! mask the relative vorticity
               DO_2D( 1, 0, 1, 0 )
                  zwz(ji,jj) = zwz(ji,jj) * fmask(ji,jj,jk)
               END_2D
            ENDIF
         CASE ( np_MET )                           !* metric term
            DO_2D( 1, 0, 1, 0 )
               zwz(ji,jj) = ( pv(ji+1,jj  ,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &       - ( pu(ji  ,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)
            END_2D
         CASE ( np_CRV )                           !* Coriolis + relative vorticity
            DO_2D( 1, 0, 1, 0 )
               zwz(ji,jj) = ff_f(ji,jj) + ( ( e2v(ji+1,jj  ) * pv(ji+1,jj  ,jk) - e2v(ji,jj) * pv(ji,jj,jk) ) &! add () for NP repro
                  &                        -( e1u(ji  ,jj+1) * pu(ji  ,jj+1,jk) - e1u(ji,jj) * pu(ji,jj,jk) )  ) * r1_e1e2f(ji,jj)
            END_2D
            IF( ln_dynvor_msk ) THEN                     ! mask the relative vorticity (NOT the Coriolis term)
               DO_2D( 1, 0, 1, 0 )
                  zwz(ji,jj) = ( zwz(ji,jj) - ff_f(ji,jj) ) * fmask(ji,jj,jk) + ff_f(ji,jj)
               END_2D
            ENDIF
         CASE ( np_CME )                           !* Coriolis + metric
            DO_2D( 1, 0, 1, 0 )
               zwz(ji,jj) = ff_f(ji,jj) + ( pv(ji+1,jj  ,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &                     - ( pu(ji  ,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)
            END_2D
         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','dyn_vor: wrong value for kvor'  )
         END SELECT
         !
         !
#if defined key_qco   ||   defined key_linssh
         DO_2D( 1, 0, 1, 0 )                 !==  potential vorticity  ==!   (key_qco)
            zwz(ji,jj) = zwz(ji,jj) / e3f_vor(ji,jj,jk)
         END_2D
#else
         SELECT CASE( nn_e3f_typ )           !==  potential vorticity  ==!
         CASE ( 0 )                                   ! original formulation  (masked averaging of e3t divided by 4)
            DO_2D( 1, 0, 1, 0 )
               ze3f = (  ( e3t(ji  ,jj+1,jk,Kmm)*tmask(ji  ,jj+1,jk)     &   ! need additional () for
                  &    +   e3t(ji+1,jj  ,jk,Kmm)*tmask(ji+1,jj  ,jk) )   &   ! reproducibility around NP
                  &    + ( e3t(ji  ,jj  ,jk,Kmm)*tmask(ji  ,jj  ,jk)     &
                  &    +   e3t(ji+1,jj+1,jk,Kmm)*tmask(ji+1,jj+1,jk) )  )
               IF( ze3f /= 0._wp ) THEN   ;   zwz(ji,jj) = zwz(ji,jj) * 4._wp / ze3f
               ELSE                       ;   zwz(ji,jj) = 0._wp
               ENDIF
            END_2D
         CASE ( 1 )                                   ! new formulation  (masked averaging of e3t divided by the sum of mask)
            DO_2D( 1, 0, 1, 0 )
               ze3f = (   ( e3t(ji  ,jj+1,jk,Kmm)*tmask(ji  ,jj+1,jk)     &   ! need additional () for
                  &     +   e3t(ji+1,jj  ,jk,Kmm)*tmask(ji+1,jj  ,jk) )   &   ! reproducibility around NP
                  &     + ( e3t(ji  ,jj  ,jk,Kmm)*tmask(ji  ,jj  ,jk)     &
                  &     +   e3t(ji+1,jj+1,jk,Kmm)*tmask(ji+1,jj+1,jk) )   )
               zmsk = (   (tmask(ji,jj+1,jk) + tmask(ji+1,jj  ,jk))   &
                  &     + (tmask(ji,jj  ,jk) + tmask(ji+1,jj+1,jk))   )
               IF( ze3f /= 0._wp ) THEN   ;   zwz(ji,jj) = zwz(ji,jj) * zmsk / ze3f
               ELSE                       ;   zwz(ji,jj) = 0._wp
               ENDIF
            END_2D
         END SELECT
#endif
         !                                   !==  horizontal fluxes  ==!
         DO_2D( 1, 1, 1, 1 )
            zwx(ji,jj) = e2u(ji,jj) * e3u(ji,jj,jk,Kmm) * pu(ji,jj,jk)
            zwy(ji,jj) = e1v(ji,jj) * e3v(ji,jj,jk,Kmm) * pv(ji,jj,jk)
         END_2D
         !
         !                                   !==  compute and add the vorticity term trend  =!
         DO_2D( 0, 0, 0, 0 )
            zuav = r1_8 * r1_e1u(ji,jj) * (  ( zwy(ji  ,jj-1) + zwy(ji+1,jj  ) )  &   ! need additional () for
               &                           + ( zwy(ji  ,jj  ) + zwy(ji+1,jj-1) )  )   ! reproducibility around NP
            zvau =-r1_8 * r1_e2v(ji,jj) * (  ( zwx(ji-1,jj  ) + zwx(ji  ,jj+1) )  &
               &                           + ( zwx(ji  ,jj  ) + zwx(ji-1,jj+1) )  )
            pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zuav * ( zwz(ji  ,jj-1) + zwz(ji,jj) )
            pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zvau * ( zwz(ji-1,jj  ) + zwz(ji,jj) )
         END_2D
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
   END SUBROUTINE vor_ens


   SUBROUTINE vor_een( kt, Kmm, kvor, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE vor_een  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      and the Arakawa and Lamb (1980) flux form formulation : conserves
      !!      both the horizontal kinetic energy and the potential enstrophy
      !!      when horizontal divergence is zero (see the NEMO documentation)
      !!      Add this trend to the general momentum trend (pu_rhs,pv_rhs).
      !!
      !! ** Action : - Update (pu_rhs,pv_rhs) with the now vorticity term trend
      !!
      !! References : Arakawa and Lamb 1980, Mon. Wea. Rev., 109, 18-36
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt          ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kmm              ! ocean time level index
      INTEGER                         , INTENT(in   ) ::   kvor        ! total, planetary, relative, or metric
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu, pv    ! now velocities
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs    ! total v-trend
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   ierr         ! local integer
      REAL(wp) ::   zua, zva     ! local scalars
      REAL(wp) ::   zmsk, ze3f   ! local scalars
      REAL(wp), DIMENSION(T2D(1)) ::   z1_e3f
      REAL(wp), DIMENSION(T2D(1)) ::   zwx , zwy
      REAL(wp), DIMENSION(T2D(1)) ::   ztnw, ztne, ztsw, ztse
      REAL(wp), DIMENSION(T2D(1)) ::   zwz   ! 3D workspace, jpkm1 -> jpkm1 -> avoid lbc_lnk on jpk that is not defined
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:vor_een : vorticity term: energy and enstrophy conserving scheme'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         ENDIF
      ENDIF
      !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
#if defined key_qco   ||   defined key_linssh
         DO_2D( 1, 1, 1, 1 )                 ! == reciprocal of e3 at F-point (key_qco)
            z1_e3f(ji,jj) = 1._wp / e3f_vor(ji,jj,jk)
         END_2D
#else
         SELECT CASE( nn_e3f_typ )           ! == reciprocal of e3 at F-point
         CASE ( 0 )                                   ! original formulation  (masked averaging of e3t divided by 4)
            DO_2D( 1, 1, 1, 1 )
               ze3f = (  ( e3t(ji  ,jj+1,jk,Kmm)*tmask(ji  ,jj+1,jk)     &   ! need additional () for
                  &    +   e3t(ji+1,jj  ,jk,Kmm)*tmask(ji+1,jj  ,jk) )   &   ! reproducibility around NP
                  &    + ( e3t(ji  ,jj  ,jk,Kmm)*tmask(ji  ,jj  ,jk)     &
                  &    +   e3t(ji+1,jj+1,jk,Kmm)*tmask(ji+1,jj+1,jk) )  )
               IF( ze3f /= 0._wp ) THEN   ;   z1_e3f(ji,jj) = 4._wp / ze3f
               ELSE                       ;   z1_e3f(ji,jj) = 0._wp
               ENDIF
            END_2D
         CASE ( 1 )                                   ! new formulation  (masked averaging of e3t divided by the sum of mask)
            DO_2D( 1, 1, 1, 1 )
               ze3f = (  ( e3t(ji  ,jj+1,jk,Kmm)*tmask(ji  ,jj+1,jk)     &   ! need additional () for
                  &    +   e3t(ji+1,jj  ,jk,Kmm)*tmask(ji+1,jj  ,jk) )   &   ! reproducibility around NP
                  &    + ( e3t(ji  ,jj  ,jk,Kmm)*tmask(ji  ,jj  ,jk)     &
                  &    +   e3t(ji+1,jj+1,jk,Kmm)*tmask(ji+1,jj+1,jk) )  )
               zmsk = ((tmask(ji,jj+1,jk) + tmask(ji+1,jj  ,jk))   &
                  &  + (tmask(ji,jj  ,jk) + tmask(ji+1,jj+1,jk))  )
               IF( ze3f /= 0._wp ) THEN   ;   z1_e3f(ji,jj) = zmsk / ze3f
               ELSE                       ;   z1_e3f(ji,jj) = 0._wp
               ENDIF
            END_2D
         END SELECT
#endif
         !
         SELECT CASE( kvor )                 !==  vorticity considered  ==!
         !
         CASE ( np_COR )                           !* Coriolis (planetary vorticity)
            DO_2D( 1, 1, 1, 1 )
               zwz(ji,jj) = ff_f(ji,jj) * z1_e3f(ji,jj)
            END_2D
         CASE ( np_RVO )                           !* relative vorticity
            DO_2D( 1, 1, 1, 1 )
               zwz(ji,jj) = ( ( e2v(ji+1,jj  ) * pv(ji+1,jj,jk) - e2v(ji,jj) * pv(ji,jj,jk) )  &   ! add () for NP repro
                  &         - ( e1u(ji  ,jj+1) * pu(ji,jj+1,jk) - e1u(ji,jj) * pu(ji,jj,jk) )  ) * r1_e1e2f(ji,jj)*z1_e3f(ji,jj)
            END_2D
            IF( ln_dynvor_msk ) THEN                     ! mask the relative vorticity
               DO_2D( 1, 1, 1, 1 )
                  zwz(ji,jj) = zwz(ji,jj) * fmask(ji,jj,jk)
               END_2D
            ENDIF
         CASE ( np_MET )                           !* metric term
            DO_2D( 1, 1, 1, 1 )
               zwz(ji,jj) = (   ( pv(ji+1,jj,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &           - ( pu(ji,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)   ) * z1_e3f(ji,jj)
            END_2D
         CASE ( np_CRV )                           !* Coriolis + relative vorticity
            IF( ln_dynvor_msk ) THEN 
               DO_2D( 1, 1, 1, 1 )
                  zwz(ji,jj) = (  ff_f(ji,jj) + ( ( e2v(ji+1,jj  ) * pv(ji+1,jj,jk) - e2v(ji,jj) * pv(ji,jj,jk) )   & ! add () for
                     &                          - ( e1u(ji  ,jj+1) * pu(ji,jj+1,jk) - e1u(ji,jj) * pu(ji,jj,jk) )   & ! NP repro
                     &                          ) * r1_e1e2f(ji,jj) * fmask(ji,jj,jk)   ) * z1_e3f(ji,jj)
               END_2D
            ELSE
               DO_2D( 1, 1, 1, 1 )
                  zwz(ji,jj) = (  ff_f(ji,jj) + ( ( e2v(ji+1,jj  ) * pv(ji+1,jj,jk) - e2v(ji,jj) * pv(ji,jj,jk) )   & ! add () for
                     &                          - ( e1u(ji  ,jj+1) * pu(ji,jj+1,jk) - e1u(ji,jj) * pu(ji,jj,jk) )   & ! NP repro
                     &                          ) * r1_e1e2f(ji,jj)   ) * z1_e3f(ji,jj)
               END_2D
            ENDIF
         CASE ( np_CME )                           !* Coriolis + metric
            DO_2D( 1, 1, 1, 1 )
               zwz(ji,jj) = (   ff_f(ji,jj) + ( pv(ji+1,jj  ,jk) + pv(ji,jj,jk) ) * di_e2v_2e1e2f(ji,jj)   &
                  &                         - ( pu(ji  ,jj+1,jk) + pu(ji,jj,jk) ) * dj_e1u_2e1e2f(ji,jj)   ) * z1_e3f(ji,jj)
            END_2D
         CASE DEFAULT                                             ! error
            CALL ctl_stop('STOP','dyn_vor: wrong value for kvor'  )
         END SELECT
         !
         !                                   !==  horizontal fluxes  ==!
         DO_2D( 1, 1, 1, 1 )
            zwx(ji,jj) = e2u(ji,jj) * e3u(ji,jj,jk,Kmm) * pu(ji,jj,jk)
            zwy(ji,jj) = e1v(ji,jj) * e3v(ji,jj,jk,Kmm) * pv(ji,jj,jk)
         END_2D
         !
         !                                   !==  compute and add the vorticity term trend  =!
         DO_2D( 0, 1, 0, 1 )
            ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
            ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
            ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
            ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
         END_2D
         !
         DO_2D( 0, 0, 0, 0 )   ! add () for NP repro
            zua = + r1_12 * r1_e1u(ji,jj) * (  ( ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )   & ! add () for
               &                             + ( ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztnw(ji+1,jj) * zwy(ji+1,jj  ) ) )   ! NP repro
            zva = - r1_12 * r1_e2v(ji,jj) * (  ( ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )   &
               &                             + ( ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztse(ji,jj+1) * zwx(ji  ,jj+1) ) )
            pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zua
            pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zva
         END_2D
      END DO
      !                                                ! ===============
      !                                                !   End of slab
      !                                                ! ===============
   END SUBROUTINE vor_een


   SUBROUTINE vor_weno( kt, Kmm, pu, pv, pu_rrd, pv_rrd )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE vor_weno  ***
      !!
      !! ** Purpose :   Compute the relative vorticity term of the vector invariant
      !!      momentum equations with an upwind (WENO) reconstruction of the relative
      !!      vorticity, and add it to the general trend:
      !!                    du/dt = + vh * zeta_u ,   dv/dt = - uh * zeta_v
      !!      The Coriolis term is computed separately by the chosen vorticity scheme.
      !!
      !! ** Method  :   This is the vorticity term of Oceananigans WENOVectorInvariant
      !!      (Silvestri et al. 2024). The relative vorticity zeta is computed at
      !!      f-points (masked by fmask if ln_dynvor_msk=T). At a u-point the transverse
      !!      velocity vh is the 4-point average of e1v*v divided by e1u, and zeta_u is
      !!      reconstructed along j from the f-points upstream of vh (WENO of order
      !!      nn_dynvor_weno); at a v-point uh is the average of e2u*u divided by e2v,
      !!      and zeta_v is reconstructed along i from the f-points upstream of uh.
      !!      The WENO weights are computed from the smoothness of the velocities
      !!      interpolated at the f-points of the stencil (u averaged along j, v along i).
      !!      Near land the order is reduced (5 -> 3 -> 1) so that the stencil only
      !!      uses active f-points (f-points with at least one ocean t-point around).
      !!        The reconstruction needs nn_dynvor_weno/2 + 1 points of halo.
      !!
      !! ** Action : - Update (pu_rrd,pv_rrd) with the now relative vorticity term trend
      !!
      !! References : Silvestri, S., G. L. Wagner, et al., 2024, J. Adv. Model. Earth Syst., 16, e2023MS004130.
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt               ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kmm              ! ocean time level index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pu, pv           ! now velocities
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rrd, pv_rrd   ! total v-trend
      !
      INTEGER  ::   ji, jj, jk, ib    ! dummy loop indices, half width of the stencil
      REAL(wp) ::   zvel              ! local scalar
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zwz        ! relative vorticity at f-points
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zuf, zvf   ! u (averaged along j) and v (averaged along i) at f-points
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zmf        ! =1 at active f-points
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:vor_weno : vorticity term: WENO upwind relative vorticity, order ', nn_dynvor_weno
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
         ENDIF
      ENDIF
      !
      ib = ( nn_dynvor_weno + 1 ) / 2               ! f-points needed on each side of a velocity point
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO_2D( ib, ib-1, ib, ib-1 )                   !==  f-point quantities  ==!
            zwz(ji,jj) = (  ( e2v(ji+1,jj  ) * pv(ji+1,jj,jk) - e2v(ji,jj) * pv(ji,jj,jk) )  &   ! add () for NP repro
               &          - ( e1u(ji  ,jj+1) * pu(ji,jj+1,jk) - e1u(ji,jj) * pu(ji,jj,jk) )  ) * r1_e1e2f(ji,jj)
            zuf(ji,jj) = 0.5_wp * ( pu(ji,jj,jk) + pu(ji  ,jj+1,jk) )
            zvf(ji,jj) = 0.5_wp * ( pv(ji,jj,jk) + pv(ji+1,jj  ,jk) )
            zmf(ji,jj) = MAX( tmask(ji,jj,jk), tmask(ji+1,jj,jk), tmask(ji,jj+1,jk), tmask(ji+1,jj+1,jk) )
         END_2D
         IF( ln_dynvor_msk ) THEN                      ! mask the relative vorticity
            DO_2D( ib, ib-1, ib, ib-1 )
               zwz(ji,jj) = zwz(ji,jj) * fmask(ji,jj,jk)
            END_2D
         ENDIF
         !
         IF( nn_dynvor_weno == 5 ) THEN                !==  vorticity trend  ==!
            DO_2D( 0, 0, 0, 0 )
               zvel = r1_4 * (  ( e1v(ji,jj  ) * pv(ji,jj  ,jk) + e1v(ji+1,jj  ) * pv(ji+1,jj  ,jk) )   &   ! vh at u-point
                  &           + ( e1v(ji,jj-1) * pv(ji,jj-1,jk) + e1v(ji+1,jj-1) * pv(ji+1,jj-1,jk) )  ) * r1_e1u(ji,jj)
               pu_rrd(ji,jj,jk) = pu_rrd(ji,jj,jk) + zvel * weno_vor5( zwz(ji,jj-3:jj+2), zuf(ji,jj-3:jj+2),   &
                  &                                                    zvf(ji,jj-3:jj+2), zmf(ji,jj-3:jj+2), zvel )
               zvel = r1_4 * (  ( e2u(ji-1,jj  ) * pu(ji-1,jj  ,jk) + e2u(ji,jj  ) * pu(ji,jj  ,jk) )   &   ! uh at v-point
                  &           + ( e2u(ji-1,jj+1) * pu(ji-1,jj+1,jk) + e2u(ji,jj+1) * pu(ji,jj+1,jk) )  ) * r1_e2v(ji,jj)
               pv_rrd(ji,jj,jk) = pv_rrd(ji,jj,jk) - zvel * weno_vor5( zwz(ji-3:ji+2,jj), zuf(ji-3:ji+2,jj),   &
                  &                                                    zvf(ji-3:ji+2,jj), zmf(ji-3:ji+2,jj), zvel )
            END_2D
         ELSE
            DO_2D( 0, 0, 0, 0 )
               zvel = r1_4 * (  ( e1v(ji,jj  ) * pv(ji,jj  ,jk) + e1v(ji+1,jj  ) * pv(ji+1,jj  ,jk) )   &   ! vh at u-point
                  &           + ( e1v(ji,jj-1) * pv(ji,jj-1,jk) + e1v(ji+1,jj-1) * pv(ji+1,jj-1,jk) )  ) * r1_e1u(ji,jj)
               pu_rrd(ji,jj,jk) = pu_rrd(ji,jj,jk) + zvel * weno_vor3( zwz(ji,jj-2:jj+1), zuf(ji,jj-2:jj+1),   &
                  &                                                    zvf(ji,jj-2:jj+1), zmf(ji,jj-2:jj+1), zvel )
               zvel = r1_4 * (  ( e2u(ji-1,jj  ) * pu(ji-1,jj  ,jk) + e2u(ji,jj  ) * pu(ji,jj  ,jk) )   &   ! uh at v-point
                  &           + ( e2u(ji-1,jj+1) * pu(ji-1,jj+1,jk) + e2u(ji,jj+1) * pu(ji,jj+1,jk) )  ) * r1_e2v(ji,jj)
               pv_rrd(ji,jj,jk) = pv_rrd(ji,jj,jk) - zvel * weno_vor3( zwz(ji-2:ji+1,jj), zuf(ji-2:ji+1,jj),   &
                  &                                                    zvf(ji-2:ji+1,jj), zmf(ji-2:ji+1,jj), zvel )
            END_2D
         ENDIF
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
   END SUBROUTINE vor_weno


   SUBROUTINE dyn_vor_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_vor_init  ***
      !!
      !! ** Purpose :   Control the consistency between cpp options for
      !!              tracer advection schemes
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk    ! dummy loop indices
      INTEGER ::   ioptio, ios   ! local integer
      REAL(wp) ::   zmsk    ! local scalars
      !!
      NAMELIST/namdyn_vor/ ln_dynvor_ens, ln_dynvor_ene, ln_dynvor_enT,   &
         &                 ln_dynvor_een, nn_e3f_typ   , ln_dynvor_msk,   &
         &                 ln_dynvor_weno, nn_dynvor_weno
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_vor_init : vorticity term : read namelist and control the consistency'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      READ_NML_REF(numnam,namdyn_vor)
      READ_NML_CFG(numnam,namdyn_vor)
      IF(lwm) WRITE ( numond, namdyn_vor )
      !
      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*) '   Namelist namdyn_vor : choice of the vorticity term scheme'
         WRITE(numout,*) '      enstrophy conserving scheme                    ln_dynvor_ens = ', ln_dynvor_ens
         WRITE(numout,*) '      f-point energy conserving scheme               ln_dynvor_ene = ', ln_dynvor_ene
         WRITE(numout,*) '      t-point energy conserving scheme               ln_dynvor_enT = ', ln_dynvor_enT
         WRITE(numout,*) '      enstrophy and energy conserving scheme         ln_dynvor_een = ', ln_dynvor_een
         WRITE(numout,*) '         e3f = averaging /4 (=0) or /sum(tmask) (=1)    nn_e3f_typ = ', nn_e3f_typ
         WRITE(numout,*) '      masked (=T) or unmasked(=F) vorticity          ln_dynvor_msk = ', ln_dynvor_msk
         WRITE(numout,*) '      WENO upwind relative vorticity               ln_dynvor_weno = ', ln_dynvor_weno
         WRITE(numout,*) '         order of the WENO reconstruction          nn_dynvor_weno = ', nn_dynvor_weno
      ENDIF

!!gm  this should be removed when choosing a unique strategy for fmask at the coast
      ! If energy, enstrophy or mixed advection of momentum in vector form change the value for masks
      ! at angles with three ocean points and one land point
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '      change fmask value in the angles (T)           ln_vorlat = ', ln_vorlat
      IF( ln_vorlat .AND. ( ln_dynvor_ene .OR. ln_dynvor_ens ) ) THEN
         DO_3D( 1, 0, 1, 0, 1, jpk )
            IF(    tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk)              &
               & + tmask(ji,jj  ,jk) + tmask(ji+1,jj  ,jk) == 3._wp )   fmask(ji,jj,jk) = 1._wp
         END_3D
         !
         CALL lbc_lnk( 'dynvor', fmask, 'F', 1._wp )      ! Lateral boundary conditions on fmask
         !
      ENDIF
!!gm end

      ioptio = 0                     ! type of scheme for vorticity (set nvor_scheme)
      IF( ln_dynvor_ens ) THEN   ;   ioptio = ioptio + 1   ;   nvor_scheme = np_ENS   ;   ENDIF
      IF( ln_dynvor_ene ) THEN   ;   ioptio = ioptio + 1   ;   nvor_scheme = np_ENE   ;   ENDIF
      IF( ln_dynvor_enT ) THEN   ;   ioptio = ioptio + 1   ;   nvor_scheme = np_ENT   ;   ENDIF
      IF( ln_dynvor_een ) THEN   ;   ioptio = ioptio + 1   ;   nvor_scheme = np_EEN   ;   ENDIF
      !
      IF( ioptio /= 1 ) CALL ctl_stop( ' use ONE and ONLY one vorticity scheme' )
      !
      IF( ln_dynvor_weno ) THEN      ! WENO relative vorticity: vector form only, the chosen scheme computes the Coriolis term
         IF( n_dynadv /= np_VEC_c2 )   CALL ctl_stop( 'dyn_vor_init: ln_dynvor_weno requires the vector form (ln_dynadv_vec = T)' )
         IF( nn_dynvor_weno /= 3 .AND. nn_dynvor_weno /= 5 )   CALL ctl_stop( 'dyn_vor_init: nn_dynvor_weno must be 3 or 5' )
         IF( nn_hls < ( nn_dynvor_weno + 1 ) / 2 )   &
            &   CALL ctl_stop( 'dyn_vor_init: 5th order WENO relative vorticity requires nn_hls = 3 (namelist nammpp)' )
         IF( ln_sdw .AND. ln_vortex_force )   CALL ctl_warn( 'dyn_vor_init: with ln_dynvor_weno, the vortex force uses the chosen scheme' )
      ENDIF
      !
      IF(lwp) WRITE(numout,*)        ! type of calculated vorticity (set ncor, nrvm, ntot)
      ncor = np_COR                       ! planetary vorticity
      SELECT CASE( n_dynadv )
      CASE( np_LIN_dyn )
         IF(lwp) WRITE(numout,*) '   ==>>>   linear dynamics : total vorticity = Coriolis'
         nrvm = np_COR        ! planetary vorticity
         ntot = np_COR        !     -         -
      CASE( np_VEC_c2  )
         IF(lwp) WRITE(numout,*) '   ==>>>   vector form dynamics : total vorticity = Coriolis + relative vorticity'
         nrvm = np_RVO        ! relative vorticity
         ntot = np_CRV        ! relative + planetary vorticity
      CASE( np_FLX_c2 , np_FLX_up3 )
         IF(lwp) WRITE(numout,*) '   ==>>>   flux form dynamics : total vorticity = Coriolis + metric term'
         nrvm = np_MET        ! metric term
         ntot = np_CME        ! Coriolis + metric term
         !
         SELECT CASE( nvor_scheme )    ! pre-computed gradients for the metric term:
         CASE( np_ENT )                      !* T-point metric term :   pre-compute di(e2u)/2 and dj(e1v)/2
            ALLOCATE( di_e2u_2(jpi,jpj), dj_e1v_2(jpi,jpj) )
            DO_2D( 0, 0, 0, 0 )
               di_e2u_2(ji,jj) = ( e2u(ji,jj) - e2u(ji-1,jj  ) ) * 0.5_wp
               dj_e1v_2(ji,jj) = ( e1v(ji,jj) - e1v(ji  ,jj-1) ) * 0.5_wp
            END_2D
            CALL lbc_lnk( 'dynvor', di_e2u_2, 'T', -1.0_wp , dj_e1v_2, 'T', -1.0_wp )   ! Lateral boundary conditions
            !
         CASE DEFAULT                        !* F-point metric term :   pre-compute di(e2u)/(2*e1e2f) and dj(e1v)/(2*e1e2f)
            ALLOCATE( di_e2v_2e1e2f(jpi,jpj), dj_e1u_2e1e2f(jpi,jpj) )
            DO_2D( 0, 0, 0, 0 )
               di_e2v_2e1e2f(ji,jj) = ( e2v(ji+1,jj  ) - e2v(ji,jj) )  * 0.5 * r1_e1e2f(ji,jj)
               dj_e1u_2e1e2f(ji,jj) = ( e1u(ji  ,jj+1) - e1u(ji,jj) )  * 0.5 * r1_e1e2f(ji,jj)
            END_2D
            CALL lbc_lnk( 'dynvor', di_e2v_2e1e2f, 'F', -1.0_wp , dj_e1u_2e1e2f, 'F', -1.0_wp )   ! Lateral boundary conditions
         END SELECT
         !
      END SELECT
#if defined key_qco   ||   defined key_linssh
      SELECT CASE( nvor_scheme )    ! qco or linssh cases : pre-computed a specific e3f_0 for some vorticity schemes
      CASE( np_ENS , np_ENE , np_EEN )
         !
         ALLOCATE( e3f_0vor(jpi,jpj,jpk) )
         !
         SELECT CASE( nn_e3f_typ )
         CASE ( 0 )                        ! original formulation  (masked averaging of e3t divided by 4)
            DO_3D( 0, 0, 0, 0, 1, jpk )
               e3f_0vor(ji,jj,jk) = (   ( e3t_0(ji  ,jj+1,jk)*tmask(ji  ,jj+1,jk)     &   ! need additional () for
                  &                   +   e3t_0(ji+1,jj  ,jk)*tmask(ji+1,jj  ,jk) )   &   ! reproducibility around NP
                  &                   + ( e3t_0(ji  ,jj  ,jk)*tmask(ji  ,jj  ,jk)     &
                  &                   +   e3t_0(ji+1,jj+1,jk)*tmask(ji+1,jj+1,jk) )   ) * 0.25_wp
            END_3D
         CASE ( 1 )                        ! new formulation  (masked averaging of e3t divided by the sum of mask)
            DO_3D( 0, 0, 0, 0, 1, jpk )
               zmsk = (tmask(ji,jj+1,jk) +tmask(ji+1,jj+1,jk)   &
                  &  + tmask(ji,jj  ,jk) +tmask(ji+1,jj  ,jk)  )
               !
               IF( zmsk /= 0._wp ) THEN
                  e3f_0vor(ji,jj,jk) = (   ( e3t_0(ji  ,jj+1,jk)*tmask(ji  ,jj+1,jk)     &   ! need additional () for
                     &                   +   e3t_0(ji+1,jj  ,jk)*tmask(ji+1,jj  ,jk) )   &   ! reproducibility around NP
                     &                   + ( e3t_0(ji  ,jj  ,jk)*tmask(ji  ,jj  ,jk)     &
                     &                   +   e3t_0(ji+1,jj+1,jk)*tmask(ji+1,jj+1,jk) )   ) / zmsk
               ELSE ; e3f_0vor(ji,jj,jk) = 0._wp
               ENDIF
            END_3D
         END SELECT
         !
         CALL lbc_lnk( 'dynvor', e3f_0vor, 'F', 1._wp )
         !                                 ! insure e3f_0vor /= 0
# if defined key_vco_1d
         DO jk = 1, jpk
            WHERE( e3f_0vor(:,:,jk) == 0._wp )   e3f_0vor(:,:,jk) = e3f_0(:,:,jk)
         END DO
# else
         WHERE( e3f_0vor(:,:,:) == 0._wp )   e3f_0vor(:,:,:) = e3f_0(:,:,:)
# endif
         !
      END SELECT
      !
#endif
      IF(lwp) THEN                   ! Print the choice
         WRITE(numout,*)
         SELECT CASE( nvor_scheme )
         CASE( np_ENS )   ;   WRITE(numout,*) '   ==>>>   enstrophy conserving scheme (ENS)'
         CASE( np_ENE )   ;   WRITE(numout,*) '   ==>>>   energy conserving scheme (Coriolis at F-points) (ENE)'
         CASE( np_ENT )   ;   WRITE(numout,*) '   ==>>>   energy conserving scheme (Coriolis at T-points) (ENT)'
                              IF( ln_dynadv_vec )   CALL ctl_warn('dyn_vor_init: ENT scheme may not work in vector form')
         CASE( np_EEN )   ;   WRITE(numout,*) '   ==>>>   energy and enstrophy conserving scheme (EEN)'
         END SELECT
         IF( ln_dynvor_weno )   WRITE(numout,*) '   ==>>>   used for the Coriolis term only, relative vorticity by WENO of order ', &
            &                                   nn_dynvor_weno
      ENDIF
      !
   END SUBROUTINE dyn_vor_init

   !! * WENO reconstructions (weno_vor5, weno_vor3, ...)
#  include "weno_kernels.h90"

   !!==============================================================================
END MODULE dynvor
