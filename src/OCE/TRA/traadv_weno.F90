MODULE traadv_weno
   !!======================================================================
   !!                       ***  MODULE  traadv_weno  ***
   !! Ocean  tracers:  horizontal & vertical advective trend
   !!======================================================================
   !! History :  5.x  !  2026-10  (G. Wagner)  3rd and 5th order WENO-Z scheme, following Oceananigans.jl
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv_weno  : update the tracer trend with the horizontal and
   !!                   vertical advection trends using a WENO scheme
   !!   weno5_face, weno5, weno3 : WENO-Z reconstructions (weno_kernels.h90)
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE trc_oce        ! share passive tracers/Ocean variables
   USE dom_oce        ! ocean space and time domain
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! tracers trends manager
   USE diaptr         ! poleward transport diagnostics
   USE diaar5         ! AR5 diagnostics
   USE domutl, ONLY : lbnd_ij
   !
   USE iom            ! XIOS library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_weno       ! routine called by traadv.F90

   LOGICAL  ::   l_trd   ! flag to compute trends
   LOGICAL  ::   l_ptr   ! flag to compute poleward transport
   LOGICAL  ::   l_hst   ! flag to compute heat/salt transport

   REAL(wp), PUBLIC, PARAMETER :: pp_stb_thres_tra_weno = 0.85_wp  ! starting Courant number threshold for adaptive implicit vertical advection
   REAL(wp), PUBLIC, PARAMETER :: pp_stb_cstra_tra_weno = 1.00_wp  ! stability constraint

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_weno( kt, kit000, cdtype, pU, pV, pW,             &
      &                     Kmm, pt, kjpt, Krhs, kn_weno_h, kn_weno_v   )
      !!
      INTEGER                                  , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                                  , INTENT(in   ) ::   Kmm, Krhs       ! ocean time level indices
      INTEGER                                  , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                         , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                                  , INTENT(in   ) ::   kjpt            ! number of tracers
      INTEGER                                  , INTENT(in   ) ::   kn_weno_h       ! order of the horizontal scheme (3 or 5)
      INTEGER                                  , INTENT(in   ) ::   kn_weno_v       ! order of the vertical   scheme (3 or 5)
      REAL(wp), DIMENSION(:,:,:               ), INTENT(in   ) ::   pU, pV, pW      ! 3 ocean volume flux components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt,jpt), INTENT(inout) ::   pt              ! tracers and RHS of tracer equation
      !!
      CALL tra_adv_weno_t( kt, kit000, cdtype, pU, pV, pW, lbnd_ij(pU),    &
         &                 Kmm, pt, kjpt, Krhs, kn_weno_h, kn_weno_v      )
   END SUBROUTINE tra_adv_weno


   SUBROUTINE tra_adv_weno_t( kt, kit000, cdtype, pU, pV, pW, ktpuvw,   &
      &                       Kmm, pt, kjpt, Krhs, kn_weno_h, kn_weno_v )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_weno  ***
      !!
      !! ** Purpose :   Compute the now trend due to the advection of tracers
      !!              and add it to the general trend of tracer equations.
      !!
      !! ** Method  :   Weighted Essentially Non-Oscillatory (WENO) scheme.
      !!      The tracer value at a cell face is reconstructed from the cells
      !!      upstream of the face as a convex combination of the reconstructions
      !!      on each candidate stencil (2 stencils of 2 cells at 3rd order,
      !!      3 stencils of 3 cells at 5th order). The weights are the WENO-Z
      !!      weights of Borges et al. (2008), which reduce to the optimal
      !!      (3rd or 5th order upstream) weights where the tracer is smooth
      !!      and favour the smoothest stencils near fronts. The face values
      !!      are multiplied by the volume fluxes, and the advective trend is
      !!      the divergence of these fluxes (flux form, conservative).
      !!        The order is chosen separately in the horizontal (kn_weno_h)
      !!      and in the vertical (kn_weno_v). Where the upstream stencil reaches
      !!      a land point (or the surface/bottom) the reconstruction is reduced
      !!      to the highest order whose stencil is in the ocean: 5th, 3rd, then
      !!      the 1st order upstream value, as Oceananigans does near immersed
      !!      boundaries with WENO(order=5, minimum_buffer_upwind_order=1).
      !!        The whole scheme is evaluated at the RK3 stage time level Kmm.
      !!        The horizontal stencil spans kn_weno_h/2 + 1 points on each side
      !!      of a face: 3rd order needs nn_hls >= 2, 5th order nn_hls >= 3.
      !!
      !! ** Action : - update pt(:,:,:,:,Krhs)  with the now advective tracer trends
      !!             - send trends to trdtra module for further diagnostcs (l_trdtra=T)
      !!             - poleward advective heat and salt transport (ln_diaptr=T)
      !!
      !! References : Jiang, G.-S., C.-W. Shu, 1996, J. Comput. Phys., 126, 202-228.
      !!              Borges, R., M. Carmona, B. Costa, W. S. Don, 2008, J. Comput. Phys., 227, 3191-3211.
      !!              Silvestri, S., G. L. Wagner, et al., 2024, J. Adv. Model. Earth Syst., 16, e2023MS004130.
      !!----------------------------------------------------------------------
      INTEGER,  DIMENSION(2)                   , INTENT(in   ) ::   ktpuvw
      INTEGER                                  , INTENT(in   ) ::   kt              ! ocean time-step index
      INTEGER                                  , INTENT(in   ) ::   Kmm, Krhs       ! ocean time level indices
      INTEGER                                  , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                         , INTENT(in   ) ::   cdtype          ! =TRA or TRC (tracer indicator)
      INTEGER                                  , INTENT(in   ) ::   kjpt            ! number of tracers
      INTEGER                                  , INTENT(in   ) ::   kn_weno_h       ! order of the horizontal scheme (3 or 5)
      INTEGER                                  , INTENT(in   ) ::   kn_weno_v       ! order of the vertical   scheme (3 or 5)
      REAL(wp), DIMENSION(AB2D(ktpuvw),JPK    ), INTENT(in   ) ::   pU, pV, pW      ! 3 ocean volume flux components
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt,jpt), INTENT(inout) ::   pt              ! tracers and RHS of tracer equation
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      INTEGER  ::   ikm3, ikm2, ikp2 ! local integers
      LOGICAL  ::   ll_pos           ! upstream side of the face
      LOGICAL  ::   ll_v5            ! 5th order in the vertical
      REAL(wp) ::   ztra, zm3, zm5                       ! local scalars
      REAL(wp) ::   zt_m2, zt_m1, zt_0, zt_p1, zt_p2     !   -      -
      REAL(wp), DIMENSION(T2D(1))     ::   ztFu, ztFv   ! 2D workspace
      REAL(wp), DIMENSION(T2D(0),jpk) ::   ztFw         ! 3D workspace
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdx, ztrdy, ztrdz
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == kit000 )  THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'tra_adv_weno : WENO advection scheme on ', cdtype
            IF(lwp) WRITE(numout,*) '               horizontal order ', kn_weno_h, ' vertical order ', kn_weno_v
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
         ENDIF
         !
         l_trd = .FALSE.
         l_hst = .FALSE.
         l_ptr = .FALSE.
         IF( ( cdtype == 'TRA' .AND. l_trdtra ) .OR. ( cdtype == 'TRC' .AND. l_trdtrc ) )                      l_trd = .TRUE.
         IF(  l_diaptr .AND. cdtype == 'TRA' .AND. ( iom_use( 'sophtadv' )  .OR. iom_use( 'sopstadv'  ) )  )   l_ptr = .TRUE.
         IF(  l_diaar5 .AND. cdtype == 'TRA' .AND. ( iom_use("uadv_heattr") .OR. iom_use("vadv_heattr") .OR. &
            &                                        iom_use("uadv_salttr") .OR. iom_use("vadv_salttr") )  )   l_hst = .TRUE.
      ENDIF
      !
      IF( l_trd .OR. l_hst .OR. l_ptr )  THEN
         ALLOCATE( ztrdx(T2D(1),jpk), ztrdy(T2D(1),jpk), ztrdz(T2D(0),jpk) )
         ztrdx(:,:,jpk) = 0._wp   ;   ztrdy(:,:,jpk) = 0._wp
         ztrdz(:,:,1  ) = 0._wp   ;   ztrdz(:,:,jpk) = 0._wp
      ENDIF
      !
      ll_v5 = kn_weno_v == 5
      ztFw(:,:,jpk) = 0._wp
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         DO jk = 1, jpkm1
            !
            IF( kn_weno_h == 5 ) THEN
               DO_2D( 1, 0, 1, 0 )   !==  Horizontal advective fluxes: 5th order  ==!
                  !                               ! i-direction: upstream cell ji (pU>=0) or ji+1 (pU<0)
                  ll_pos = pU(ji,jj,jk) >= 0._wp
                  zt_m2 = MERGE( pt(ji-2,jj,jk,jn,Kmm), pt(ji+3,jj,jk,jn,Kmm), ll_pos )   ! 3rd upstream cell
                  zt_m1 = MERGE( pt(ji-1,jj,jk,jn,Kmm), pt(ji+2,jj,jk,jn,Kmm), ll_pos )   ! 2nd upstream cell
                  zt_0  = MERGE( pt(ji  ,jj,jk,jn,Kmm), pt(ji+1,jj,jk,jn,Kmm), ll_pos )   ! upstream cell
                  zt_p1 = MERGE( pt(ji+1,jj,jk,jn,Kmm), pt(ji  ,jj,jk,jn,Kmm), ll_pos )   ! downstream cell
                  zt_p2 = MERGE( pt(ji+2,jj,jk,jn,Kmm), pt(ji-1,jj,jk,jn,Kmm), ll_pos )   ! 2nd downstream cell
                  zm3 = MERGE( tmask(ji-1,jj,jk), tmask(ji+2,jj,jk), ll_pos )
                  zm5 = zm3 * MERGE( tmask(ji-2,jj,jk) * tmask(ji+2,jj,jk), tmask(ji+3,jj,jk) * tmask(ji-1,jj,jk), ll_pos )
                  ztFu(ji,jj) = pU(ji,jj,jk) * weno5_face( zt_m2, zt_m1, zt_0, zt_p1, zt_p2, zm5, zm3 )
                  !                               ! j-direction: upstream cell jj (pV>=0) or jj+1 (pV<0)
                  ll_pos = pV(ji,jj,jk) >= 0._wp
                  zt_m2 = MERGE( pt(ji,jj-2,jk,jn,Kmm), pt(ji,jj+3,jk,jn,Kmm), ll_pos )
                  zt_m1 = MERGE( pt(ji,jj-1,jk,jn,Kmm), pt(ji,jj+2,jk,jn,Kmm), ll_pos )
                  zt_0  = MERGE( pt(ji,jj  ,jk,jn,Kmm), pt(ji,jj+1,jk,jn,Kmm), ll_pos )
                  zt_p1 = MERGE( pt(ji,jj+1,jk,jn,Kmm), pt(ji,jj  ,jk,jn,Kmm), ll_pos )
                  zt_p2 = MERGE( pt(ji,jj+2,jk,jn,Kmm), pt(ji,jj-1,jk,jn,Kmm), ll_pos )
                  zm3 = MERGE( tmask(ji,jj-1,jk), tmask(ji,jj+2,jk), ll_pos )
                  zm5 = zm3 * MERGE( tmask(ji,jj-2,jk) * tmask(ji,jj+2,jk), tmask(ji,jj+3,jk) * tmask(ji,jj-1,jk), ll_pos )
                  ztFv(ji,jj) = pV(ji,jj,jk) * weno5_face( zt_m2, zt_m1, zt_0, zt_p1, zt_p2, zm5, zm3 )
               END_2D
            ELSE
               DO_2D( 1, 0, 1, 0 )   !==  Horizontal advective fluxes: 3rd order  ==!
                  !                               ! i-direction: upstream cell ji (pU>=0) or ji+1 (pU<0)
                  ll_pos = pU(ji,jj,jk) >= 0._wp
                  zt_m1 = MERGE( pt(ji-1,jj,jk,jn,Kmm), pt(ji+2,jj,jk,jn,Kmm), ll_pos )   ! 2nd upstream cell
                  zt_0  = MERGE( pt(ji  ,jj,jk,jn,Kmm), pt(ji+1,jj,jk,jn,Kmm), ll_pos )   ! upstream cell
                  zt_p1 = MERGE( pt(ji+1,jj,jk,jn,Kmm), pt(ji  ,jj,jk,jn,Kmm), ll_pos )   ! downstream cell
                  zm3 = MERGE( tmask(ji-1,jj,jk), tmask(ji+2,jj,jk), ll_pos )
                  ztFu(ji,jj) = pU(ji,jj,jk) * MERGE( weno3( zt_m1, zt_0, zt_p1 ), zt_0, zm3 > 0.5_wp )
                  !                               ! j-direction: upstream cell jj (pV>=0) or jj+1 (pV<0)
                  ll_pos = pV(ji,jj,jk) >= 0._wp
                  zt_m1 = MERGE( pt(ji,jj-1,jk,jn,Kmm), pt(ji,jj+2,jk,jn,Kmm), ll_pos )
                  zt_0  = MERGE( pt(ji,jj  ,jk,jn,Kmm), pt(ji,jj+1,jk,jn,Kmm), ll_pos )
                  zt_p1 = MERGE( pt(ji,jj+1,jk,jn,Kmm), pt(ji,jj  ,jk,jn,Kmm), ll_pos )
                  zm3 = MERGE( tmask(ji,jj-1,jk), tmask(ji,jj+2,jk), ll_pos )
                  ztFv(ji,jj) = pV(ji,jj,jk) * MERGE( weno3( zt_m1, zt_0, zt_p1 ), zt_0, zm3 > 0.5_wp )
               END_2D
            ENDIF
            !
            DO_2D( 0, 0, 0, 0 )   !==  add the horizontal advective trend  ==!
               ztra = - (  ( ztFu(ji,jj) - ztFu(ji-1,jj  ) )   &   ! add () for NP reproducibility
                  &      + ( ztFv(ji,jj) - ztFv(ji  ,jj-1) ) ) * r1_e1e2t(ji,jj)
               pt(ji,jj,jk,jn,Krhs) = pt(ji,jj,jk,jn,Krhs) + ztra * tmask(ji,jj,jk)
            END_2D
            !
            IF( l_trd .OR. l_hst .OR. l_ptr ) THEN   ! trend diagnostics // heat/salt transport
               DO_2D( 1, 0, 1, 0 )
                  ztrdx(ji,jj,jk) = ztFu(ji,jj)
                  ztrdy(ji,jj,jk) = ztFv(ji,jj)
               END_2D
            ENDIF
            !
         END DO
         !
         !                       !==  vertical advective fluxes  ==!   (k increases downward, pW>0 upward)
         !                                 ! out-of-column indices are clamped and masked
         DO_3D( 0, 0, 0, 0, 2, jpkm1 )
            ikm3 = MAX( jk-3, 1 )
            ikm2 = MAX( jk-2, 1 )
            ikp2 = MIN( jk+2, jpk )        ! tmask(jpk) = 0
            ll_pos = pW(ji,jj,jk) >= 0._wp        ! upstream cell jk (pW>=0) or jk-1 (pW<0)
            zt_m2 = MERGE( pt(ji,jj,ikp2,jn,Kmm), pt(ji,jj,ikm3,jn,Kmm), ll_pos )
            zt_m1 = MERGE( pt(ji,jj,jk+1,jn,Kmm), pt(ji,jj,ikm2,jn,Kmm), ll_pos )
            zt_0  = MERGE( pt(ji,jj,jk  ,jn,Kmm), pt(ji,jj,jk-1,jn,Kmm), ll_pos )
            zt_p1 = MERGE( pt(ji,jj,jk-1,jn,Kmm), pt(ji,jj,jk  ,jn,Kmm), ll_pos )
            zt_p2 = MERGE( pt(ji,jj,ikm2,jn,Kmm), pt(ji,jj,jk+1,jn,Kmm), ll_pos )
            zm3 = MERGE( tmask(ji,jj,jk+1), MERGE( tmask(ji,jj,ikm2), 0._wp, jk > 2 ), ll_pos )
            zm5 = zm3 * MERGE( tmask(ji,jj,ikp2) * MERGE( tmask(ji,jj,ikm2), 0._wp, jk > 2 ),   &
               &               MERGE( tmask(ji,jj,ikm3), 0._wp, jk > 3 ) * tmask(ji,jj,jk+1), ll_pos )
            zm5 = MERGE( zm5, 0._wp, ll_v5 )
            ztFw(ji,jj,jk) = pW(ji,jj,jk) * weno5_face( zt_m2, zt_m1, zt_0, zt_p1, zt_p2, zm5, zm3 ) * wmask(ji,jj,jk)
         END_3D
         !
         IF( lk_linssh ) THEN                ! top ocean value (only in linear free surface as ztFw has been w-masked)
            IF( ln_isfcav ) THEN                   ! top of the ice-shelf cavities and at the ocean surface
               DO_2D( 0, 0, 0, 0 )
                  ztFw(ji,jj, mikt(ji,jj) ) = pW(ji,jj,mikt(ji,jj)) * pt(ji,jj,mikt(ji,jj),jn,Kmm)   ! linear free surface
               END_2D
            ELSE                                   ! no cavities: only at the ocean surface
               DO_2D( 0, 0, 0, 0 )
                  ztFw(ji,jj,1) = pW(ji,jj,1) * pt(ji,jj,1,jn,Kmm)
               END_2D
            ENDIF
         ELSE
            DO_2D( 0, 0, 0, 0 )
               ztFw(ji,jj,1) = 0._wp
            END_2D
         ENDIF
         !
         DO_3D( 0, 0, 0, 0, 1, jpkm1 )   !==  add the vertical advective trend  ==!
            ztra = - ( ztFw(ji,jj,jk) - ztFw(ji,jj,jk+1) ) * r1_e1e2t(ji,jj)
            pt(ji,jj,jk,jn,Krhs) = pt(ji,jj,jk,jn,Krhs) + ztra * tmask(ji,jj,jk)
         END_3D
         !
         IF( l_trd ) THEN              ! trend diagnostics
            DO_3D( 0, 0, 0, 0, 2, jpkm1 )
               ztrdz(ji,jj,jk) = ztFw(ji,jj,jk)
            END_3D
            CALL trd_tra( kt, Kmm, Krhs, cdtype, jn, jptra_xad, ztrdx, pU, pt(:,:,:,jn,Kmm) )
            CALL trd_tra( kt, Kmm, Krhs, cdtype, jn, jptra_yad, ztrdy, pV, pt(:,:,:,jn,Kmm) )
            CALL trd_tra( kt, Kmm, Krhs, cdtype, jn, jptra_zad, ztrdz, pW, pt(:,:,:,jn,Kmm) )
         ENDIF
         IF( l_hst )   CALL dia_ar5_hst( jn, 'adv', ztrdx(:,:,:), ztrdy(:,:,:) ) ! heat/salt transport
         IF( l_ptr )   CALL dia_ptr_hst( jn, 'adv', ztrdy(:,:,:) )            ! "Poleward" transports
         !
      END DO
      !
      IF( l_trd .OR. l_hst .OR. l_ptr )   DEALLOCATE( ztrdx, ztrdy, ztrdz )
      !
   END SUBROUTINE tra_adv_weno_t


   !! * WENO reconstructions (weno5_face, weno5, weno3, ...)
#  include "weno_kernels.h90"

   !!======================================================================
END MODULE traadv_weno
