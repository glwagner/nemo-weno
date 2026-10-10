MODULE dynkeg
   !!======================================================================
   !!                       ***  MODULE  dynkeg  ***
   !! Ocean dynamics:  kinetic energy gradient trend
   !!======================================================================
   !! History :  1.0  !  1987-09  (P. Andrich, M.-A. Foujols)  Original code
   !!            7.0  !  1997-05  (G. Madec)  Split dynber into dynkeg and dynhpg
   !!  NEMO      1.0  !  2002-07  (G. Madec)  F90: Free form and module
   !!            3.6  !  2015-05  (N. Ducousso, G. Madec)  add Hollingsworth scheme as an option
   !!            4.5  !  2022-06  (S. Techene, G, Madec) refactorization to reduce local memory usage
   !!            5.x  !  2026-10  (G. Wagner)  WENO upwind kinetic energy gradient (dyn_keg_weno), following Oceananigans.jl
   !!----------------------------------------------------------------------
   
   !!----------------------------------------------------------------------
   !!   dyn_keg      : update the momentum trend with the horizontal tke
   !!   dyn_keg_weno : idem with an upwind (WENO) reconstruction of the gradient (vector invariant form)
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE trd_oce         ! trends: ocean variables
   USE trddyn          ! trend manager: dynamics
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE prtctl          ! Print control
   USE timing          ! Timing
   USE bdy_oce         ! ocean open boundary conditions

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_keg        ! routine called by stp2d.F90 and dynadv.F90
   PUBLIC   dyn_keg_weno   ! routine called by stp2d.F90 and dynadv.F90
   
   INTEGER, PARAMETER, PUBLIC  ::   nkeg_C2  = 0   !: 2nd order centered scheme (standard scheme)
   INTEGER, PARAMETER, PUBLIC  ::   nkeg_HW  = 1   !: Hollingsworth et al., QJRMS, 1983
   !
   REAL(wp) ::   r1_48 = 1._wp / 48._wp   !: =1/(4*2*6)
   
   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_keg( kt, kscheme, Kmm, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_keg  ***
      !!
      !! ** Purpose :   Compute the now momentum trend due to the horizontal
      !!      gradient of the horizontal kinetic energy and add it to the 
      !!      general momentum trend.
      !!
      !! ** Method  : * kscheme = nkeg_C2 : 2nd order centered scheme that 
      !!      conserve kinetic energy. Compute the now horizontal kinetic energy 
      !!         zhke = 1/2 [ mi-1( un^2 ) + mj-1( vn^2 ) ]
      !!              * kscheme = nkeg_HW : Hollingsworth correction following
      !!      Arakawa (2001). The now horizontal kinetic energy is given by:
      !!         zhke = 1/6 [ mi-1(  2 * un^2 + ((u(j+1)+u(j-1))/2)^2  )
      !!                    + mj-1(  2 * vn^2 + ((v(i+1)+v(i-1))/2)^2  ) ]
      !!      
      !!      Take its horizontal gradient and add it to the general momentum
      !!      trend.
      !!         u(rhs) = u(rhs) - 1/e1u di[ zhke ]
      !!         v(rhs) = v(rhs) - 1/e2v dj[ zhke ]
      !!
      !! ** Action : - Update the (puu(:,:,:,Krhs), pvv(:,:,:,Krhs)) with the hor. ke gradient trend
      !!             - send this trends to trd_dyn (l_trddyn=T) for post-processing
      !!
      !! ** References : Arakawa, A., International Geophysics 2001.
      !!                 Hollingsworth et al., Quart. J. Roy. Meteor. Soc., 1983.
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::   kt          ! ocean time-step index
      INTEGER                             , INTENT(in   ) ::   kscheme     ! =0/1   type of KEG scheme 
      INTEGER                             , INTENT(in   ) ::   Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::   puu, pvv    ! ocean velocities and RHS of momentum equation
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zu, zv       ! local scalars
#if ! defined key_PSYCLONE_2p5p0
      REAL(wp), DIMENSION(:,:  ) , ALLOCATABLE ::   zhke
#else
      REAL(wp), DIMENSION(T2D(1))              ::   zhke
#endif
      REAL(wp), DIMENSION(:,:,:) , ALLOCATABLE ::   zu_trd, zv_trd
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_keg')
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn_keg : kinetic energy gradient trend, scheme number=', kscheme
            IF(lwp) WRITE(numout,*) '~~~~~~~'
         ENDIF
      ENDIF
      !
      IF( l_trddyn ) THEN                       ! Save the input trends
         ALLOCATE( zu_trd(T2D(0),jpk), zv_trd(T2D(0),jpk) )
         zu_trd(:,:,:) = puu(T2D(0),:,Krhs)
         zv_trd(:,:,:) = pvv(T2D(0),:,Krhs)
      ENDIF
      !
      SELECT CASE ( kscheme )
      !
      CASE ( nkeg_C2 )                    !==  Standard scheme  ==!
#if ! defined key_PSYCLONE_2p5p0
         ALLOCATE( zhke(T2D(1)) )
#endif
         DO jk = 1, jpkm1
            DO_2D( 0, 1, 0, 1 )                 !* Horizontal kinetic energy at T-point
               zu =    puu(ji-1,jj  ,jk,Kmm) * puu(ji-1,jj  ,jk,Kmm)   &
                  &  + puu(ji  ,jj  ,jk,Kmm) * puu(ji  ,jj  ,jk,Kmm)
               zv =    pvv(ji  ,jj-1,jk,Kmm) * pvv(ji  ,jj-1,jk,Kmm)   &
                  &  + pvv(ji  ,jj  ,jk,Kmm) * pvv(ji  ,jj  ,jk,Kmm)
               zhke(ji,jj) = 0.25_wp * ( zv + zu )
            END_2D
            !
            DO_2D( 0, 0, 0, 0 )                  !* grad( KE ) added to the general momentum trends
               puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - ( zhke(ji+1,jj  ) - zhke(ji,jj) ) * r1_e1u(ji,jj)
               pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - ( zhke(ji  ,jj+1) - zhke(ji,jj) ) * r1_e2v(ji,jj)
            END_2D
         END DO
#if ! defined key_PSYCLONE_2p5p0
         DEALLOCATE( zhke )
#endif
         !
      CASE ( nkeg_HW )                           !* Hollingsworth scheme
#if ! defined key_PSYCLONE_2p5p0
         ALLOCATE( zhke(T2D(1)) )
#endif
         DO jk = 1, jpkm1
            DO_2D( 0, 1, 0, 1 )
               ! round brackets added to fix the order of floating point operations
               ! needed to ensure halo 1 - halo 2 compatibility
               zu =   (   puu(ji-1,jj  ,jk,Kmm) * puu(ji-1,jj  ,jk,Kmm)               &
                  &     + puu(ji  ,jj  ,jk,Kmm) * puu(ji  ,jj  ,jk,Kmm)   ) * 8._wp   &
                  & + ( ( puu(ji-1,jj-1,jk,Kmm) + puu(ji-1,jj+1,jk,Kmm) ) * ( puu(ji-1,jj-1,jk,Kmm) + puu(ji-1,jj+1,jk,Kmm) )   &
                  & +   ( puu(ji  ,jj-1,jk,Kmm) + puu(ji  ,jj+1,jk,Kmm) ) * ( puu(ji  ,jj-1,jk,Kmm) + puu(ji  ,jj+1,jk,Kmm) )   &
                  &   )                                                    ! bracket for halo 1 - halo 2 compatibility
               zv =   (   pvv(ji  ,jj-1,jk,Kmm) * pvv(ji  ,jj-1,jk,Kmm)               &
                  &     + pvv(ji  ,jj  ,jk,Kmm) * pvv(ji  ,jj  ,jk,Kmm)   ) * 8._wp   &
                  & + ( ( pvv(ji-1,jj-1,jk,Kmm) + pvv(ji+1,jj-1,jk,Kmm) ) * ( pvv(ji-1,jj-1,jk,Kmm) + pvv(ji+1,jj-1,jk,Kmm) )   &
                  & +   ( pvv(ji-1,jj  ,jk,Kmm) + pvv(ji+1,jj  ,jk,Kmm) ) * ( pvv(ji-1,jj  ,jk,Kmm) + pvv(ji+1,jj  ,jk,Kmm) )   &
                  &   )                                                    ! bracket for halo 1 - halo 2 compatibility
               zhke(ji,jj) = r1_48 * ( zv + zu )
            END_2D
            !
            DO_2D( 0, 0, 0, 0 )                  !* grad( KE ) added to the general momentum trends
               puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - ( zhke(ji+1,jj  ) - zhke(ji,jj) ) * r1_e1u(ji,jj)
               pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - ( zhke(ji  ,jj+1) - zhke(ji,jj) ) * r1_e2v(ji,jj)
            END_2D
         END DO
#if ! defined key_PSYCLONE_2p5p0
         DEALLOCATE( zhke )
#endif
         !
      END SELECT
      !
      IF( l_trddyn ) THEN                        ! save the Kinetic Energy trends for diagnostic
         zu_trd(:,:,:) = puu(T2D(0),:,Krhs) - zu_trd(:,:,:)
         zv_trd(:,:,:) = pvv(T2D(0),:,Krhs) - zv_trd(:,:,:)
         CALL trd_dyn( zu_trd, zv_trd, jpdyn_keg, kt, Kmm )
         DEALLOCATE( zu_trd, zv_trd )
      ENDIF
      !
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' keg  - Ua: ', mask1=umask,   &
         &                                  tab3d_2=pvv(:,:,:,Krhs), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_keg')
      !
   END SUBROUTINE dyn_keg


   SUBROUTINE dyn_keg_weno( kt, korder, Kmm, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_keg_weno  ***
      !!
      !! ** Purpose :   Compute the now momentum trend due to the horizontal
      !!      gradient of the horizontal kinetic energy with an upwind (WENO)
      !!      reconstruction, and add it to the general momentum trend.
      !!
      !! ** Method  :   This is the kinetic energy gradient of Oceananigans
      !!      WENOVectorInvariant with its default OnlySelfUpwinding treatment
      !!      (Silvestri et al. 2024). At a u-point the gradient of u^2/2 is
      !!      reconstructed from its differences at the surrounding t-points,
      !!      upstream of u (WENO of order korder, with weights computed from the
      !!      smoothness of u averaged at the same t-points), while the gradient
      !!      of v^2/2 along i is interpolated from the f-points along j with a
      !!      centered scheme (4th order for korder = 5, 2nd order otherwise):
      !!         u(rhs) = u(rhs) - 1/e1u [ WENO_i( di[u^2/2] ) + mj( di[v^2/2] ) ]
      !!         v(rhs) = v(rhs) - 1/e2v [ WENO_j( dj[v^2/2] ) + mi( dj[u^2/2] ) ]
      !!      Near land the order of the WENO reconstruction is reduced
      !!      (5 -> 3 -> 1) so that only ocean t-points are used, and the
      !!      centered interpolation falls back to 2nd order if one of its
      !!      f-points has no ocean t-point around.
      !!        The reconstruction needs (korder+1)/2 points of halo.
      !!
      !! ** Action : - Update the (puu(:,:,:,Krhs), pvv(:,:,:,Krhs)) with the hor. ke gradient trend
      !!             - send this trends to trd_dyn (l_trddyn=T) for post-processing
      !!
      !! References : Silvestri, S., G. L. Wagner, et al., 2024, J. Adv. Model. Earth Syst., 16, e2023MS004130.
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::   kt          ! ocean time-step index
      INTEGER                             , INTENT(in   ) ::   korder      ! order of the WENO reconstruction (=3/5)
      INTEGER                             , INTENT(in   ) ::   Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::   puu, pvv    ! ocean velocities and RHS of momentum equation
      !
      INTEGER  ::   ji, jj, jk, ib    ! dummy loop indices, half width of the stencil
      REAL(wp) ::   zu, zv, zku, zkv  ! local scalars
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zdxu2, zdyv2   ! di[u^2/2] and dj[v^2/2] at t-points
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zus  , zvs     ! u averaged along i and v averaged along j at t-points
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zdxv2, zdyu2   ! di[v^2/2] and dj[u^2/2] at f-points
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zmf            ! =1 at active f-points
      REAL(wp), DIMENSION(:,:,:) , ALLOCATABLE ::   zu_trd, zv_trd
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_keg_weno')
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn_keg_weno : kinetic energy gradient trend, WENO upwind reconstruction of order ', korder
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
         ENDIF
      ENDIF
      !
      IF( l_trddyn ) THEN                       ! Save the input trends
         ALLOCATE( zu_trd(T2D(0),jpk), zv_trd(T2D(0),jpk) )
         zu_trd(:,:,:) = puu(T2D(0),:,Krhs)
         zv_trd(:,:,:) = pvv(T2D(0),:,Krhs)
      ENDIF
      !
      ib = ( korder + 1 ) / 2                           ! t-points needed on each side of a velocity point
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO_2D( ib-1, ib, ib-1, ib )                   !==  t-point quantities  ==!
            zdxu2(ji,jj) = 0.5_wp * ( puu(ji,jj,jk,Kmm) * puu(ji,jj,jk,Kmm) - puu(ji-1,jj,jk,Kmm) * puu(ji-1,jj,jk,Kmm) )
            zdyv2(ji,jj) = 0.5_wp * ( pvv(ji,jj,jk,Kmm) * pvv(ji,jj,jk,Kmm) - pvv(ji,jj-1,jk,Kmm) * pvv(ji,jj-1,jk,Kmm) )
            zus  (ji,jj) = 0.5_wp * ( puu(ji,jj,jk,Kmm) + puu(ji-1,jj  ,jk,Kmm) )
            zvs  (ji,jj) = 0.5_wp * ( pvv(ji,jj,jk,Kmm) + pvv(ji  ,jj-1,jk,Kmm) )
         END_2D
         DO_2D( 2, 1, 2, 1 )                           !==  f-point quantities  ==!
            zdxv2(ji,jj) = 0.5_wp * ( pvv(ji+1,jj,jk,Kmm) * pvv(ji+1,jj,jk,Kmm) - pvv(ji,jj,jk,Kmm) * pvv(ji,jj,jk,Kmm) )
            zdyu2(ji,jj) = 0.5_wp * ( puu(ji,jj+1,jk,Kmm) * puu(ji,jj+1,jk,Kmm) - puu(ji,jj,jk,Kmm) * puu(ji,jj,jk,Kmm) )
            zmf  (ji,jj) = MAX( tmask(ji,jj,jk), tmask(ji+1,jj,jk), tmask(ji,jj+1,jk), tmask(ji+1,jj+1,jk) )
         END_2D
         !
         IF( korder == 5 ) THEN                        !==  grad( KE ) added to the general momentum trends  ==!
            DO_2D( 0, 0, 0, 0 )
               zu  = puu(ji,jj,jk,Kmm)
               zku = weno_up5( zdxu2(ji-2:ji+3,jj), zus(ji-2:ji+3,jj), tmask(ji-2:ji+3,jj,jk), zu )
               zkv = cen4    ( zdxv2(ji,jj-2:jj+1), zmf(ji,jj-2:jj+1) )
               puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - ( zku + zkv ) * r1_e1u(ji,jj)
               zv  = pvv(ji,jj,jk,Kmm)
               zkv = weno_up5( zdyv2(ji,jj-2:jj+3), zvs(ji,jj-2:jj+3), tmask(ji,jj-2:jj+3,jk), zv )
               zku = cen4    ( zdyu2(ji-2:ji+1,jj), zmf(ji-2:ji+1,jj) )
               pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - ( zkv + zku ) * r1_e2v(ji,jj)
            END_2D
         ELSE
            DO_2D( 0, 0, 0, 0 )
               zu  = puu(ji,jj,jk,Kmm)
               zku = weno_up3( zdxu2(ji-1:ji+2,jj), zus(ji-1:ji+2,jj), tmask(ji-1:ji+2,jj,jk), zu )
               zkv = 0.5_wp * ( zdxv2(ji,jj-1) + zdxv2(ji,jj) )
               puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - ( zku + zkv ) * r1_e1u(ji,jj)
               zv  = pvv(ji,jj,jk,Kmm)
               zkv = weno_up3( zdyv2(ji,jj-1:jj+2), zvs(ji,jj-1:jj+2), tmask(ji,jj-1:jj+2,jk), zv )
               zku = 0.5_wp * ( zdyu2(ji-1,jj) + zdyu2(ji,jj) )
               pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - ( zkv + zku ) * r1_e2v(ji,jj)
            END_2D
         ENDIF
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      IF( l_trddyn ) THEN                        ! save the Kinetic Energy trends for diagnostic
         zu_trd(:,:,:) = puu(T2D(0),:,Krhs) - zu_trd(:,:,:)
         zv_trd(:,:,:) = pvv(T2D(0),:,Krhs) - zv_trd(:,:,:)
         CALL trd_dyn( zu_trd, zv_trd, jpdyn_keg, kt, Kmm )
         DEALLOCATE( zu_trd, zv_trd )
      ENDIF
      !
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' keg  - Ua: ', mask1=umask,   &
         &                                  tab3d_2=pvv(:,:,:,Krhs), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_keg_weno')
      !
   END SUBROUTINE dyn_keg_weno

   !!----------------------------------------------------------------------
   !! * WENO reconstructions (weno_up5, weno_up3, cen4, ...)
   !!----------------------------------------------------------------------
#  include "weno_kernels.h90"

   !!======================================================================
END MODULE dynkeg
