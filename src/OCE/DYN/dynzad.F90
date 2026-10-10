MODULE dynzad
   !!======================================================================
   !!                       ***  MODULE  dynzad  ***
   !! Ocean dynamics : vertical advection trend
   !!======================================================================
   !! History :  OPA  !  1991-01  (G. Madec) Original code
   !!   NEMO     0.5  !  2002-07  (G. Madec) Free form, F90
   !!            4.5  !  2021-01  (S. Techene, G. Madec) memory optimization
   !!            5.x  !  2026-10  (G. Wagner)  WENO upwind vertical advection (dyn_zad_weno), following Oceananigans.jl
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zad       : vertical advection momentum trend
   !!   dyn_zad_weno  : idem with upwind (WENO) reconstructions (vector invariant form)
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! surface boundary condition: ocean
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   USE sbcwave, ONLY: wsd   ! Surface Waves (add vertical Stokes-drift)
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zad       ! routine called by stp2d.F90 and dynadv.F90
   PUBLIC   dyn_zad_weno  ! routine called by stp2d.F90 and dynadv.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_zad ( kt, Kmm, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dynzad  ***
      !!
      !! ** Purpose :   Compute the now vertical momentum advection trend and
      !!      add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The now vertical advection of momentum is given by:
      !!         w dz(u) = u(rhs) + 1/(e1e2u*e3u) mk+1[ mi(e1e2t*ww) dk(u) ]
      !!         w dz(v) = v(rhs) + 1/(e1e2v*e3v) mk+1[ mj(e1e2t*ww) dk(v) ]
      !!      Add this trend to the general trend (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)):
      !!         (u(rhs),v(rhs)) = (u(rhs),v(rhs)) + w dz(u,v)
      !!
      !! ** Action  : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the vert. momentum adv. trends
      !!              - Send the trends to trddyn for diagnostics (l_trddyn=T)
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::    kt, Kmm, Krhs   ! ocean time-step & time-level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::    puu, pvv        ! ocean velocities and RHS of momentum equation
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zWf, zWfi, zzWfu, zzWdzU   ! local scalars
      REAL(wp) ::        zWfj, zzWfv, zzWdzV   !   -      -
      REAL(wp), DIMENSION(T2D(0))             ::   zWdzU, zWdzV   ! 2D inner workspace
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdu, ztrdv   ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_zad')
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn_zad : 2nd order vertical advection scheme'
         ENDIF
      ENDIF
      !
      IF( l_trddyn )   THEN           ! Save puu(:,:,:,Krhs) and pvv(:,:,:,Krhs) trends
         ALLOCATE( ztrdu(T2D(0),jpk), ztrdv(T2D(0),jpk) )
         ztrdu(:,:,:) = puu(T2D(0),:,Krhs)
         ztrdv(:,:,:) = pvv(T2D(0),:,Krhs)
      ENDIF
      !
      !                                !==  vertical momentum advection ==!   at u- and v-points
      !
      zWdzU(T2D(0)) = 0._wp                  ! set surface (jk=1) vertical advection to zero
      zWdzV(T2D(0)) = 0._wp
      !        
      DO_3D( 0, 0, 0, 0 , 1, jpk-2 )   !=  surface to jpk-2 vertical advection
         !                                ! vertical transport at jk+1 uw/vw-level (x2): 2*mi/j[e1e2t*(We)]
         IF( ln_vortex_force ) THEN             ! We = ww+wsd
            zWf  = e1e2t(ji  ,jj  ) * ( ww(ji  ,jj  ,jk+1) + wsd(ji  ,jj  ,jk+1) )
            zWfi = e1e2t(ji+1,jj  ) * ( ww(ji+1,jj  ,jk+1) + wsd(ji+1,jj  ,jk+1) )
            zWfj = e1e2t(ji  ,jj+1) * ( ww(ji  ,jj+1,jk+1) + wsd(ji  ,jj+1,jk+1) )
         ELSE                                   ! We = ww
            zWf  = e1e2t(ji  ,jj  ) *   ww(ji  ,jj  ,jk+1)
            zWfi = e1e2t(ji+1,jj  ) *   ww(ji+1,jj  ,jk+1)
            zWfj = e1e2t(ji  ,jj+1) *   ww(ji  ,jj+1,jk+1)
         ENDIF
         zzWfu = zWfi + zWf                     ! averaging at uw- and vw-points (x2)
         zzWfv = zWfj + zWf 
         !                                ! vertical advection at jk+1 uw-level (x4): zzWfu/v*dk+1[u/v] 
         zzWdzU = zzWfu * ( puu(ji,jj,jk,Kmm) - puu(ji,jj,jk+1,Kmm) )
         zzWdzV = zzWfv * ( pvv(ji,jj,jk,Kmm) - pvv(ji,jj,jk+1,Kmm) )
         !
         !                                ! vertical advection  at jk u/v-level 
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - 0.25_wp * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)   &
              &                                            * ( zWdzU(ji,jj) + zzWdzU )
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - 0.25_wp * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)   &
              &                                            * ( zWdzV(ji,jj) + zzWdzV )
         !
         zWdzU(ji,jj) = zzWdzU            ! save for next level computation 
         zWdzV(ji,jj) = zzWdzV
      END_3D
      !
      jk = jpkm1
      DO_2D( 0, 0, 0, 0 )           !=  bottom vertical advection  at jpkm1
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - 0.25_wp * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)   &
            &                                              * zWdzU(ji,jj)
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - 0.25_wp * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)   &
            &                                              * zWdzV(ji,jj)
      END_2D
      !
      IF( l_trddyn ) THEN             ! save the vertical advection trends for diagnostic
         ztrdu(:,:,:) = puu(T2D(0),:,Krhs) - ztrdu(:,:,:)
         ztrdv(:,:,:) = pvv(T2D(0),:,Krhs) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_zad, kt, Kmm )
         DEALLOCATE( ztrdu, ztrdv )
      ENDIF
      !                               ! Control print
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' zad  - Ua: ', mask1=umask,   &
         &                                  tab3d_2=pvv(:,:,:,Krhs), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_zad')
      !
   END SUBROUTINE dyn_zad


   SUBROUTINE dyn_zad_weno( kt, korder, Kmm, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zad_weno  ***
      !!
      !! ** Purpose :   Compute the now vertical momentum advection trend with
      !!      upwind (WENO) reconstructions and add it to the general trend of
      !!      momentum equation.
      !!
      !! ** Method  :   This is the vertical advection of Oceananigans
      !!      WENOVectorInvariant with its default OnlySelfUpwinding treatment
      !!      (Silvestri et al. 2024). Using the continuity equation, the vertical
      !!      advection w dk[u] is split into a flux form vertical advection and
      !!      a divergence flux:
      !!         w dz(u) = u(rhs) - 1/(e1e2u*e3u) [ dk( Wu * u_k ) + u ( di[U] + dj[V] + e1e2t dt(e3t) ) ]
      !!      where Wu is the vertical transport e1e2t*ww interpolated at the
      !!      uw-point (4th order centered for korder = 5, 2nd order otherwise),
      !!      u_k is u reconstructed at the uw-point upstream of Wu (WENO of order
      !!      korder with the usual smoothness indicators) and U = e2u*e3u*u,
      !!      V = e1v*e3v*v. In the divergence flux, di[U] is reconstructed at
      !!      the u-point upstream of u (WENO of order korder, with weights
      !!      computed from the smoothness of the horizontal divergence
      !!      di[U] + dj[V]), while the rest is interpolated with the centered
      !!      scheme. The latter, dj[V] + e1e2t dt(e3t), is obtained from the
      !!      discrete continuity equation as  - e1e2t dk[ww] - di[U]  so that it
      !!      is consistent with the vertical velocity ww whatever the vertical
      !!      coordinate (the grid motion term of Oceananigans is the dt(e3t)
      !!      term of the NEMO continuity equation, which is included in ww by wzv).
      !!      The same holds for the v-component. The vertical fluxes vanish at
      !!      the bottom; the surface flux is Wu * u(jk=1) (it cancels the surface
      !!      vertical transport of the divergence flux).
      !!        Near land the order of the reconstructions is reduced (5 -> 3 -> 1
      !!      along i and j, and in the vertical) so that only ocean points are
      !!      used. The vertical velocity ww is needed over (korder+1)/2 - 1
      !!      points of halo (lbc_lnk done by the caller).
      !!
      !! ** Action  : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the vert. momentum adv. trends
      !!              - Send the trends to trddyn for diagnostics (l_trddyn=T)
      !!
      !! References : Silvestri, S., G. L. Wagner, et al., 2024, J. Adv. Model. Earth Syst., 16, e2023MS004130.
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::   kt, Kmm, Krhs   ! ocean time-step & time-level indices
      INTEGER                             , INTENT(in   ) ::   korder          ! order of the WENO reconstructions (=3/5)
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::   puu, pvv        ! ocean velocities and RHS of momentum equation
      !
      INTEGER  ::   ji, jj, jk, ib         ! dummy loop indices, half width of the stencil
      INTEGER  ::   ik1, ik2, ik3, ik4, ik5, ik6   ! vertical indices of the stencil
      REAL(wp) ::   zu, zv, zwu, zwv, zdu, zdv   ! local scalars
      REAL(wp), DIMENSION(6) ::   zst, zsm              ! vertical stencil and its mask
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zdxU, zdyV   ! di[U] and dj[V] at t-points
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zdiv         ! horizontal divergence di[U] + dj[V] at t-points (smoothness)
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zcru, zcrv   ! - e1e2t dk[ww] - di[U]  and  - e1e2t dk[ww] - dj[V] at t-points
      REAL(wp), DIMENSION(T2D(nn_hls)) ::   zwt          ! vertical transport e1e2t*ww at t-points
      REAL(wp), DIMENSION(T2D(0))      ::   zfut, zfvt   ! vertical fluxes at the top of the level
      REAL(wp), DIMENSION(T2D(0))      ::   zfub, zfvb   ! vertical fluxes at the bottom of the level
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdu, ztrdv   ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_zad_weno')
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn_zad_weno : vertical advection scheme, WENO upwind reconstructions of order ', korder
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
         ENDIF
      ENDIF
      !
      IF( l_trddyn )   THEN           ! Save puu(:,:,:,Krhs) and pvv(:,:,:,Krhs) trends
         ALLOCATE( ztrdu(T2D(0),jpk), ztrdv(T2D(0),jpk) )
         ztrdu(:,:,:) = puu(T2D(0),:,Krhs)
         ztrdv(:,:,:) = pvv(T2D(0),:,Krhs)
      ENDIF
      !
      ib = ( korder + 1 ) / 2                          ! t-points needed on each side of a velocity point
      !
      !                                !==  surface vertical fluxes  ==!   Wu * u(jk=1)  (1st order)
      DO_2D( ib-1, ib, ib-1, ib )
         zwt(ji,jj) = e1e2t(ji,jj) * ww(ji,jj,1)
      END_2D
      DO_2D( 0, 0, 0, 0 )
         zfut(ji,jj) = w_at_u( zwt, ji, jj, 1 ) * puu(ji,jj,1,Kmm) * wumask(ji,jj,1)
         zfvt(ji,jj) = w_at_v( zwt, ji, jj, 1 ) * pvv(ji,jj,1,Kmm) * wvmask(ji,jj,1)
      END_2D
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !                                   !==  vertical fluxes at the bottom of the level (jk+1 w-level)  ==!
         IF( jk < jpkm1 ) THEN
            DO_2D( ib-1, ib, ib-1, ib )
               zwt(ji,jj) = e1e2t(ji,jj) * ww(ji,jj,jk+1)
            END_2D
            ik1 = MIN( jk+3, jpk )   ;   ik2 = MIN( jk+2, jpk )   ;   ik3 = jk+1    ! stencil ordered upward, the jk+1 w-level
            ik4 = jk                 ;   ik5 = MAX( jk-1, 1 )     ;   ik6 = MAX( jk-2, 1 )   ! lying between ik3 and ik4
            DO_2D( 0, 0, 0, 0 )
               zsm(1) = MERGE( umask(ji,jj,ik1), 0._wp, jk+3 <= jpkm1 )   ;   zsm(2) = MERGE( umask(ji,jj,ik2), 0._wp, jk+2 <= jpkm1 )
               zsm(3) = umask(ji,jj,ik3)                                  ;   zsm(4) = umask(ji,jj,ik4)
               zsm(5) = MERGE( umask(ji,jj,ik5), 0._wp, jk-1 >= 1 )       ;   zsm(6) = MERGE( umask(ji,jj,ik6), 0._wp, jk-2 >= 1 )
               zst(1) = puu(ji,jj,ik1,Kmm)   ;   zst(2) = puu(ji,jj,ik2,Kmm)   ;   zst(3) = puu(ji,jj,ik3,Kmm)
               zst(4) = puu(ji,jj,ik4,Kmm)   ;   zst(5) = puu(ji,jj,ik5,Kmm)   ;   zst(6) = puu(ji,jj,ik6,Kmm)
               zwu = w_at_u( zwt, ji, jj, jk+1 )
               IF( korder == 5 ) THEN   ;   zu = weno_up5( zst     , zst     , zsm     , zwu )
               ELSE                     ;   zu = weno_up3( zst(2:5), zst(2:5), zsm(2:5), zwu )
               ENDIF
               zfub(ji,jj) = zwu * zu * wumask(ji,jj,jk+1)
               !
               zsm(1) = MERGE( vmask(ji,jj,ik1), 0._wp, jk+3 <= jpkm1 )   ;   zsm(2) = MERGE( vmask(ji,jj,ik2), 0._wp, jk+2 <= jpkm1 )
               zsm(3) = vmask(ji,jj,ik3)                                  ;   zsm(4) = vmask(ji,jj,ik4)
               zsm(5) = MERGE( vmask(ji,jj,ik5), 0._wp, jk-1 >= 1 )       ;   zsm(6) = MERGE( vmask(ji,jj,ik6), 0._wp, jk-2 >= 1 )
               zst(1) = pvv(ji,jj,ik1,Kmm)   ;   zst(2) = pvv(ji,jj,ik2,Kmm)   ;   zst(3) = pvv(ji,jj,ik3,Kmm)
               zst(4) = pvv(ji,jj,ik4,Kmm)   ;   zst(5) = pvv(ji,jj,ik5,Kmm)   ;   zst(6) = pvv(ji,jj,ik6,Kmm)
               zwv = w_at_v( zwt, ji, jj, jk+1 )
               IF( korder == 5 ) THEN   ;   zv = weno_up5( zst     , zst     , zsm     , zwv )
               ELSE                     ;   zv = weno_up3( zst(2:5), zst(2:5), zsm(2:5), zwv )
               ENDIF
               zfvb(ji,jj) = zwv * zv * wvmask(ji,jj,jk+1)
            END_2D
         ELSE                                !* no flux through the bottom (ww(jpk) = 0)
            zfub(:,:) = 0._wp
            zfvb(:,:) = 0._wp
         ENDIF
         !
         !                                   !==  divergence flux  ==!
         DO_2D( ib-1, ib, ib-1, ib )                   ! t-point quantities
            zdxU(ji,jj) = e2u(ji,jj) * e3u(ji,jj,jk,Kmm) * puu(ji,jj,jk,Kmm) - e2u(ji-1,jj) * e3u(ji-1,jj,jk,Kmm) * puu(ji-1,jj,jk,Kmm)
            zdyV(ji,jj) = e1v(ji,jj) * e3v(ji,jj,jk,Kmm) * pvv(ji,jj,jk,Kmm) - e1v(ji,jj-1) * e3v(ji,jj-1,jk,Kmm) * pvv(ji,jj-1,jk,Kmm)
            zdiv(ji,jj) = zdxU(ji,jj) + zdyV(ji,jj)
            zcru(ji,jj) = e1e2t(ji,jj) * ( ww(ji,jj,jk+1) - ww(ji,jj,jk) ) - zdxU(ji,jj)   ! = dj[V] + e1e2t dt(e3t)
            zcrv(ji,jj) = e1e2t(ji,jj) * ( ww(ji,jj,jk+1) - ww(ji,jj,jk) ) - zdyV(ji,jj)   ! = di[U] + e1e2t dt(e3t)
         END_2D
         !
         IF( korder == 5 ) THEN                        !==  trends  ==!
            DO_2D( 0, 0, 0, 0 )
               zu  = puu(ji,jj,jk,Kmm)
               zdu = weno_up5( zdxU(ji-2:ji+3,jj), zdiv(ji-2:ji+3,jj), tmask(ji-2:ji+3,jj,jk), zu )
               zdv = cen4    ( zcru(ji-1:ji+2,jj), tmask(ji-1:ji+2,jj,jk) )
               puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - (  ( zfut(ji,jj) - zfub(ji,jj) ) + zu * ( zdv + zdu )  )   &
                  &                                      * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)
               zv  = pvv(ji,jj,jk,Kmm)
               zdv = weno_up5( zdyV(ji,jj-2:jj+3), zdiv(ji,jj-2:jj+3), tmask(ji,jj-2:jj+3,jk), zv )
               zdu = cen4    ( zcrv(ji,jj-1:jj+2), tmask(ji,jj-1:jj+2,jk) )
               pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - (  ( zfvt(ji,jj) - zfvb(ji,jj) ) + zv * ( zdu + zdv )  )   &
                  &                                      * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)
            END_2D
         ELSE
            DO_2D( 0, 0, 0, 0 )
               zu  = puu(ji,jj,jk,Kmm)
               zdu = weno_up3( zdxU(ji-1:ji+2,jj), zdiv(ji-1:ji+2,jj), tmask(ji-1:ji+2,jj,jk), zu )
               zdv = 0.5_wp * ( zcru(ji,jj) + zcru(ji+1,jj) )
               puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - (  ( zfut(ji,jj) - zfub(ji,jj) ) + zu * ( zdv + zdu )  )   &
                  &                                      * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)
               zv  = pvv(ji,jj,jk,Kmm)
               zdv = weno_up3( zdyV(ji,jj-1:jj+2), zdiv(ji,jj-1:jj+2), tmask(ji,jj-1:jj+2,jk), zv )
               zdu = 0.5_wp * ( zcrv(ji,jj) + zcrv(ji,jj+1) )
               pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - (  ( zfvt(ji,jj) - zfvb(ji,jj) ) + zv * ( zdu + zdv )  )   &
                  &                                      * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)
            END_2D
         ENDIF
         !
         DO_2D( 0, 0, 0, 0 )                           ! bottom fluxes become the top fluxes of the next level
            zfut(ji,jj) = zfub(ji,jj)
            zfvt(ji,jj) = zfvb(ji,jj)
         END_2D
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      IF( l_trddyn ) THEN             ! save the vertical advection trends for diagnostic
         ztrdu(:,:,:) = puu(T2D(0),:,Krhs) - ztrdu(:,:,:)
         ztrdv(:,:,:) = pvv(T2D(0),:,Krhs) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_zad, kt, Kmm )
         DEALLOCATE( ztrdu, ztrdv )
      ENDIF
      !                               ! Control print
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' zad  - Ua: ', mask1=umask,   &
         &                                  tab3d_2=pvv(:,:,:,Krhs), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_zad_weno')
      !
   CONTAINS
      !
      PURE FUNCTION w_at_u( pwt, ki, kj, kk ) RESULT( pw )
         !! vertical transport at the uw-point (ki,kj) of the kk w-level, interpolated from the t-points
         REAL(wp), DIMENSION(T2D(nn_hls)), INTENT(in) ::   pwt          ! vertical transport at t-points
         INTEGER                         , INTENT(in) ::   ki, kj, kk   ! indices
         REAL(wp)                                     ::   pw
         IF( korder == 5 ) THEN   ;   pw = cen4( pwt(ki-1:ki+2,kj), tmask(ki-1:ki+2,kj,kk) )
         ELSE                     ;   pw = 0.5_wp * ( pwt(ki,kj) + pwt(ki+1,kj) )
         ENDIF
      END FUNCTION w_at_u
      !
      PURE FUNCTION w_at_v( pwt, ki, kj, kk ) RESULT( pw )
         !! vertical transport at the vw-point (ki,kj) of the kk w-level, interpolated from the t-points
         REAL(wp), DIMENSION(T2D(nn_hls)), INTENT(in) ::   pwt          ! vertical transport at t-points
         INTEGER                         , INTENT(in) ::   ki, kj, kk   ! indices
         REAL(wp)                                     ::   pw
         IF( korder == 5 ) THEN   ;   pw = cen4( pwt(ki,kj-1:kj+2), tmask(ki,kj-1:kj+2,kk) )
         ELSE                     ;   pw = 0.5_wp * ( pwt(ki,kj) + pwt(ki,kj+1) )
         ENDIF
      END FUNCTION w_at_v
      !
   END SUBROUTINE dyn_zad_weno

   !!----------------------------------------------------------------------
   !! * WENO reconstructions (weno_up5, weno_up3, cen4, ...)
   !!----------------------------------------------------------------------
#  include "weno_kernels.h90"

   !!======================================================================
END MODULE dynzad
