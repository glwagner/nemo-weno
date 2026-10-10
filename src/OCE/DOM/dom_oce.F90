MODULE dom_oce
   !!======================================================================
   !!                       ***  MODULE dom_oce  ***
   !! ** Purpose :   Define in memory all the ocean space domain variables
   !!======================================================================
   !! History :  1.0  ! 2005-10  (A. Beckmann, G. Madec)  reactivate s-coordinate
   !!            3.3  ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !!            3.4  ! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!            3.5  ! 2012     (S. Mocavero, I. Epicoco) Add arrays associated
   !!                             to the optimization of BDY communications
   !!            3.7  ! 2015-11  (G. Madec) introduce surface and scale factor ratio
   !!             -   ! 2015-11  (G. Madec, A. Coward)  time varying zgr by default
   !!            4.1  ! 2019-08  (A. Coward, D. Storkey) rename prognostic variables in preparation for new time scheme.
   !!            4.x  ! 2020-02  (G. Madec, S. Techene) introduce ssh to h0 ratio
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   Agrif_Root    : dummy function used when lk_agrif=F
   !!   Agrif_Fixed   : dummy function used when lk_agrif=F
   !!   Agrif_CFixed  : dummy function used when lk_agrif=F
   !!   dom_oce_alloc : dynamical allocation of dom_oce arrays
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters

   IMPLICIT NONE
   PUBLIC             ! allows the acces to par_oce when dom_oce is used (exception to coding rules)

   PUBLIC dom_oce_alloc    ! Called from nemogcm.F90
   PUBLIC dom_oce_dealloc  ! Called from nemogcm.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! time & space domain namelist
   !! ----------------------------
   LOGICAL , PUBLIC ::   ln_meshmask    !: =T  create a mesh-mask file (mesh_mask.nc)
   REAL(wp), PUBLIC ::   rn_Dt          !: time step for the dynamics and tracer
   REAL(wp), PUBLIC ::   rDt_e          !: Barotropic time step
   LOGICAL , PUBLIC ::   ln_c1d         !: =T  single column domain (1x1 pt)
   LOGICAL , PUBLIC ::   ln_shuman      !: =T  shuman averaging

   !                                          !**  time level indices  **!
   INTEGER, PUBLIC ::   Nbb, Nnn, Naa, Nrhs   !: used by nemo_init

   !! Free surface parameters
   !! =======================
   LOGICAL , PUBLIC :: ln_dynspg_exp    !: Explicit free surface flag
   LOGICAL , PUBLIC :: ln_dynspg_ts     !: Split-Explicit free surface flag

   !! Time splitting parameters
   !! =========================
   LOGICAL,  PUBLIC :: ln_bt_auto       !: Set number of barotropic iterations automatically
   INTEGER,  PUBLIC :: nn_bt_flt        !: Choice of dissipation method for barotropic sub-stepping
   INTEGER,  PUBLIC :: nn_bt_frc        !: order of the Adams-Bashforth extrapolation of the advective forcing of the barotropic mode
   INTEGER,  PUBLIC :: nn_e          !: Number of barotropic iterations during one baroclinic step (rn_Dt)
   REAL(wp), PUBLIC :: rn_bt_cmax       !: Maximum allowed courant number (used if ln_bt_auto=T)
   REAL(wp), PUBLIC :: rn_bt_alpha      !: Time stepping diffusion parameter


   !                                   !!! associated variables
   REAL(wp), PUBLIC ::   rDt, r1_Dt     !: Current model timestep and reciprocal

   !!----------------------------------------------------------------------
   !! space domain parameters
   !!----------------------------------------------------------------------
   LOGICAL         , PUBLIC ::   l_Iperio, l_Jperio   ! i- j-periodicity
   LOGICAL         , PUBLIC ::   l_NFold              ! North Pole folding
   CHARACTER(len=1), PUBLIC ::   c_NFtype             ! type of North pole Folding: T or F point

   ! Tiling namelist
   LOGICAL, PUBLIC ::   ln_tile
   INTEGER         ::   nn_ltile_i, nn_ltile_j

   ! Domain tiling
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) ::   ntsi_a       !: start of internal part of tile domain
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) ::   ntsj_a       !
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) ::   ntei_a       !: end of internal part of tile domain
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) ::   ntej_a       !
   LOGICAL, PUBLIC  ::   l_istiled = .FALSE.   ! whether tiling is currently active or not, default definition needed by timing.F90

   !                             !: domain MPP decomposition parameters
   INTEGER              , PUBLIC ::   nimpp, njmpp     !: i- & j-indexes for mpp-subdomain left bottom
   INTEGER              , PUBLIC ::   narea            !: number for local area (starting at 1) = MPI rank + 1
   INTEGER              , PUBLIC ::   nimpi, njmpi     !: i and j position in the MPI domain decomposition
   INTEGER,               PUBLIC ::   nidom      !: IOIPSL things...

   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   mig        !: local ==> global domain, i-index
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   mjg        !: local ==> global domain, j-index
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   mi0, mi1   !: global ==> local domain, i-index
   !                                                                !:    (mi0=1 and mi1=0 if global index not in local domain)
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   mj0, mj1   !: global ==> local domain, j-index
   !                                                                !:    (mj0=1 and mj1=0 if global index not in local domain)
   INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) ::   nfimpp, nfproc, nfjpi, nfni_0

   !!----------------------------------------------------------------------
   !! horizontal curvilinear coordinate and scale factors
   !! ---------------------------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE,         DIMENSION(:,:) ::   glamt , glamu, glamv , glamf    !: longitude at t, u, v, f-points [degree]
   REAL(wp), PUBLIC, ALLOCATABLE,         DIMENSION(:,:) ::   gphit , gphiu, gphiv , gphif    !: latitude  at t, u, v, f-points [degree]
   REAL(wp), PUBLIC, ALLOCATABLE, TARGET, DIMENSION(:,:) ::   e1t   , e2t  , r1_e1t, r1_e2t   !: t-point horizontal scale factors    [m]
   REAL(wp), PUBLIC, ALLOCATABLE, TARGET, DIMENSION(:,:) ::   e1u   , e2u  , r1_e1u, r1_e2u   !: horizontal scale factors at u-point [m]
   REAL(wp), PUBLIC, ALLOCATABLE, TARGET, DIMENSION(:,:) ::   e1v   , e2v  , r1_e1v, r1_e2v   !: horizontal scale factors at v-point [m]
   REAL(wp), PUBLIC, ALLOCATABLE, TARGET, DIMENSION(:,:) ::   e1f   , e2f  , r1_e1f, r1_e2f   !: horizontal scale factors at f-point [m]
   !
   REAL(wp), PUBLIC, ALLOCATABLE,         DIMENSION(:,:) ::   e1e2t , r1_e1e2t                !: associated metrics at t-point
   REAL(wp), PUBLIC, ALLOCATABLE,         DIMENSION(:,:) ::   e1e2u , r1_e1e2u , e2_e1u       !: associated metrics at u-point
   REAL(wp), PUBLIC, ALLOCATABLE,         DIMENSION(:,:) ::   e1e2v , r1_e1e2v , e1_e2v       !: associated metrics at v-point
   REAL(wp), PUBLIC, ALLOCATABLE,         DIMENSION(:,:) ::   e1e2f , r1_e1e2f                !: associated metrics at f-point
   !
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   ff_f  , ff_t                    !: Coriolis factor at f- & t-points  [1/s]
   
   !!----------------------------------------------------------------------
   !! vertical coordinate and scale factors
   !! ---------------------------------------------------------------------
#if defined key_qco
   LOGICAL, PUBLIC, PARAMETER ::   lk_qco    = .TRUE.   !: qco key flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_qco    = .FALSE.  !: qco key flag
#endif
#if defined key_linssh
   LOGICAL, PUBLIC ::   lk_linssh = .TRUE.   !: linssh key flag
#else
   LOGICAL, PUBLIC ::   lk_linssh = .FALSE.  !: linssh key flag
#endif
#if defined key_ALE
   LOGICAL, PUBLIC, PARAMETER ::   lk_ALE    = .TRUE.   !: ALE key flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_ALE    = .FALSE.  !: ALE key flag
#endif
#if defined key_vco_1d
   LOGICAL, PUBLIC, PARAMETER ::   lk_vco_1d   = .TRUE.   !: 1d key flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_vco_1d   = .FALSE.  !: 1d key flag
#endif
#if defined key_vco_1d3d
   LOGICAL, PUBLIC, PARAMETER ::   lk_vco_1d3d = .TRUE.   !: 1d3d key flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_vco_1d3d = .FALSE.  !: 1d3d key flag
#endif
#if defined key_vco_3d
   LOGICAL, PUBLIC, PARAMETER ::   lk_vco_3d   = .TRUE.   !: 3d key flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_vco_3d   = .FALSE.  !: 3d key flag
#endif

   LOGICAL, PUBLIC ::   l_zco       !: z-coordinate - full step
   LOGICAL, PUBLIC ::   l_zps       !: z-coordinate - partial step
   LOGICAL, PUBLIC ::   l_sco       !: s-coordinate or hybrid z-s coordinate

   !        !-----------------------------------!
   !        !  split of time & space variation  !      coord(i,j,k,t) = coord_0(i,j,k) * (1+rt(i,j,t))
   !        !-----------------------------------!

   !                    !==  reference thickness of ocean water column and its inverse  ==!   (all coord. except ALE/z-tilde)
   !
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   ht_0, r1_ht_0   !: t-depth   [m] and [1/m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   hu_0, r1_hu_0   !: u-depth   [m] and [1/m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   hv_0, r1_hv_0   !: v-depth   [m] and [1/m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   hf_0, r1_hf_0   !: f-depth   [m] and [1/m]
   

   !                    !==  1D reference coordinate  ==!   used in zco and zps cases (z- or z-partial cell coord.)
   !
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)     ::   gdept_1d, e3t_1d   !: reference depth & scale factor of T-level points [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)     ::   gdepw_1d, e3w_1d   !: reference depth & scale factor of W-level points [m]


   !                    !==  3D reference coordinate  ==!   used in zps and sco cases (z-partial cell ans s- coord.)
   !
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   gdept_3d   !: t- depth                [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   gdepw_3d   !: w- depth                [m]
   !
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::     e3t_3d   !: t- vert. scale factor   [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::     e3u_3d   !: u- vert. scale factor   [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::     e3v_3d   !: v- vert. scale factor   [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::     e3f_3d   !: f- vert. scale factor   [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::     e3w_3d   !: w- vert. scale factor   [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::    e3uw_3d   !: uw-vert. scale factor   [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::    e3vw_3d   !: vw-vert. scale factor   [m]

   !                    !==  time varying factor  ==!   used in qco case (quasi-eulerian coordinate) 
   !                                                                           !!!  time-dependent ratio ssh / h_0   (domqco)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)   ::   r3t, r3u, r3v         !: time-dependent    ratio at t-, u- and v-point [-]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)     ::   r3f                   !: mid-time-level    ratio at f-point            [-]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)     ::   r3t_f, r3u_f, r3v_f   !: now time-filtered ratio at t-, u- and v-point [-]
   
   !        !----------------------------------------!
   !        !  combined of time & space variations   !      coord(i,j,k,t)
   !        !----------------------------------------!

   !                    !==  z-tilde or ALE case  ==!   default : time-dependent coord.   (currently use of domvvl old staff)
   !
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   ht          !: t-points           [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   hu, r1_hu   !: u-depth            [m] and [1/m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   hv, r1_hv   !: v-depth            [m] and [1/m]
   !
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:) ::   gdept     !: t- depth                [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:) ::   gdepw     !: w- depth                [m]
   !
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:) ::       e3t   !: t- vert. scale factor   [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:) ::       e3u   !: u- vert. scale factor   [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:) ::       e3v   !: v- vert. scale factor   [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)   ::       e3f   !: f- vert. scale factor   [m]   (only need at Nmm)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:) ::       e3w   !: w- vert. scale factor   [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:) ::      e3uw   !: uw-vert. scale factor   [m]
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:,:) ::      e3vw   !: vw-vert. scale factor   [m]

   
!!gm this is to be replaced by a 2D field for sco  
   INTEGER, PUBLIC ::   nla10              !: deepest    W level Above  ~10m (nlb10 - 1)
   INTEGER, PUBLIC ::   nlb10              !: shallowest W level Bellow ~10m (nla10 + 1)

!!gm bathy should be removed  use ht_0 or ht-ssh
#if defined key_isf
   LOGICAL, PUBLIC, PARAMETER ::   lk_isf    = .TRUE.   !: isf key flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_isf    = .FALSE.  !: isf key flag
#endif
   LOGICAL, PUBLIC            ::   ln_isfcav = .FALSE.  !: absence of ISF (init for debug mode)


   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   risfdep, bathy

   !!----------------------------------------------------------------------
   !! masks, top and bottom ocean point position
   !! ---------------------------------------------------------------------
!!gm Proposition of new name for top/bottom vertical indices
!   INTEGER , PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   mtk_t, mtk_u, mtk_v   !: top    first wet T-, U-, and V-level (ISF)
!   INTEGER , PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   mbk_t, mbk_u, mbk_v   !: bottom last  wet T-, U-, and V-level
!!gm
   INTEGER , PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   mbkt, mbku, mbkv, mbkf   !: bottom last wet T-, U-, V- and F-level
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   tmask_i                  !: interior (excluding halos+duplicated points) domain T-point mask

   INTEGER , PUBLIC, ALLOCATABLE, DIMENSION(:,:) ::   mikt, miku, mikv, mikf   !: top first wet T-, U-, V-, F-level           (ISF)

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)           ::   smask0                              !: surface mask at T-pts on inner domain
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)           ::   smask0_i                            !: equivalent of tmask_i for inner domain
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)           ::   ssmask, ssumask, ssvmask, ssfmask   !: surface mask at T-,U-, V- and F-pts
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:), TARGET ::   tmask, umask, vmask, wmask, fmask   !: land/ocean mask at T-, U-, V-, W- and F-pts
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:), TARGET ::   wumask, wvmask                      !: land/ocean mask at WU- and WV-pts
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:), TARGET ::   fe3mask                             !: land/ocean mask at F-pts
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: tmask_upd, umask_upd, vmask_upd                 !: land/ocean mask at F-pts
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: tmask_agrif                !: agrif mask at T-points excluding ghosts and updated areas 

   !!----------------------------------------------------------------------
   !! calendar variables
   !! ---------------------------------------------------------------------
   INTEGER , PUBLIC ::   nyear         !: current year
   INTEGER , PUBLIC ::   nmonth        !: current month
   INTEGER , PUBLIC ::   nday          !: current day of the month
   INTEGER , PUBLIC ::   nhour0        !: time of day at the start of the run: hour
   INTEGER , PUBLIC ::   nminute0      !:    and minute
   INTEGER , PUBLIC ::   ndastp        !: time step date in yyyymmdd format
   INTEGER , PUBLIC ::   nday_year     !: current day counted from jan 1st of the current year
   INTEGER , PUBLIC ::   nsec_year     !: seconds between 00h jan 1st of the current  year and half of the current time step
   INTEGER , PUBLIC ::   nsec_month    !: seconds between 00h 1st day of the current month and half of the current time step
   INTEGER , PUBLIC ::   nsec_monday   !: seconds between 00h         of the last Monday   and half of the current time step
   INTEGER , PUBLIC ::   nsec_day      !: seconds between 00h         of the current   day and half of the current time step
   REAL(wp), PUBLIC ::   fjulday       !: current julian day
   REAL(wp), PUBLIC ::   fjulstartyear !: first day of the current year in julian days
   REAL(wp), PUBLIC ::   adatrj        !: number of elapsed days since the begining of the whole simulation
   !                                   !: (cumulative duration of previous runs that may have used different time-step size)
   INTEGER , PUBLIC, DIMENSION(  0: 2) ::   nyear_len     !: length in days of the previous/current/next year
   INTEGER , PUBLIC, DIMENSION(-11:25) ::   nmonth_len    !: length in days of the months of the current year
   INTEGER , PUBLIC, DIMENSION(-11:25) ::   nmonth_beg    !: second since Jan 1st 0h of the current year and the half of the months
   INTEGER , PUBLIC                  ::   nsec1jan000     !: second since Jan 1st 0h of nit000 year and Jan 1st 0h the current year
   INTEGER , PUBLIC                  ::   nsec000_1jan000   !: second since Jan 1st 0h of nit000 year and nit000
   INTEGER , PUBLIC                  ::   nsecend_1jan000   !: second since Jan 1st 0h of nit000 year and nitend

   !!----------------------------------------------------------------------
   !! variable defined here to avoid circular dependencies...
   !! ---------------------------------------------------------------------
   INTEGER, PUBLIC ::   nbasin         ! number of basin to be considered in diaprt (glo, atl, pac, ind, ipc)

   !!----------------------------------------------------------------------
   !! agrif domain
   !!----------------------------------------------------------------------
#if defined key_agrif
   LOGICAL, PUBLIC, PARAMETER ::   lk_agrif = .TRUE.    !: agrif flag
   LOGICAL, PUBLIC            ::   lk_south, lk_north, lk_west, lk_east !: Child grid boundaries (interpolation or not)
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_agrif = .FALSE.   !: agrif flag
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#if ! defined key_agrif
   !!----------------------------------------------------------------------
   !! NOT 'key_agrif'      dummy function                     No AGRIF zoom
   !!----------------------------------------------------------------------
   LOGICAL FUNCTION Agrif_Root()
      Agrif_Root = .TRUE.
   END FUNCTION Agrif_Root

   INTEGER FUNCTION Agrif_Fixed()
      Agrif_Fixed = 0
   END FUNCTION Agrif_Fixed

   CHARACTER(len=3) FUNCTION Agrif_CFixed()
      Agrif_CFixed = '0'
   END FUNCTION Agrif_CFixed
#endif

   INTEGER FUNCTION dom_oce_alloc()
      !!----------------------------------------------------------------------
      INTEGER                ::   ii
      INTEGER, DIMENSION(30) :: ierr
      !!----------------------------------------------------------------------
      ii = 0   ;   ierr(:) = 0
      !
      ii = ii+1
      ALLOCATE( glamt(jpi,jpj) ,    glamu(jpi,jpj) ,  glamv(jpi,jpj) ,  glamf(jpi,jpj) ,     &
         &      gphit(jpi,jpj) ,    gphiu(jpi,jpj) ,  gphiv(jpi,jpj) ,  gphif(jpi,jpj) ,     &
         &       e1t (jpi,jpj) ,     e2t (jpi,jpj) , r1_e1t(jpi,jpj) , r1_e2t(jpi,jpj) ,     &
         &       e1u (jpi,jpj) ,     e2u (jpi,jpj) , r1_e1u(jpi,jpj) , r1_e2u(jpi,jpj) ,     &
         &       e1v (jpi,jpj) ,     e2v (jpi,jpj) , r1_e1v(jpi,jpj) , r1_e2v(jpi,jpj) ,     &
         &       e1f (jpi,jpj) ,     e2f (jpi,jpj) , r1_e1f(jpi,jpj) , r1_e2f(jpi,jpj) ,     &
         &      e1e2t(jpi,jpj) , r1_e1e2t(jpi,jpj)                                     ,     &
         &      e1e2u(jpi,jpj) , r1_e1e2u(jpi,jpj) , e2_e1u(jpi,jpj)                   ,     &
         &      e1e2v(jpi,jpj) , r1_e1e2v(jpi,jpj) , e1_e2v(jpi,jpj)                   ,     &
         &      e1e2f(jpi,jpj) , r1_e1e2f(jpi,jpj)                                     ,     &
         &      ff_f (jpi,jpj) ,    ff_t (jpi,jpj)                                     , STAT=ierr(ii) )
      !
      !           !=====================!
      !           !==  vertical mesh  ==!   gdep and e3 arrays allocation
      !           !=====================!
      !                                   !-----------------------------------!
      IF( lk_qco .OR. lk_linssh ) THEN    !-  split time & space variations  -!    qco or linear ssh
         !                                !-----------------------------------!
         ii = ii+1
         ALLOCATE( ht_0(jpi,jpj) ,    hu_0(jpi,jpj)    ,    hv_0(jpi,jpj)     , hf_0(jpi,jpj) ,       &
            &   r1_ht_0(jpi,jpj) , r1_hu_0(jpi,jpj) ,    r1_hv_0(jpi,jpj),   r1_hf_0(jpi,jpj) ,   STAT=ierr(ii)  )
         !
         IF( lk_qco ) THEN                !* qco only :   time variation factor (ssh/h ratio)
            ii = ii+1
            ALLOCATE( r3t(jpi,jpj,jpt) , r3u(jpi,jpj,jpt) , r3v(jpi,jpj,jpt) , r3f(jpi,jpj) ,   STAT=ierr(ii) )
         ENDIF
         !
         ii = ii+1
         IF(     lk_vco_1d ) THEN         !* zco :   allocate 1d vertical arrays for all gdep and e3 fields
            !
            ALLOCATE(        gdept_1d(jpk) ,        gdepw_1d(jpk) ,       &
               &               e3t_1d(jpk) ,          e3w_1d(jpk) ,   STAT=ierr(ii) )
         ELSEIF( lk_vco_1d3d ) THEN
            !                             !* zps :   allocate 1d vertical arrays for gdep and w-level e3 fields and t-level e3 fields
            ALLOCATE(        gdept_1d(jpk) ,        gdepw_1d(jpk) ,       &
               &               e3t_1d(jpk) ,          e3w_1d(jpk) ,       &
               &       e3t_3d(jpi,jpj,jpk) ,  e3u_3d(jpi,jpj,jpk) ,       &
               &       e3v_3d(jpi,jpj,jpk) ,  e3f_3d(jpi,jpj,jpk) ,   STAT=ierr(ii) )
         ELSEIF( lk_vco_3d ) THEN
            !                             !* sco :   allocate 3d vertical arrays for all gdep and e3 fields (no more _1d)
            ALLOCATE(        gdept_1d(jpk) ,        gdepw_1d(jpk) ,       &
               &               e3t_1d(jpk) ,          e3w_1d(jpk) ,       &
               &     gdept_3d(jpi,jpj,jpk) ,gdepw_3d(jpi,jpj,jpk) ,       &
               &       e3t_3d(jpi,jpj,jpk) ,  e3u_3d(jpi,jpj,jpk) ,       &
               &       e3v_3d(jpi,jpj,jpk) ,  e3f_3d(jpi,jpj,jpk) ,       &
               &       e3w_3d(jpi,jpj,jpk) , e3uw_3d(jpi,jpj,jpk) ,       &
               &                             e3vw_3d(jpi,jpj,jpk) ,   STAT=ierr(ii) )
         ENDIF
         !                                !-------------------------------------!
      ELSEIF( lk_ALE ) THEN               !-  combine time & space variations  -!   (vertical ALE coordinate)
         !                                !-------------------------------------!      NOT yet implemented
         ii = ii+1
         ALLOCATE(    ht(jpi,jpj,jpt) ,    hu(jpi,jpj,jpt) ,      hv(jpi,jpj,jpt) ,       &
            &                           r1_hu(jpi,jpj,jpt) , r1_hv  (jpi,jpj,jpt) ,   STAT=ierr(ii)  )
         !
         ii = ii+1
         ALLOCATE( gdept(jpi,jpj,jpk,jpt) , gdepw(jpi,jpj,jpk,jpt) ,       &
            &        e3t(jpi,jpj,jpk,jpt) ,   e3u(jpi,jpj,jpk,jpt) ,       &
            &        e3v(jpi,jpj,jpk,jpt) ,   e3f(jpi,jpj,jpk)     ,       &
            &        e3w(jpi,jpj,jpk,jpt) ,  e3uw(jpi,jpj,jpk,jpt) ,       &
            &                                e3vw(jpi,jpj,jpk,jpt) ,   STAT=ierr(ii) )
      ENDIF
         !
      ii = ii+1
      ALLOCATE( risfdep(jpi,jpj) , bathy(jpi,jpj) , STAT=ierr(ii)  )
         !
      ii = ii+1
      ALLOCATE( tmask_i(jpi,jpj) , smask0(A2D(0)) , smask0_i(A2D(0)) , &
         &      ssmask (jpi,jpj) , ssumask(jpi,jpj) , ssvmask(jpi,jpj) , ssfmask(jpi,jpj) ,     &
         &      mbkt   (jpi,jpj) , mbku   (jpi,jpj) , mbkv   (jpi,jpj) , mbkf(jpi,jpj)    , STAT=ierr(ii) )
         !
      ii = ii+1
      ALLOCATE( mikt(jpi,jpj), miku(jpi,jpj), mikv(jpi,jpj), mikf(jpi,jpj), STAT=ierr(ii) )
         !
      ii = ii+1
      ALLOCATE( tmask(jpi,jpj,jpk) , umask(jpi,jpj,jpk) ,     &
         &      vmask(jpi,jpj,jpk) , fmask(jpi,jpj,jpk) , fe3mask(jpi,jpj,jpk), STAT=ierr(ii) )
         !
      ii = ii+1
      ALLOCATE( wmask(jpi,jpj,jpk) , wumask(jpi,jpj,jpk), wvmask(jpi,jpj,jpk) , STAT=ierr(ii) )
         !
#if defined key_agrif
      ii = ii+1
      ALLOCATE( tmask_upd(jpi,jpj) , umask_upd(jpi,jpj), vmask_upd(jpi,jpj), & 
         &      tmask_agrif(jpi,jpj), STAT=ierr(ii) )
#endif
      !
      dom_oce_alloc = MAXVAL(ierr)
      !
   END FUNCTION dom_oce_alloc


   SUBROUTINE dom_oce_dealloc()

      DEALLOCATE(  &
         &   glamt , glamu, glamv , glamf ,   &
         &   gphit , gphiu, gphiv , gphif ,   &
         &   e1t   , e2t  , r1_e1t, r1_e2t,   &
         &   e1u   , e2u  , r1_e1u, r1_e2u,   &
         &   e1v   , e2v  , r1_e1v, r1_e2v,   &
         &   e1f   , e2f  , r1_e1f, r1_e2f,   &
         &   e1e2t , r1_e1e2t             ,   &
         &   e1e2u , r1_e1e2u , e2_e1u    ,   &
         &   e1e2v , r1_e1e2v , e1_e2v    ,   &
         &   e1e2f , r1_e1e2f             ,   &
         &   ff_f  , ff_t  )
      
      IF( ALLOCATED(ht_0    ) )   DEALLOCATE( ht_0, r1_ht_0, hu_0, r1_hu_0, hv_0, r1_hv_0, hf_0, r1_hf_0 )
      IF( ALLOCATED(gdept_3d) )   DEALLOCATE( gdept_3d, gdepw_3d)
      IF( ALLOCATED(e3t_3d  ) )   DEALLOCATE( e3t_3d, e3u_3d, e3v_3d, e3f_3d  )
      IF( ALLOCATED(e3w_3d  ) )   DEALLOCATE( e3w_3d, e3uw_3d, e3vw_3d )
      IF( ALLOCATED(r3t     ) )   DEALLOCATE( r3t, r3u, r3v, r3f )
      IF( ALLOCATED(r3t_f   ) )   DEALLOCATE( r3t_f, r3u_f, r3v_f )
      IF( ALLOCATED(ht      ) )   DEALLOCATE( ht, hu, hv, r1_hu, r1_hv )
      IF( ALLOCATED(gdept   ) )   DEALLOCATE( gdept, gdepw, e3t, e3u, e3v, e3f, e3w,e3uw, e3vw )

      DEALLOCATE(                                      &
         &   risfdep , bathy   ,                       &
         &   tmask_i , smask0  , smask0_i,             &
         &   ssmask  , ssumask , ssvmask , ssfmask ,   &
         &   mbkt    , mbku    , mbkv    , mbkf    ,   &
         &   mikt    , miku    , mikv    , mikf    ,   &
         &   tmask   , umask   , vmask   , fmask   ,   &
         &   wmask   , wumask  , wvmask  , fe3mask )

#if defined key_agrif
      DEALLOCATE( tmask_upd , umask_upd, vmask_upd, tmask_agrif )
#endif
 
   END SUBROUTINE dom_oce_dealloc
   
   !!======================================================================
END MODULE dom_oce
