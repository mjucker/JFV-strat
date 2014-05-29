
module hs_forcing_mod

!-----------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
use mpp_mod, only: input_nml_file
#else
use fms_mod, only: open_namelist_file
#endif

use     constants_mod, only: KAPPA, CP_AIR, GRAV, PI, SECONDS_PER_DAY, RDGAS

use           fms_mod, only: error_mesg, FATAL, file_exist,       &
                             open_namelist_file, set_domain,      &
			     read_data, check_nml_error,          &
                             mpp_pe, mpp_root_pe, close_file,     &
                             write_version_number, stdlog,        &
                             uppercase,&  !pjk
                             mpp_clock_id,mpp_clock_begin,mpp_clock_end,CLOCK_COMPONENT!,mpp_chksum

use  time_manager_mod, only: time_type, get_time

use  diag_manager_mod, only: register_diag_field, send_data

use  field_manager_mod, only: MODEL_ATMOS, parse
use tracer_manager_mod, only: query_method, get_number_tracers
use   interpolator_mod, only: interpolate_type, interpolator_init, &
                              interpolator, interpolator_end, &
                              CONSTANT, INTERP_WEIGHTED_P
implicit none
private

!-----------------------------------------------------------------------
!---------- interfaces ------------

   public :: hs_forcing, hs_forcing_init, hs_forcing_end

   type(interpolate_type),save         ::  heating_source_interp
   type(interpolate_type),save         ::  u_interp
   type(interpolate_type),save         ::  v_interp
   type(interpolate_type),save         ::  temp_interp
   type(interpolate_type),save         ::  tau_interp

!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

   logical :: no_forcing = .false.
   logical :: surface_forcing_input = .false.

   real :: t_zero=315., t_strat=200., delh=60., delv=10., eps=10., sigma_b=0.7
   real :: P00 = 1.e5
   real :: eps_sc !mj eps in seasonal cycle

   real :: ka = -40. !  negative values are damping time in days
   real :: k_strat = -40. ! tjr
   real :: ks =  -4., kf = -1.

   logical :: do_conserve_energy = .true.

   real :: trflux = 1.e-5   !  surface flux for optional tracer
   real :: trsink = -0.     !  damping time for tracer

   character(len=256) :: local_heating_option='' ! Valid options are 'from_file' and 'Isidoro'. Local heating not done otherwise.
   character(len=256) :: local_heating_file=''   ! Name of file relative to $work/INPUT  Used only when local_heating_option='from_file'
   real :: local_heating_srfamp=0.0              ! Degrees per day.   Used only when local_heating_option='Isidoro'
   real :: local_heating_xwidth=10.              ! degrees longitude  Used only when local_heating_option='Isidoro'
   real :: local_heating_ywidth=10.              ! degrees latitude   Used only when local_heating_option='Isidoro'
   real :: local_heating_xcenter=180.            ! degrees longitude  Used only when local_heating_option='Isidoro'
   real :: local_heating_ycenter=45.             ! degrees latitude   Used only when local_heating_option='Isidoro'
   real :: local_heating_vert_decay=1.e4         ! pascals            Used only when local_heating_option='Isidoro'

   logical :: relax_to_specified_wind = .false.
   character(len=256) :: u_wind_file='u', v_wind_file='v' ! Name of files relative to $work/INPUT  Used only when relax_to_specified_wind=.true.

   character(len=256) :: equilibrium_t_option = 'Held_Suarez'
   character(len=256) :: equilibrium_t_file='temp'  ! Name of file relative to $work/INPUT  Used only when equilibrium_t_option='from_file'
   character(len=256) :: equilibrium_tau_option='Held_Suarez' !mj: choose 'Held-Suarez','profile','from_file'
   character(len=256) :: equilibrium_tau_file='tau'  ! Name of file relative to $work/INPUT  Used only when equilibrium_tau_option='from_file'

!
!  standard atmosphere (sat), polar vortex (pv) and TOA sponge layer tjr
!
   logical :: pv_sat_flag   = .true.  ! flag for SAT strat + polar vortex
   logical :: sat_only_flag = .false. ! mj flag for SAT strat without  polar vortex
   real :: pv_phi0  = 60.    !  polar vortex edge location (in degrees)
   real :: pv_dphi  = 10.    !  polar vortex edge width (in degrees)
   real :: pv_gamma = -1.e-3 !  polar vortex lapse rate (in degK/m)
!
   logical :: sponge_flag = .true. !flag for sponge at top of model
   real :: sponge_pbottom = 1.e2    !bottom of sponge layer, where damping is zero (Pa)
   real :: sponge_tau_days  = 1.0   !damping time scale for the sponge (days)

   real :: p_tropopause = 0.1       !tropopause pressure divided by reference pressure

   real :: scaife_damp  = 2.5      !damping time scale for the scaife damp (days)
   logical :: scaife_flag = .false. !flag for scaife damping

   real :: sigma_strat1=-1,sigma_strat2=-2 ! levels defining transition from tka to tk_strat tjr

   logical :: sc_flag = .false. ! flag for seasonal cycle in stratospheric te
   real ::  sc_phi0n,sc_phi0s,sc_dphin,sc_dphis !generalizations, for seasonal cycle option, of pv_phi0 and pv_dphi above
   
!-----------------------------------------------------------------------

   namelist /hs_forcing_nml/  no_forcing, surface_forcing_input,             &
   	    		      t_zero, t_strat, delh, delv, eps,              &
                              sigma_b, ka, ks, kf, do_conserve_energy,       &
                              trflux, trsink, local_heating_srfamp,          &
                              local_heating_xwidth,  local_heating_ywidth,   &
                              local_heating_xcenter, local_heating_ycenter,  &
                              local_heating_vert_decay, local_heating_option,&
                              local_heating_file, relax_to_specified_wind,   &
                              u_wind_file, v_wind_file, equilibrium_t_option,&
                              equilibrium_t_file,                            &
                              pv_sat_flag, pv_phi0, pv_dphi, pv_gamma,       &  ! tjr
                              sponge_flag,sponge_pbottom,sponge_tau_days,    &  ! tjr
                              p_tropopause,scaife_damp, scaife_flag,         &  ! hmchen
                              sigma_strat1,sigma_strat2,k_strat,             &  ! tjr
                              sc_flag,sc_phi0n,sc_phi0s,sc_dphin,sc_dphis,   &
                              equilibrium_tau_option,equilibrium_tau_file     !mj

!-----------------------------------------------------------------------

   character(len=128) :: version='$Id: hs_forcing.f90, 2012/05/24 mj $'
   character(len=128) :: tagname='$Name: riga_201012_mj $'

   real :: tka, tks, vkf
   real :: scdamp
   real :: trdamp, twopi
   real :: tk_strat

   integer :: id_teq, id_tdt, id_udt, id_vdt, id_tdt_diss, id_diss_heat, id_local_heating, id_newtonian_damping
   real    :: missing_value = -1.e10
   real    :: xwidth, ywidth, xcenter, ycenter ! namelist values converted from degrees to radians
   real    :: srfamp ! local_heating_srfamp converted from deg/day to deg/sec
   character(len=14) :: mod_name = 'hs_forcing'

   logical :: module_is_initialized = .false.
   integer :: id_newt_damp1,id_newt_damp2,id_newt_damp3
!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine hs_forcing ( is, ie, js, je, dt, Time, lon, lat, p_half, p_full, &
                         u, v, t, r, um, vm, tm, rm, udt, vdt, tdt, rdt,&
                         mask, kbot )

!-----------------------------------------------------------------------
   integer, intent(in)                        :: is, ie, js, je
      real, intent(in)                        :: dt
 type(time_type), intent(in)                  :: Time
      real, intent(in),    dimension(:,:)     :: lon, lat
      real, intent(in),    dimension(:,:,:)   :: p_half, p_full
      real, intent(in),    dimension(:,:,:)   :: u, v, t, um, vm, tm
      real, intent(in),    dimension(:,:,:,:) :: r, rm
      real, intent(inout), dimension(:,:,:)   :: udt, vdt, tdt
      real, intent(inout), dimension(:,:,:,:) :: rdt

      real, intent(in),    dimension(:,:,:), optional :: mask
   integer, intent(in),    dimension(:,:)  , optional :: kbot
!-----------------------------------------------------------------------
   real, dimension(size(t,1),size(t,2))           :: ps, diss_heat, surface_forcing
   real, dimension(size(t,1),size(t,2),size(t,3)) :: ttnd, utnd, vtnd, teq, pmass
   real, dimension(size(r,1),size(r,2),size(r,3)) :: rst, rtnd
   integer :: i, j, k, kb, n, num_tracers
   logical :: used
   real    :: flux, sink, value
   character(len=128) :: scheme, params

!-----------------------------------------------------------------------
     if (no_forcing) return

     if (.not.module_is_initialized) call error_mesg ('hs_forcing','hs_forcing_init has not been called', FATAL)

!-----------------------------------------------------------------------
!     surface pressure

     if (present(kbot)) then
         do j=1,size(p_half,2)
         do i=1,size(p_half,1)
            kb = kbot(i,j)
            ps(i,j) = p_half(i,j,kb+1)
         enddo
         enddo
     else
            ps(:,:) = p_half(:,:,size(p_half,3))
     endif

!-----------------------------------------------------------------------
!     rayleigh damping of wind components near the surface

      call rayleigh_damping ( Time, ps, p_full, p_half, u, v, utnd, vtnd, mask=mask )

      if (do_conserve_energy) then
         ttnd = -((um+.5*utnd*dt)*utnd + (vm+.5*vtnd*dt)*vtnd)/CP_AIR
         tdt = tdt + ttnd
         if (id_tdt_diss > 0) used = send_data ( id_tdt_diss, ttnd, Time, is, js)
       ! vertical integral of ke dissipation
         if ( id_diss_heat > 0 ) then
          do k = 1, size(t,3)
            pmass(:,:,k) = p_half(:,:,k+1)-p_half(:,:,k)
          enddo
          diss_heat = CP_AIR/GRAV * sum( ttnd*pmass, 3)
          used = send_data ( id_diss_heat, diss_heat, Time, is, js)
         endif
      endif

      udt = udt + utnd
      vdt = vdt + vtnd

      if (id_udt > 0) used = send_data ( id_udt, utnd, Time, is, js)
      if (id_vdt > 0) used = send_data ( id_vdt, vtnd, Time, is, js)

!-----------------------------------------------------------------------
!     thermal forcing for held & suarez (1994) benchmark calculation
      call newtonian_damping ( Time, lat, ps, p_full, p_half, t, ttnd, teq, mask,surface_forcing )

      tdt = tdt + ttnd
      if (id_newtonian_damping > 0) used = send_data(id_newtonian_damping, ttnd, Time, is, js)

      if(trim(local_heating_option) /= '') then
        call local_heating ( Time, is, js, lon, lat, ps, p_full, p_half, ttnd )
        tdt = tdt + ttnd
        if (id_local_heating > 0) used = send_data ( id_local_heating, ttnd, Time, is, js)
      endif

      if (id_tdt > 0) used = send_data ( id_tdt, tdt, Time, is, js)
      if (id_teq > 0) used = send_data ( id_teq, teq, Time, is, js)

!-----------------------------------------------------------------------
!     -------- tracers -------

      call get_number_tracers(MODEL_ATMOS, num_tracers=num_tracers)
      
      if(num_tracers == size(rdt,4)) then
        do n = 1, size(rdt,4)
           flux = trflux
           sink = trsink
           if (query_method('tracer_sms', MODEL_ATMOS, n, scheme, params)) then
              if (uppercase(trim(scheme)) == 'NONE') cycle
              if (uppercase(trim(scheme)) == 'OFF') then
                 flux = 0.; sink = 0.
              else
                 if (parse(params,'flux',value) == 1) flux = value
                 if (parse(params,'sink',value) == 1) sink = value
              endif
           endif
           rst = rm(:,:,:,n) + dt*rdt(:,:,:,n)
           call tracer_source_sink ( flux, sink, p_half, rst, rtnd, kbot )
           rdt(:,:,:,n) = rdt(:,:,:,n) + rtnd
        enddo
      else
        call error_mesg('hs_forcing','size(rdt,4) not equal to num_tracers', FATAL)
      endif

!-----------------------------------------------------------------------

 end subroutine hs_forcing

!#######################################################################

 subroutine hs_forcing_init ( axes, Time, lonb, latb )

!-----------------------------------------------------------------------
!
!           routine for initializing the model with an
!              initial condition at rest (u & v = 0)
!
!-----------------------------------------------------------------------

           integer, intent(in) :: axes(4)
   type(time_type), intent(in) :: Time
   real, intent(in), optional, dimension(:,:) :: lonb, latb
   

!-----------------------------------------------------------------------
   integer  unit, io, ierr

!     ----- read namelist -----

#ifdef INTERNAL_FILE_NML
     read (input_nml_file, nml=hs_forcing_nml, iostat=io)
     ierr = check_nml_error(io, 'hs_forcing_nml')
#else
      if (file_exist('input.nml')) then
         unit = open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=hs_forcing_nml, iostat=io, end=10)
            ierr = check_nml_error (io, 'hs_forcing_nml')
         enddo
  10     call close_file (unit)
      endif
#endif

     if(trim(equilibrium_t_option) == 'from_file')then
        pv_sat_flag=.false.
        sat_only_flag=.false.
     endif

!     ----- write version info and namelist to log file -----

      call write_version_number (version,tagname)
      if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=hs_forcing_nml)

      if (no_forcing) return

      twopi = 2*PI

!     ----- convert local heating variables from degrees to radians -----

      xwidth  = local_heating_xwidth*PI/180.
      ywidth  = local_heating_ywidth*PI/180.
      xcenter = local_heating_xcenter*PI/180.
      ycenter = local_heating_ycenter*PI/180.

!     ----- Make sure xcenter falls in the range zero to 2*PI -----

      xcenter = xcenter - twopi*floor(xcenter/twopi)

!     ----- convert local_heating_srfamp from deg/day to deg/sec ----

      srfamp = local_heating_srfamp/SECONDS_PER_DAY

!     ----- compute coefficients -----

! If positive, damping time units are (1/s),  value is the inverse of damping time.
! If negative, damping time units are (days), value is the damping time. It is converted to (1/s)
      
      if (ka < 0.) then
        tka = -1./(86400*ka)
      else
        tka = ka
      endif
      if (ks < 0.) then
        tks = -1./(86400*ks)
      else
        tks = ks
      endif
      if (kf < 0.) then
        vkf = -1./(86400*kf)
      else
        vkf = kf
      endif

!     ----- for tracers -----

      if (trsink < 0.) trsink = -86400.*trsink
      trdamp = 0.; if (trsink > 0.) trdamp = 1./trsink

!     ----- register diagnostic fields -----

      id_teq = register_diag_field ( mod_name, 'teq', axes(1:3), Time, &
                      'equilibrium temperature (deg K)', 'deg_K'   , &
                      missing_value=missing_value, range=(/50.,400./) )

      id_newtonian_damping = register_diag_field ( mod_name, 'tdt_ndamp', axes(1:3), Time, &
                      'Heating due to newtonian damping (deg/sec)', 'deg/sec' ,    &
                       missing_value=missing_value     )

      id_tdt = register_diag_field ( mod_name, 'tdt', axes(1:3), Time, &
                      'Total heating: newtonian damping + local heating (deg/sec)', 'deg/sec' ,    &
                       missing_value=missing_value     )

      if(trim(local_heating_option) /= '') then
        id_local_heating=register_diag_field ( mod_name, 'local_heating', axes(1:3), Time, &
                        'Local heating (deg/sec)', 'deg/sec' ,    &
                         missing_value=missing_value     )
      endif

      id_udt = register_diag_field ( mod_name, 'udt_rdamp', axes(1:3), Time, &
                      'zonal wind tendency due to rayleigh damping (m/s2)', 'm/s2',       &
                       missing_value=missing_value     )

      id_vdt = register_diag_field ( mod_name, 'vdt_rdamp', axes(1:3), Time, &
                      'meridional wind tendency due to rayleigh damping (m/s2)', 'm/s2',  &
                       missing_value=missing_value     )

      if (do_conserve_energy) then
         id_tdt_diss = register_diag_field ( mod_name, 'tdt_diss_rdamp', axes(1:3), &
                   Time, 'Dissipative heating from Rayleigh damping (deg/sec)', 'deg/sec',&
                   missing_value=missing_value     )

         id_diss_heat = register_diag_field ( mod_name, 'diss_heat_rdamp', axes(1:2), &
                   Time, 'Vertically integrated dissipative heating from Rayleigh damping (W/m2)', 'W/m2')
      endif


     if(trim(local_heating_option) == 'from_file') then
       call interpolator_init(heating_source_interp, trim(local_heating_file)//'.nc', lonb, latb, data_out_of_bounds=(/CONSTANT/))
     endif
     if(trim(equilibrium_t_option) == 'from_file') then
       call interpolator_init (temp_interp, trim(equilibrium_t_file)//'.nc', lonb, latb, data_out_of_bounds=(/CONSTANT/))
     endif
     if(relax_to_specified_wind) then
       call interpolator_init (u_interp,    trim(u_wind_file)//'.nc', lonb, latb, data_out_of_bounds=(/CONSTANT/))
       call interpolator_init (v_interp,    trim(v_wind_file)//'.nc', lonb, latb, data_out_of_bounds=(/CONSTANT/))
     endif
!
! mj read Newtonian time scale from file
!
     if(trim(equilibrium_tau_option) == 'from_file' ) then
        call interpolator_init (tau_interp, trim(equilibrium_tau_file)//'.nc', lonb, latb, data_out_of_bounds=(/CONSTANT/))
     endif
     


     module_is_initialized  = .true.

 end subroutine hs_forcing_init

!#######################################################################

 subroutine hs_forcing_end 

!-----------------------------------------------------------------------
!
!       routine for terminating held-suarez benchmark module
!             (this routine currently does nothing)
!
!-----------------------------------------------------------------------

 if(trim(local_heating_option) == 'from_file') then
   call interpolator_end(heating_source_interp)
 endif

 if(trim(equilibrium_t_option) == 'from_file') then
   call interpolator_end(temp_interp)
 endif

 if(trim(equilibrium_tau_option) == 'from_file') then
   call interpolator_end(tau_interp)
 endif


 if(relax_to_specified_wind) then
   call interpolator_end(u_interp)
   call interpolator_end(v_interp)
 endif
 
 module_is_initialized = .false.

 end subroutine hs_forcing_end

!#######################################################################

 subroutine newtonian_damping ( Time, lat, ps, p_full, p_half, t, tdt, teq, mask,surface_forcing) 
!-----------------------------------------------------------------------
!
!   routine to compute thermal forcing for held & suarez (1994)
!   benchmark calculation.
!
!-----------------------------------------------------------------------

type(time_type), intent(in)         :: Time
real, intent(in),  dimension(:,:)   :: lat, ps, surface_forcing
real, intent(in),  dimension(:,:,:) :: p_full, t, p_half
real, intent(out), dimension(:,:,:) :: tdt, teq
real, intent(in),  dimension(:,:,:), optional :: mask

!-----------------------------------------------------------------------

          real, dimension(size(t,1),size(t,2)) :: &
     sin_lat, sin_lat_2, cos_lat_2, t_star, cos_lat_4, &
     tstr, sigma, the, tfactr, rps, p_norm, dQdistr

       real, dimension(size(t,1),size(t,2),size(t,3)) :: tdamp, dQdamp

       real, dimension(size(t,2),size(t,3)) :: tz,tauz !mj
       real    ::   k0,k1,k2,k3,k4,k5  !mj (pzg corrected)

       integer :: i, j, k, l
       real    :: tcoeff, pref, tcoeff_strat       ! tjr 06/27/03

!----------- US Standard Atmospheric Temperature 1976 ------------------ t.r.
   integer, parameter :: sat_levs = 8
   real, parameter, dimension(sat_levs) :: sat_p = (/ 1.0000000000000, &
                                                      0.2233611050922, &
                                                      0.0540329501078, &
                                                      0.0085666783593, &
                                                      0.0010945601338, &
                                                      0.0006606353133, &
                                                      0.0000390468337, &
                                                      0.0000000000001 /), &
                                           sat_t = (/ 288.150, &
                                                      216.650, &
                                                      216.650, &
                                                      228.650, &
                                                      270.650, &
                                                      270.650, &
                                                      214.650, &
                                                      186.946 /), &
                                           sat_g = (/ -6.5e-3, &
                                                       0.0e-3, &
                                                       1.0e-3, &
                                                       2.8e-3, &
                                                       0.0e-3, &
                                                      -2.8e-3, &
                                                      -2.0e-3, &
                                                       0.0e-3 /)

   real :: t_tropopause =  216.650
   real :: pif = 3.14159265358979/180.
   real, dimension(size(lat,1),size(lat,2)) :: lat_wgt, t_sat, t_pv
   real :: t0n,t0s,t_days,en,es
   integer :: days,seconds
   integer :: hits
! from namelist: delth [degrees], latitudinal extent of temperature vatiation
!                p0 [p_tropopause], < 1, defines upper boundary of temperature variation
!                delT [K], temperature variation
!-----------------------------------------------------------------------
!------------latitudinal constants--------------------------------------

      sin_lat  (:,:) = sin(lat(:,:))
      sin_lat_2(:,:) = sin_lat(:,:)*sin_lat(:,:)
      cos_lat_2(:,:) = 1.0-sin_lat_2(:,:)
      cos_lat_4(:,:) = cos_lat_2(:,:)*cos_lat_2(:,:)
!mj seasonal cycle in eps
      if(sc_flag)then
         t0n=0
         t0s=180
!         if (.not.present(Time)) call error_mesg('newtonian_damping','sc_flag true but time not present',FATAL)
         call get_time(Time,seconds,days)
         t_days = days+seconds/86400
         es = max(0.0,sin((t_days-t0s)*2*acos(-1.)/360));
         en = max(0.0,sin((t_days-t0n)*2*acos(-1.)/360));
         eps_sc = eps*(en - es)
!mj
      else
         eps_sc = eps
      endif

      t_star(:,:) = t_zero - delh*sin_lat_2(:,:) - eps_sc*sin_lat(:,:)
      if ( .not. pv_sat_flag) then
         tstr  (:,:) = t_strat - eps_sc*sin_lat(:,:)
      else
         tstr  (:,:) = t_tropopause
      endif

!-----------------------------------------------------------------------
      if(trim(equilibrium_t_option) == 'from_file') then
         call get_zonal_mean_temp(Time, p_half, tz)
      endif
      if(trim(equilibrium_tau_option) == 'from_file') then !mj
         call get_zonal_mean_tau(Time, p_half, tauz)
      endif
      tcoeff = (tks-tka)/(1.0-sigma_b)
      pref = P00
      rps  = 1./ps

!begin tjr 06/27/03
      tcoeff_strat = (tka-tk_strat)/(sigma_strat1 - sigma_strat2)
!end tjr 06/27/03

      do k = 1, size(t,3)

!  ----- compute equilibrium temperature (teq) -----

      if(equilibrium_t_option == 'from_file') then
         do i=1, size(t,1)
         do j=1, size(t,2)
           teq(i,j,k)=tz(j,k)
         enddo
         enddo
      else if(trim(equilibrium_t_option) == 'Held_Suarez') then
         p_norm(:,:) = p_full(:,:,k)/pref
         the   (:,:) = t_star(:,:) - delv*cos_lat_2(:,:)*log(p_norm(:,:))
         teq(:,:,k) = the(:,:)*(p_norm(:,:))**KAPPA
         teq(:,:,k) = max( teq(:,:,k), tstr(:,:) )
      else if(trim(equilibrium_t_option) == 'Constant') then
      	 the(:,:) = t_zero
	 teq(:,:,k) = the(:,:) - delh*abs(lat(:,:))/maxval(abs(lat(:,:)))
	 teq(:,:,k) = max( teq(:,:,k), tstr(:,:) )
      else
         call error_mesg ('hs_forcing_nml', &
         '"'//trim(equilibrium_t_option)//'"  is not a valid value for equilibrium_t_option',FATAL)
      endif

!  ----- compute damping -----
      sigma(:,:) = p_full(:,:,k)*rps(:,:)
      if(trim(equilibrium_tau_option) == 'from_file') then !mj
         do i=1, size(t,1)
         do j=1, size(t,2)
           tdamp(i,j,k)=1./tauz(j,k)
         enddo
         enddo
      elseif(trim(equilibrium_tau_option) == 'profile') then
         where (sigma(:,:) <= 1.0 .and. sigma(:,:) > sigma_b)
            tfactr(:,:) = tcoeff*(sigma(:,:)-sigma_b)
            tdamp(:,:,k) = tka + cos_lat_4(:,:)*tfactr(:,:)
         elsewhere(sigma_strat1 <= sigma(:,:) .and. sigma(:,:) < sigma_b )
            tdamp(:,:,k) = tka
         elsewhere(sigma_strat2 <= sigma(:,:) .and. sigma(:,:) < sigma_strat1)
            tfactr(:,:) = tcoeff_strat * (sigma(:,:) - sigma_strat2)
            tdamp(:,:,k) = tk_strat + tfactr(:,:)
         elsewhere
            tdamp(:,:,k) = tk_strat*(k0 + k1*p_full(:,:,k) + k2*p_full(:,:,k)**2 + k3*p_full(:,:,k)**3 + k4*p_full(:,:,k)**4 + k5*p_full(:,:,k)**5)
         endwhere
      elseif(trim(equilibrium_tau_option) == 'Held_Suarez')then
         where (sigma(:,:) <= 1.0 .and. sigma(:,:) > sigma_b)
            tfactr(:,:) = tcoeff*(sigma(:,:)-sigma_b)
            tdamp(:,:,k) = tka + cos_lat_4(:,:)*tfactr(:,:)
         elsewhere(sigma_strat1 <= sigma(:,:) .and. sigma(:,:) < sigma_b )
            tdamp(:,:,k) = tka
         elsewhere(sigma_strat2 <= sigma(:,:) .and. sigma(:,:) < sigma_strat1)
            tfactr(:,:) = tcoeff_strat * (sigma(:,:) - sigma_strat2)
            tdamp(:,:,k) = tk_strat + tfactr(:,:)
         elsewhere
            tdamp(:,:,k) = tk_strat
         endwhere
      else
         call error_mesg ('hs_forcing_nml', &
         '"'//trim(equilibrium_tau_option)//'"  is not a valid value for equilibrium_tau_option',FATAL)         
      endif


      enddo
!  ----- add a US_SAT_1976 stratosphere and a polar vortex to teq ----- tjr
      if ( pv_sat_flag ) then

         if (sc_flag) then
            lat_wgt = 0.5*(1-tanh((lat/pif-sc_phi0n)/sc_dphin))*en + 0.5*(1-tanh((lat/pif-sc_phi0s)/sc_dphis))*es
         else
            lat_wgt = 0.5 * ( 1. - tanh((lat/pif-pv_phi0)/pv_dphi) );
         end if

         do k = 1, size(t,3)
            p_norm(:,:) = p_full(:,:,k)/pref

            do l = 1,sat_levs - 1
               where (p_norm(:,:) < p_tropopause .and. p_norm(:,:) > sat_p(l+1) .and. p_norm(:,:) <= sat_p(l)) 
                  t_sat(:,:) = sat_t(l) * (p_norm(:,:)/sat_p(l))**(-rdgas*sat_g(l)/grav);
               end where
            end do
            
         end do
      end if

      do k=1,size(t,3)
         tdt(:,:,k) = -tdamp(:,:,k)*(t(:,:,k)-teq(:,:,k))
      enddo

      if (present(mask)) then
         tdt = tdt * mask
         teq = teq * mask
      endif

!-----------------------------------------------------------------------

 end subroutine newtonian_damping

!#######################################################################

 subroutine rayleigh_damping ( Time, ps, p_full, p_half, u, v, udt, vdt, mask )

!-----------------------------------------------------------------------
!
!           rayleigh damping of wind components near surface
!
!-----------------------------------------------------------------------

type(time_type), intent(in)         :: Time
real, intent(in),  dimension(:,:)   :: ps
real, intent(in),  dimension(:,:,:) :: p_full, p_half, u, v
real, intent(out), dimension(:,:,:) :: udt, vdt
real, intent(in),  dimension(:,:,:), optional :: mask

!-----------------------------------------------------------------------

real, dimension(size(u,1),size(u,2)) :: sigma, vfactr, sfactr, rps

integer :: i,j,k
real    :: vcoeff
real, dimension(size(u,2),size(u,3)) :: uz, vz, um, vm
real :: umean, vmean

real    :: sponge_coeff !used if sponge_flag
real    :: scdamp
!-----------------------------------------------------------------------
!----------------compute damping----------------------------------------

      if(relax_to_specified_wind) then
        call get_zonal_mean_flow(Time, p_half, uz, vz)
      endif

      vcoeff = -vkf/(1.0-sigma_b)
      rps = 1./ps

      do k = 1, size(u,3)
      if (relax_to_specified_wind) then
         do j=1, size(u,2)
            umean=sum(u(:,j,k))/size(u,1)
            vmean=sum(v(:,j,k))/size(v,1)
            udt(:,j,k) = (uz(j,k)-umean)*vkf
            vdt(:,j,k) = (vz(j,k)-vmean)*vkf
         enddo
      else

         sigma(:,:) = p_full(:,:,k)*rps(:,:)

         where (sigma(:,:) <= 1.0 .and. sigma(:,:) > sigma_b)
            vfactr(:,:) = vcoeff*(sigma(:,:)-sigma_b)
            udt(:,:,k)  = vfactr(:,:)*u(:,:,k)
            vdt(:,:,k)  = vfactr(:,:)*v(:,:,k)
         elsewhere
            udt(:,:,k) = 0.0
            vdt(:,:,k) = 0.0
         endwhere

         if (sponge_flag) then   ! t.r.
            sponge_coeff = 1./sponge_tau_days/86400.
            where (p_full(:,:,k) < sponge_pbottom)
               vfactr(:,:) = -sponge_coeff*(sponge_pbottom-p_full(:,:,k))**2/(sponge_pbottom)**2
               udt(:,:,k) = udt(:,:,k) + vfactr(:,:)*u(:,:,k)
               vdt(:,:,k) = vdt(:,:,k) + vfactr(:,:)*v(:,:,k)
            endwhere
         endif

         if (scaife_flag) then
                scdamp = 1./scaife_damp/86400.
                where( (log(0.1) > log(sigma(:,:))) .and. (log(sigma(:,:)) > log(0.01)))
                   sfactr(:,:) = -scdamp * (log(0.1) - log(sigma(:,:)))
                   udt(:,:,k) = udt(:,:,k) + sfactr(:,:)*u(:,:,k)
                   vdt(:,:,k) = vdt(:,:,k) + sfactr(:,:)*v(:,:,k)
                elsewhere(log(sigma(:,:)) <= log(0.01))
                   sfactr(:,:) = -scdamp
                   udt(:,:,k) = udt(:,:,k) + sfactr(:,:)*u(:,:,k)
                   vdt(:,:,k) = vdt(:,:,k) + sfactr(:,:)*v(:,:,k)
                elsewhere
                   sfactr(:,:) = 0.
                   udt(:,:,k) = udt(:,:,k) + sfactr(:,:)*u(:,:,k)
                   vdt(:,:,k) = vdt(:,:,k) + sfactr(:,:)*v(:,:,k)
                endwhere
         endif
      endif
      enddo

      if (present(mask)) then
          udt = udt * mask
          vdt = vdt * mask
      endif

!-----------------------------------------------------------------------

 end subroutine rayleigh_damping

!#######################################################################

 subroutine tracer_source_sink ( flux, damp, p_half, r, rdt, kbot )

!-----------------------------------------------------------------------
      real, intent(in)  :: flux, damp, p_half(:,:,:), r(:,:,:)
      real, intent(out) :: rdt(:,:,:)
   integer, intent(in), optional :: kbot(:,:)
!-----------------------------------------------------------------------
      real, dimension(size(r,1),size(r,2),size(r,3)) :: source, sink
      real, dimension(size(r,1),size(r,2))           :: pmass

      integer :: i, j, kb
      real    :: rdamp
!-----------------------------------------------------------------------

      rdamp = damp
      if (rdamp < 0.) rdamp = -86400.*rdamp   ! convert days to seconds
      if (rdamp > 0.) rdamp = 1./rdamp

!------------ simple surface source and global sink --------------------

      source(:,:,:)=0.0

   if (present(kbot)) then
      do j=1,size(r,2)
      do i=1,size(r,1)
         kb = kbot(i,j)
         pmass (i,j)    = p_half(i,j,kb+1) - p_half(i,j,kb)
         source(i,j,kb) = flux/pmass(i,j)
      enddo
      enddo
   else
         kb = size(r,3)
         pmass (:,:)    = p_half(:,:,kb+1) - p_half(:,:,kb)
         source(:,:,kb) = flux/pmass(:,:)
   endif



     sink(:,:,:) = rdamp*r(:,:,:)
     rdt(:,:,:) = source(:,:,:)-sink(:,:,:)

!-----------------------------------------------------------------------

 end subroutine tracer_source_sink

!#######################################################################

subroutine local_heating ( Time, is, js, lon, lat, ps, p_full, p_half, tdt )

type(time_type), intent(in)         :: Time
integer, intent(in)                 :: is,js
real, intent(in),  dimension(:,:)   :: lon, lat, ps
real, intent(in),  dimension(:,:,:) :: p_full
real, intent(in),  dimension(:,:,:) :: p_half
real, intent(out), dimension(:,:,:) :: tdt

integer :: i, j, k
real :: lon_temp, x_temp, p_factor
real, dimension(size(lon,1),size(lon,2)) :: lon_factor
real, dimension(size(lat,1),size(lat,2)) :: lat_factor
real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)) :: p_half2
do i=1,size(p_half,3)
  p_half2(:,:,i)=p_half(:,:,size(p_half,3)-i+1)
enddo

tdt(:,:,:)=0.

if(trim(local_heating_option) == 'from_file') then
   call interpolator( heating_source_interp, p_half, tdt, trim(local_heating_file))
else if(trim(local_heating_option) == 'Isidoro') then
   do j=1,size(lon,2)
   do i=1,size(lon,1)
     lon_temp = lon(i,j)
     ! Make sure lon_temp falls in the range zero to 2*PI
     x_temp = floor(lon_temp/twopi)
     lon_temp = lon_temp - twopi*x_temp
     lon_factor(i,j) = exp(-.5*((lon_temp-xcenter)/xwidth)**2)
     lat_factor(i,j) = exp(-.5*((lat(i,j)-ycenter)/ywidth)**2)
     do k=1,size(p_full,3)
       p_factor = exp((p_full(i,j,k)-ps(i,j))/local_heating_vert_decay)
       tdt(i,j,k) = srfamp*lon_factor(i,j)*lat_factor(i,j)*p_factor
     enddo
   enddo
   enddo
else
  call error_mesg ('hs_forcing_nml','"'//trim(local_heating_option)//'"  is not a valid value for local_heating_option',FATAL)
endif

end subroutine local_heating

!#######################################################################


!#######################################################################

subroutine get_zonal_mean_flow ( Time, p_half, uz, vz)

type(time_type), intent(in)         :: Time
real, intent(in),  dimension(:,:,:) :: p_half
real, intent(inout), dimension(:,:) :: uz,vz

integer :: j, k
real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)-1) :: uf,vf
call interpolator( u_interp, p_half, uf, trim(u_wind_file))
call interpolator( v_interp, p_half, vf, trim(v_wind_file))

do j=1,size(p_half,2)
do k=1,size(p_half,3)-1
  uz(j,k)=sum(uf(:,j,k))/size(uf,1)
  vz(j,k)=sum(vf(:,j,k))/size(vf,1)
enddo
enddo
end subroutine get_zonal_mean_flow
!#######################################################################

subroutine get_zonal_mean_temp ( Time, p_half, tm)

type(time_type), intent(in)         :: Time
real, intent(in),  dimension(:,:,:) :: p_half
real, intent(inout), dimension(:,:) :: tm

integer :: j, k
real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)-1) :: tf
!call interpolator( temp_interp, p_half, tf, trim(equilibrium_t_file))
call interpolator( temp_interp, Time, p_half, tf, trim(equilibrium_t_file))

do j=1,size(p_half,2)
do k=1,size(p_half,3)-1
  tm(j,k)=sum(tf(:,j,k))/size(tf,1)
enddo
enddo
end subroutine get_zonal_mean_temp
!#######################################################################

subroutine get_zonal_mean_tau ( Time, p_half, taum)

type(time_type), intent(in)         :: Time
real, intent(in),  dimension(:,:,:) :: p_half
real, intent(inout), dimension(:,:) :: taum

integer :: j, k
real, dimension(size(p_half,1),size(p_half,2),size(p_half,3)-1) :: tauf
!call interpolator( tau_interp, p_half, tauf, trim(equilibrium_tau_file))
call interpolator( tau_interp, Time, p_half, tauf, trim(equilibrium_tau_file))

do j=1,size(p_half,2)
do k=1,size(p_half,3)-1
  taum(j,k)=sum(tauf(:,j,k))/size(tauf,1)
enddo
enddo
end subroutine get_zonal_mean_tau
!#######################################################################

end module hs_forcing_mod
