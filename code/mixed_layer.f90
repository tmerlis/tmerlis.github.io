module mixed_layer_mod

!
! Implementation of mixed layer boundary condition
!

use                  fms_mod, only: set_domain, write_version_number, &
                                    mpp_pe, mpp_root_pe, error_mesg, FATAL, WARNING

use                  fms_mod, only: stdlog, check_nml_error, close_file,&
                                    open_namelist_file, stdout, file_exist, &
                                    read_data, write_data, open_file, &
                                    nullify_domain

use            constants_mod, only: HLV, PI, RHO_CP, CP_AIR 

use         diag_manager_mod, only: register_diag_field, send_data

use       time_manager_mod,   only: time_type

use           transforms_mod, only: get_deg_lat, get_deg_lon, grid_domain

use            vert_diff_mod, only: surf_diff_type



implicit none
private
!=================================================================================================================================

character(len=128) :: version= &
'$Id: mixed_layer.f90 $'
      
character(len=128) :: tagname= &
'$Name:  $'
character(len=128), parameter :: mod_name='mixed_layer'

!=================================================================================================================================

public :: mixed_layer_init, mixed_layer, mixed_layer_end

!=================================================================================================================================

logical :: evaporation = .true.
real    :: qflux_amp = 0.0
real    :: qflux_width = 16.0  ! width of qflux region in degrees
real    :: depth = 40.0
logical :: load_qflux = .false.
logical :: do_fixed_sst = .false.
real    :: asym_amp_fixed_sst = 0.0
real    :: asym_lon0_fixed_sst = 90.0
real    :: ts_pole_fixed_sst = 273.16
real    :: ts_delta_fixed_sst = 27.0
real    :: qflux_asym_amp = 0.0 ! default value of 40 for Walker experiments
real    :: qflux_asym_width_x = 30.0
real    :: qflux_asym_width_y = 7.0
real    :: qflux_asym_lon_plus = 90.0
real    :: qflux_asym_lon_minus = 270.0

namelist/mixed_layer_nml/ evaporation, qflux_amp, depth, qflux_width, &
     load_qflux, do_fixed_sst, asym_amp_fixed_sst, ts_pole_fixed_sst, &
     ts_delta_fixed_sst, asym_lon0_fixed_sst, &
     qflux_asym_amp, qflux_asym_width_x, qflux_asym_width_y, &
     qflux_asym_lon_plus, qflux_asym_lon_minus

!=================================================================================================================================


logical :: module_is_initialized =.false.
logical :: used

integer :: iter
integer, dimension(4) :: axes
integer ::                                                                    &
     id_t_surf,            &   ! surface temperature
     id_flux_lhe,          &   ! latent heat flux at surface
     id_flux_oceanq,       &   ! oceanic Q flux 
     id_flux_t                 ! sensible heat flux at surface

real, allocatable, dimension(:,:)   ::                                        &
     ocean_qflux,           &   ! Q-flux 
     rad_lat_2d,            &   ! latitude in radians 
     rad_lon_2d                 ! latitude in radians 

real, allocatable, dimension(:)   ::                                          &
     deg_lat,               &
     deg_lon

real, allocatable, dimension(:,:)   ::                                        &
     gamma_t,               &   ! Used to calculate the implicit
     gamma_q,               &   ! correction to the diffusion in
     fn_t,                  &   ! the lowest layer
     fn_q,                  &   ! 
     en_t,                  &   !
     en_q,                  &   !
     alpha_t,               &   !
     alpha_q,               &   !
     alpha_lw,              &   !
     beta_t,                &   !
     beta_q,                &   !
     beta_lw,               &   !
     t_surf_dependence,     &   !
     corrected_flux,        &   !
     eff_heat_capacity,     &   ! Effective heat capacity
     delta_t_surf               ! Increment in surface temperature

real inv_cp_air




!=================================================================================================================================
contains
!=================================================================================================================================

subroutine mixed_layer_init(is, ie, js, je, num_levels, t_surf, axes, Time)

type(time_type), intent(in)       :: Time
real, intent(inout), dimension(:,:) :: t_surf
integer, intent(in), dimension(4) :: axes
integer, intent(in) :: is, ie, js, je, num_levels 

integer :: i,j
real    :: rad_qwidth, asym_qwidth_x, asym_qwidth_y, lonplus, lonminus
integer:: ierr, io, unit

if(module_is_initialized) return

call write_version_number(version, tagname)

unit = open_namelist_file ()
ierr=1
do while (ierr /= 0)
  read  (unit, nml=mixed_layer_nml, iostat=io, end=10)
  ierr = check_nml_error (io, 'mixed_layer_nml')
enddo
10 call close_file (unit)

if ( mpp_pe() == mpp_root_pe() )   write (stdlog(), nml=mixed_layer_nml)


allocate(rad_lat_2d              (is:ie, js:je))
allocate(rad_lon_2d              (is:ie, js:je))
allocate(ocean_qflux             (is:ie, js:je))

allocate(deg_lat                 (js:je))
allocate(deg_lon                 (is:ie))


allocate(gamma_t                 (is:ie, js:je))
allocate(gamma_q                 (is:ie, js:je))
allocate(en_t                    (is:ie, js:je))
allocate(en_q                    (is:ie, js:je))
allocate(fn_t                    (is:ie, js:je))
allocate(fn_q                    (is:ie, js:je))
allocate(alpha_t                 (is:ie, js:je))
allocate(alpha_q                 (is:ie, js:je))
allocate(alpha_lw                (is:ie, js:je))
allocate(beta_t                  (is:ie, js:je))
allocate(beta_q                  (is:ie, js:je))
allocate(beta_lw                 (is:ie, js:je))
allocate(delta_t_surf            (is:ie, js:je))
allocate(eff_heat_capacity       (is:ie, js:je))
allocate(corrected_flux          (is:ie, js:je))
allocate(t_surf_dependence       (is:ie, js:je))



!
!see if restart file exists for the surface temperature
!
if (file_exist('INPUT/mixed_layer.res.nc')) then

   call nullify_domain()
   call read_data(trim('INPUT/mixed_layer.res'), 't_surf',   t_surf, grid_domain)

else if (file_exist('INPUT/swamp.res')) then
         unit = open_file (file='INPUT/swamp.res', &
                           form='native', action='read')
         call read_data (unit, t_surf)
         call close_file (unit)
  call error_mesg('mixed_layer','mixed_layer restart file not found, using swamp restart file', WARNING)
else
  call error_mesg('mixed_layer','mixed_layer restart file not found', WARNING)
endif

id_t_surf = register_diag_field(mod_name, 't_surf',        &
                                axes(1:2), Time, 'surface temperature','K')
id_flux_t = register_diag_field(mod_name, 'flux_t',        &
                                axes(1:2), Time, 'sensible heat flux up at surface','watts/m2')
id_flux_lhe = register_diag_field(mod_name, 'flux_lhe',        &
                                 axes(1:2), Time, 'latent heat flux up at surface','watts/m2')
id_flux_oceanq = register_diag_field(mod_name, 'flux_oceanq',        &
                                 axes(1:2), Time, 'oceanic Q-flux','watts/m2')

! latitude will be needed for oceanic q flux
call get_deg_lat(deg_lat)
do j=js,je
  rad_lat_2d(:,j) = deg_lat(j)*PI/180.
enddo

! longitude needed for zonally asymmetries
call get_deg_lon(deg_lon)
do i=is,ie
  rad_lon_2d(i,:) = deg_lon(i)*PI/180.
enddo

! calculate ocean Q flux
! Tim M. note 4/25/2012
! new 1/cos(lat) factor so that it integrates to 0
rad_qwidth = qflux_width*PI/180.
ocean_qflux = qflux_amp*(1-2.*rad_lat_2d**2/rad_qwidth**2) * &
     exp(- ((rad_lat_2d)**2/(rad_qwidth)**2)) / cos(rad_lat_2d)

! ocean Q flux with Walker-esque zonal asymmetry
! added as a perturbation to zonally symmetric qflux calculated above
asym_qwidth_x = qflux_asym_width_x*PI/180.
asym_qwidth_y = qflux_asym_width_y*PI/180.
lonplus       = qflux_asym_lon_plus*PI/180.
lonminus      = qflux_asym_lon_minus*PI/180.

! if sign of lonminus, lonplus is confusing, remember:
! positive OHT divergence -> cooling
ocean_qflux = ocean_qflux + &
     qflux_asym_amp*exp( - ( (rad_lon_2d - lonminus)**2/(asym_qwidth_x**2) + &
     (rad_lat_2d)**2/(asym_qwidth_y**2) ) ) - &
     qflux_asym_amp*exp( - ( (rad_lon_2d - lonplus)**2/(asym_qwidth_x**2) + &
     (rad_lat_2d)**2/(asym_qwidth_y**2) ) )

! load Q flux 
if (load_qflux) then
  call read_data('INPUT/ocean_qflux.nc', 'ocean_qflux',  ocean_qflux)
endif

! fixed SST
if (do_fixed_sst) then
 t_surf = max(min(1.5*rad_lat_2d,pi*0.5),-pi*0.5) ! use t_surf as work array
 t_surf = ts_delta_fixed_sst*(1.-0.5*(sin(t_surf)**2+sin(t_surf)**4)) + ts_pole_fixed_sst

 t_surf = t_surf +  asym_amp_fixed_sst*cos(rad_lon_2d-asym_lon0_fixed_sst*PI/180.) &
      * cos(0.5*pi*min(max(rad_lat_2d/(pi/6),-1.),1.))**2
endif

inv_cp_air = 1.0 / CP_AIR 

module_is_initialized = .true.

return
end subroutine mixed_layer_init

!=================================================================================================================================

subroutine mixed_layer (                                               &
     Time,                                                             &
     t_surf,                                                           &
     flux_t,                                                           &
     flux_q,                                                           &
     flux_r,                                                           &
     dt,                                                               &
     net_surf_sw_down,                                                 &
     surf_lw_down,                                                     &
     Tri_surf,                                                         &
     dhdt_surf,                                                        &
     dedt_surf,                                                        &
     dedq_surf,                                                        &
     drdt_surf,                                                        &
     dhdt_atm,                                                         &
     dedq_atm)         




! ---- arguments -----------------------------------------------------------
type(time_type), intent(in)       :: Time
real, intent(in),  dimension(:,:) :: &
     net_surf_sw_down, surf_lw_down
real, intent(in), dimension(:,:) :: &
     flux_t,    flux_q,     flux_r
real, intent(inout), dimension(:,:) :: t_surf
real, intent(in), dimension(:,:) :: &
   dhdt_surf, dedt_surf, dedq_surf, &
   drdt_surf, dhdt_atm, dedq_atm  
real, intent(in) :: dt
type(surf_diff_type), intent(inout) :: Tri_surf



if(.not.module_is_initialized) then
  call error_mesg('mixed_layer','mixed_layer module is not initialized',FATAL)
endif


! Need to calculate the implicit changes to the lowest level delta_q and delta_t
! - see the discussion in vert_diff.tech.ps
                                                                                                                                    
! Care is needed to differentiate between the sensible heat flux and the
! diffusive flux of temperature
                                                                                                                                    
gamma_t = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_t + dhdt_atm * inv_cp_air))
gamma_q = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_q + dedq_atm))
                                                                                                                                 
fn_t = gamma_t * (Tri_surf%delta_t + Tri_surf%dtmass * flux_t * inv_cp_air)
fn_q = gamma_q * (Tri_surf%delta_q + Tri_surf%dtmass * flux_q)
                                                                                                                                 
en_t = gamma_t * Tri_surf%dtmass * dhdt_surf * inv_cp_air
en_q = gamma_q * Tri_surf%dtmass * dedt_surf
                                                                                                                                    
!
! Note flux_sw doesn't depend on surface or lowest layer values
! Note drdt_atm is not used - should be fixed
!
alpha_t = flux_t * inv_cp_air + dhdt_atm * inv_cp_air * fn_t
alpha_q = flux_q + dedq_atm * fn_q
alpha_lw = flux_r
                                                                                                                                 
beta_t = dhdt_surf * inv_cp_air + dhdt_atm * inv_cp_air * en_t
beta_q = dedt_surf + dedq_atm * en_q
beta_lw = drdt_surf

!
! Implement mixed layer surface boundary condition
!
corrected_flux = - net_surf_sw_down - surf_lw_down + alpha_t * CP_AIR + alpha_lw + ocean_qflux
t_surf_dependence = beta_t * CP_AIR + beta_lw


if (evaporation) then
  corrected_flux = corrected_flux + alpha_q * HLV
  t_surf_dependence = t_surf_dependence + beta_q * HLV
endif

!
! Now update the mixed layer surface temperature using an implicit step
!
if (do_fixed_sst == .false.) then

   eff_heat_capacity = depth * RHO_CP + t_surf_dependence * dt

   if (any(eff_heat_capacity .eq. 0.0))  then 
      write(*,*) 'mixed_layer: error', eff_heat_capacity
      call error_mesg('mixed_layer', 'Avoiding division by zero',fatal)
   end if

   delta_t_surf = - corrected_flux  * dt / eff_heat_capacity

   t_surf = t_surf + delta_t_surf
else
   delta_t_surf = 0.
endif ! end of do_fixed_sst false conditional

!
! Finally calculate the increments for the lowest atmospheric layer
!
Tri_surf%delta_t = fn_t + en_t * delta_t_surf
Tri_surf%delta_q = fn_q + en_q * delta_t_surf

!
! Note:
! When using an implicit step there is not a clearly defined flux for a given timestep
!
if(id_t_surf > 0) used = send_data(id_t_surf, t_surf, Time)
if(id_flux_t > 0) used = send_data(id_flux_t, flux_t, Time)
if(id_flux_lhe > 0) used = send_data(id_flux_lhe, HLV * flux_q, Time)
if(id_flux_oceanq > 0)   used = send_data(id_flux_oceanq, ocean_qflux, Time)

end subroutine mixed_layer

!=================================================================================================================================

subroutine mixed_layer_end(t_surf)

real, intent(inout), dimension(:,:) :: t_surf
integer:: unit

if(.not.module_is_initialized) return

! write a restart file for the surface temperature
call nullify_domain()
call write_data(trim('RESTART/mixed_layer.res'), 't_surf',   t_surf, grid_domain)

module_is_initialized = .false.

end subroutine mixed_layer_end

!=================================================================================================================================

end module mixed_layer_mod
