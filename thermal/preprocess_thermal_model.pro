function preprocess_thermal_model,sp

;6/3/2013 changed temperature initialization to account for subsurface heatflow
;10/21/2014 Added Ceres option to preprocess code
;8/1/2015 added generic bodies rather than specific spice calculation
;8/1/2015 sp.body_name must begin with SPICE_ to utilize old method
;12/1/2015 frame_name is now IAU_+body_name, added non-spice parameters for Pluto, Umbriel and Vesta
;1/31/2016 added 'CUSTOM' body for self defined orbits in sp.orbit - note there's problem here in that 'CUSTOM' bodies will use a martian atmosphere (can be disabled with the sp.frost_switch variable)
;1/23/2017 Added equation-of-time calculation to update LMST to LTST (including mars-specific refinement). Fixed minor bug in the eccentric-anomaly calculation (that was having no negative effects).

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Restore run parameters and calculate model inputs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

If (sp.flat_calculate NE '') then begin
 if (sp.flat_flux_file NE '') then begin
  print,'WTF... You cannot use flat terrain behavior at the same time you want to calculate it.'
  retall
 endif

 if (sp.slope NE 0.d0) then begin
  print,'WTF... You cannot calculate flat terrain behavior when the slope is non-zero.'
  retall
 endif
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Define some constants
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

au     = 149597870610.d0              ; Size of one astronomical unit in meters 
flux   = 1367.6d0                     ; Solar Flux at 1 AU
d2r    = !dpi/180.d0                  ; Convert degrees to radians
r2d    = 180.d0/!dpi                  ; Convert radians to degrees
sig    = 5.6705119d-8                 ; Stefan Boltzman constant
lfrost = 5.899d5                      ; Delta H of CO2 frost (J/kg)
gc  = 6.67259d-11                     ; Gravitational constant
sm  = 1.9891d30                       ; Solar mass
yr  = 2.d0*!dpi*sqrt(au^3/(gc*sm))    ; Length of one year in seconds


if (strmid(sp.body_name,0,5) EQ 'SPICE') then begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Restore SPICE kernels needed for insolation calculation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

body_name = (strsplit(sp.body_name,'_',/extract))[1] 

spawn,'echo $HOME',homedir
cspice_furnsh,homedir(0)+'/naif/insolation_kernels.list'
frame_name = 'IAU_'+body_name

;if (body_name EQ 'MARS')  then frame_name = 'IAU_MARS'
;if (body_name EQ 'VESTA') then frame_name = 'VESTA_FIXED'
;if (body_name EQ 'CERES') then frame_name = 'CERES_FIXED'
;if (body_name EQ 'PLUTO') then frame_name = 'IAU_PLUTO'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Setup dummy time array for ~30 years. Use retrieved position vector of Sun in the body
; frame to figure out orbital period by looking at perihelion spacing. Find rotation
; period in the pck.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

tim        = (dindgen(1d4)+0.5d0)*1d5 + sp.start_time
cspice_spkpos, 'SUN', tim, frame_name, 'LT+S', body_name, ptarg, ltime
sol_dist   = reform( sqrt(total(ptarg^2,1)) )

xx     = (shift(sol_dist,1) - sol_dist) * (sol_dist - shift(sol_dist,-1))
xx(-1) = 0.d0
xx(0)  = 0.d0
ptim   = (tim( where((xx LT 0) and (sol_dist LT mean(sol_dist))) ))[0:1]

for i=0,1 do begin
 tim        = (dindgen(3d4)-1.5d4)*10.d0 + ptim(i)
 cspice_spkpos, 'SUN', tim, frame_name, 'LT+S', body_name, ptarg, ltime
 sol_dist   = reform( sqrt(total(ptarg^2,1)) )
 ptim(i) = (tim( where(sol_dist EQ min(sol_dist)) ))[0]
endfor

for i=0,1 do begin
 tim        = (dindgen(3d4)-1.5d4)*1.d-3 + ptim(i)
 cspice_spkpos, 'SUN', tim, frame_name, 'LT+S', body_name, ptarg, ltime
 sol_dist   = reform( sqrt(total(ptarg^2,1)) )
 ptim(i) = (tim( where(sol_dist EQ min(sol_dist)) ))[0]
endfor

myr=ptim(1)-ptim(0)

cspice_BODVRD,body_name,'PM',3,pm
rot = (24.d0*3600.d0*360.d0)/pm(1)                        ; Rotational Period

endif else begin ; This elseif is for the case where the body name doesn't begin with SPICE_

body_name = sp.body_name

CASE 1 OF
(body_name EQ 'MARS'): BEGIN
  a       = 1.52366231d0           ; Orbital semi-major axis
  e       = 0.09341233d0           ; Orbital eccentricity
  lsp     = 250.87d0       * d2r   ; Ls of perihelion
  obl     = 25.19d0        * d2r   ; Obliquity
  rot     = 88775.2d0              ; Length of a solar (not sidereal) day
END

(body_name EQ 'PLUTO'): BEGIN
  a       = 39.537d0           ; Orbital semi-major axis
  e       = 0.2482d0           ; Orbital eccentricity
  lsp     = 4.3d0      * d2r   ; Ls of perihelion
  obl     = 119.6d0     * d2r  ; Obliquity
  rot     = 551854.08d0        ; Length of a solar (not sidereal) day
END

(body_name EQ 'VESTA'): BEGIN
  a       =   2.354         ; Orbital semi-major axis
  e       =   0.091         ; Orbital eccentricity
  lsp     = 202.114 * d2r   ; Ls of perihelion
  obl     =  27.47  * d2r   ; Obliquity
  rot     = 19231.662d0     ; Length of a solar (not sidereal) day
END

(body_name EQ 'CERES'): BEGIN
  a       = 2.76750591440571d0     ; Orbital semi-major axis
  e       = 0.07582276595896797d0  ; Orbital eccentricity
  lsp     = 302.11022d0    * d2r   ; Ls of perihelion
  obl     = 4.03d0         * d2r   ; Obliquity
  rot     = 32667.0d0              ; Length of a solar (not sidereal) day
END

(body_name EQ 'UMBRIEL'): BEGIN
  a       = 19.194d0     ; Orbital semi-major axis
  e       = 0.0429d0  ; Orbital eccentricity
  lsp     = 182.4d0    * d2r   ; Ls of perihelion
  obl     = 97.87d0         * d2r   ; Obliquity
  rot     = 383788.8d0              ; Length of a solar (not sidereal) day
END

(body_name EQ 'CUSTOM'): BEGIN
  a       = sp.orbit(0)        ; Orbital semi-major axis
  e       = sp.orbit(1)        ; Orbital eccentricity
  lsp     = sp.orbit(2)* d2r   ; Ls of perihelion
  obl     = sp.orbit(3)* d2r   ; Obliquity
  rot     = sp.orbit(4)        ; Length of a solar (not sidereal) day
END

ELSE: BEGIN
 print,'Body not recognized, edit the preprocess_thermal_model.pro file.'
 retall
END
ENDCASE

k3l = sqrt(gc*sm/(a*au)^3)
myr = 2.d0*!dpi/k3l                                 ; Length of a year

endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Setup subsurface arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

nclyr = max(where(sp.layer_top GT 0)) + 1L                                                                            ; Number of material layers
if (nclyr EQ 0) then nclyr = 1L                                                                                       ; Special case where material doesn't change with depth
ci = dblarr(nclyr)                                                                                                    ; Layer number where composition 'i' starts

dsd = (sp.layer_inertia(0:nclyr-1L)/(sp.layer_density(0:nclyr-1L) * sp.layer_capacity(0:nclyr-1L))) * sqrt(rot/!dpi)  ; Diurnal skin depth of each material 
asd = (sp.layer_inertia(0:nclyr-1L)/(sp.layer_density(0:nclyr-1L) * sp.layer_capacity(0:nclyr-1L))) * sqrt(myr/!dpi)  ; Annual skin depth of each material

dz0 = dsd(0) / ((1.0-sp.layer_growth^sp.daily_layers)/(1.0-sp.layer_growth))	                            ; Thickness of first layer based on diurnal skin depth of surface material
nl  = ceil(alog(1.0 - (1.0-sp.layer_growth)*(sp.annual_layers*asd(-1)/dz0))/alog(sp.layer_growth) )         ; Number of subsurface layers based on annual skin depth of deepest layer

dz   = dz0*sp.layer_growth^(dindgen(nl))
z    = total(dz,/cumulative) - dz/2.d0 
zbot = total(dz)

rho = dblarr(nl)	
cap = dblarr(nl)	
kon = dblarr(nl)	

for i=0,nclyr-1L do begin
 if (sp.layer_top(i) LT zbot) then begin
  v  = abs((z-dz/2.) - sp.layer_top(i))
  ci(i) = (where(v EQ min(v)))[0]
  
  rho(ci(i):*) = sp.layer_density(i)
  cap(ci(i):*) = sp.layer_capacity(i)
  kon(ci(i):*) = sp.layer_inertia(i)^2.0 / (sp.layer_density(i) * sp.layer_capacity(i))

  print,'Material '+string(i+1,format='(I02)')+' requested to start at '+string(sp.layer_top(i))+'m actually starts at '+string(z(ci(i))-dz(ci(i))/2.)+'m, layer '+strcompress(ci(i),/re) 
 endif
endfor

vva  = (kon(ci)   /dz(ci)   ) / ((kon(ci)/dz(ci)) + (kon(ci-1L)/dz(ci-1L)) )  
vvb  = (kon(ci-1L)/dz(ci-1L)) / ((kon(ci)/dz(ci)) + (kon(ci-1L)/dz(ci-1L)) )  

if (zbot LT max(sp.layer_top)) then print,'WARNING: Model domain does not reach the lowest material specified'

t     = rebin([sp.initial_temperature],nl)       ; Initial guess at sub-surface temperatures
t(0)  = t(0) + 0.5*dz(0)*sp.heat_flux*cos(sp.slope*d2r)/kon(0)
for i=1,nl-1L do t(i) = t(i) + (0.5*dz(i)/kon(i) + total(dz(0:i-1)/kon(0:i-1)) )*sp.heat_flux* cos(sp.slope*d2r)
ot    = t                       ; Old temperatures for explicit calculation
mco2  = 0.0d0                   ; Initial frost mass to jumpstart results
ts    = (3.0*t(0) - t(1))/2.0   ; Initial surface temperature (linear extrapolation)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Courant Test 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

cc = 1.0 / (2.0 * (1.0/(rho*cap*dz)) * (2.0*kon*shift(kon,-1)) / (shift(kon,-1)*dz + kon*shift(dz,-1)) )
cc = min(cc(0:nl-2))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Setup real time array for 1 year. First, adjust number of timesteps so that the model 
; runs for an interger number of days.  Next, setup time array and get Ls, local times 
; and solar distance.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

dt  = sp.dt
if (dt EQ 0.d0) then dt=cc        ; Make the timestep the Courant criteria if unset by the user
nt  = ceil(myr/dt)                ; Guess at the number of required timesteps

if (strmid(sp.body_name,0,5) EQ 'SPICE') then begin   ; USE THE SPICE METHOD

tim = ([0.0,nt-1L]+0.5d0)*dt + sp.start_time 
cspice_spkpos, 'SUN', tim, frame_name, 'LT+S', body_name, ptarg, ltime
ltst       = (12.d0 + (sp.longitude*d2r - atan(ptarg(1,*),ptarg(0,*)))*12.d0/!dpi) mod 24.0

if ((ltst(1)-ltst(0)) GT 12.0) then nt = nt + round(((ltst(0)-ltst(1)+24.0) mod 24.0)/24.0 * rot/dt)
if ((ltst(1)-ltst(0)) LT 12.0) then nt = nt - round(((ltst(1)-ltst(0)+24.0) mod 24.0)/24.0 * rot/dt)

tim = (dindgen(nt)+0.5d0)*dt + sp.start_time  ; Array of times - Add offset for polar studies, 0.4 Myr is ~Ls 125

ls = dblarr(n_elements(tim))
for i=0,n_elements(tim)-1L do ls(i)=cspice_lspcn( body_name, tim(i), 'LT+S' )*r2d 
cspice_spkpos, 'SUN', tim, frame_name, 'LT+S', body_name, ptarg, ltime
ltst       = reform((12.d0 + (sp.longitude*d2r - atan(ptarg(1,*),ptarg(0,*)))*12.d0/!dpi) mod 24.0)
sol_dist   = sqrt(total(ptarg^2,1)) * 1d3 / au
sin_dec    = reform(ptarg(2,*))/sqrt(total(ptarg^2,1))

cspice_kclear
endif else begin                                         ; USE THE GENERIC METHOD

; Make local times agree at the beginning and end of the year

lmst       = ([0.0,nt-1L]+0.5d0)*dt*(24.d0/rot) mod 24.d0
if ((lmst(1)-lmst(0)) GT 12.0) then nt = nt + round(((lmst(0)-lmst(1)+24.0) mod 24.0)/24.0 * rot/dt)
if ((lmst(1)-lmst(0)) LT 12.0) then nt = nt - round(((lmst(1)-lmst(0)+24.0) mod 24.0)/24.0 * rot/dt)

;
; Low-Res calculation of times at uniformly spaced true anomalies
;
n   = 1d4
ta = (dindgen(n)+0.5d)*2.d0*!dpi/double(n)         ; create evenly spaced true anomalies
ea = acos( (e+cos(ta))/(1.d0+e*cos(ta)) )          ; Calculate eccentric anomalies
ma = ea - e*sin(ea)                                ; Calculate mean anomalies
tim_rough  = ma/k3l                                        ; Time along orbital path (will be irregularly spaced)
tim_rough(where(ta GT !dpi)) = myr - tim_rough(where(ta GT !dpi))  ; correction of 2nd half of the orbit

;
; High-Res interpolation of true anomalies at uniformly spaced times
;
tim  = (dindgen(nt) + 0.5d0) * dt   ; Now use these uniformly spaced times
;ta2 = interpol(ta,tim_rough,tim)                                 ; Interpolate the True anomalies at these times
;ea2 = acos( (e+cos(ta2))/(1.d0+e*cos(ta2)) )                     ; Calculate the eccentric anomalies at these times

ta2 = interpol(ta,tim_rough,tim) mod (2.d0*!dpi)                  ; Interpolate the True anomalies at these times
ea2 = acos( (e+cos(ta2))/(1.d0+e*cos(ta2)) )                      ; Calculate the eccentric anomalies at these times
ea2(where(ta2 GT !dpi)) = 2.d0*!dpi - ea2(where(ta2 GT !dpi))     ; Correct eccentric anomalies in second half of the year
ma2 = ea2 - e*sin(ea2)                                            ; Calculate mean anomalies at these times

sol_dist = a*(1.d0 - e*cos(ea2))                                ; Solar Distance
ls       = ((ta2 + lsp) mod (2.d0*!dpi))                       ; Solar longitude
sin_dec  = sin(obl)*sin(ls)                                    ; Solar declination (or rather the sine of this)
cos_dec  = sqrt(1.d0 - sin_dec^2)                                
lmst     = ((tim/rot * 2.d0*!dpi) mod (2.d0*!dpi))             ;local mean solar time

eot = (ma2-ta2)                                                  ; Equation of time for a generic planet
if (body_name EQ 'MARS') then eot = eot + (2.861d0*sin(2.d0*Ls) - 0.071d0*sin(4.d0*Ls) + 0.002d0*sin(6.d0*Ls))*d2r  ; Adding planet-specific perturbations could improve this....

;ltst = lmst * 12.d0/!dpi                            ; I KNOW THIS IS A LITTLE WRONG BECAUSE ORBITAL SPEED VARIES, FIX SOMEDAY...
ltst = ((lmst+eot) * 12.d0/!dpi) mod 24.d0           ; someday=today...  
ls   = ls * r2d

endelse

 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate illumination angles for both flat and sloped terrain
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

lat     = sp.latitude
slope   = sp.slope
aspect  = sp.aspect

fff     = flux / sol_dist^2
cos_hr  = cos((ltst-12.d0)*!dpi/12.d0)
cos_dec = sqrt(1.d0-sin_dec^2)
cos_i   = ( sin(lat*d2r)*sin_dec + cos(lat*d2r)*cos_dec*cos_hr )>(-1.0)<1.0
sin_i   = ( sqrt(1.d0-cos_i^2) )>(-1.0)<1.0

az  = acos(( (sin_dec-sin(lat*d2r)*cos_i)/(cos(lat*d2r)*sin_i) )>(-1)<1)
az(where(ltst GT 12.d0)) = !dpi*2.d0 - az(where(ltst GT 12.d0))
if (max(az) GT  !dpi) then az(where(az GT  !dpi)) = az(where(az GT  !dpi)) - 2.d0*!dpi
if (max(az) LT -!dpi) then az(where(az LT -!dpi)) = az(where(az LT -!dpi)) + 2.d0*!dpi

cos_is  = cos_i*cos(slope*d2r) + sin_i*sin(slope*d2r)*cos(aspect*d2r-az) 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate fluxes - simple atmosphere for now...
; No attenuation of direct beam, scattered vis and emitted IR are constant percentages  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

swff     = (fff*cos_i )>0       ; Flux on a flat surface
swf      = (fff*cos_is)>0       ; Flux on a sloping surface
swf(where(cos_i LE 0)) = 0.0    ; Make sure this is zero during the night

swf = sp.atm_vis * swff * cos(slope*d2r/2.0)^2  + swf                                             ; Add atmospheric scattered visible radiation 
swff= sp.atm_vis * swff                         + swff                                            ; Add atmospheric scattered visible radiation to flat surface
lwf = sp.atm_ir  * cos(slope*d2r/2.0)^2 * (fff*(sin(lat*d2r)*sin_dec + cos(lat*d2r)*cos_dec))>0   ; atm_ir * noontime flux as downwelling IR


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Restore pre-calculated flat terrain fluxes or use simple approx...
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

cs2   = 1.0                                    ; Turn off simple approximation for upwelling fluxes scattered from surrounding terrain
cs2f  = 1.0                                    
swf2  = 0.0
swf2f = 0.0

if ((sp.flat_flux_file NE '') AND (sp.flat_flux_file NE 'simple')) then begin 
 restore,sp.flat_flux_file
 if (max(abs(timflat-tim)) GT 1e-3) then begin
  print,'ERK!!! Flat fluxes file is not compatible - some sort of timing difference...'
  retall
 endif

 swf = swf + rflat(0,*)*sin(slope*d2r/2.d0)^2.0
 lwf = lwf + rflat(1,*)*sin(slope*d2r/2.d0)^2.0
endif 

if (sp.flat_flux_file EQ 'simple') then begin
 cs2  = 1.0 - sp.emissivity      *sin(slope*d2r/2.0)^2.0 ; Use simple approximation for upwelling fluxes emitted from surrounding terrain
 cs2f = 1.0 - sp.frost_emissivity*sin(slope*d2r/2.0)^2.0 ; Use simple approximation for upwelling fluxes emitted from surrounding frosted terrain
 
 swf2  = sp.albedo      *swff*sin(slope*d2r/2.0)^2.0 ; Use simple approximation for upwelling fluxes scattered from surrounding terrain
 swf2f = sp.frost_albedo*swff*sin(slope*d2r/2.0)^2.0 ; Use simple approximation for upwelling fluxes scattered from surrounding frosted terrain
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Atmospheric and CO2 frost-point calculation 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

t_ice = dblarr(n_elements(tim))

if ((sp.body_name EQ 'MARS') OR (sp.body_name EQ 'SPICE_MARS') OR (sp.body_name EQ 'CUSTOM')) then begin
 p_atm_terms = [7.97078,-0.539781,0.468818,0.368771,-0.392702,0.0206071,-0.0224410,-0.0326866,-0.00261966,0.0145776,-0.000184519d0]
 p_atm       = p_atm_terms(0)
 for i=1,5 do p_atm = p_atm + p_atm_terms(i*2-1)*SIN((2.*!pi*i)*(ls/360.0)) + p_atm_terms(i*2)*COS((2.*!pi*i)*(ls/360.0))

 s_atm       = 10800.d0                                   ; Scale-height of the atmosphere
 p_atm       = p_atm*exp(-(sp.elevation+3627.d0)/s_atm)   ; Scale pressure from VL1 site at -3627m
 t_ice       = 3148.0/(23.102 - alog(p_atm))              ; Sublimation Point of CO2 ice
 t_ice       = t_ice * sp.frost_switch                    ; Turn frost off by setting the frost point to zero K if requested
endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Store the things needed to do the thermal model run
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

oname = sp.run_number+'_preprocessed.sav'
save,sp,rho,cap,kon,nl,z,dz,ci,vva,vvb,nt,tim,dt,ls,ltst,swf,swf2,swf2f,lwf,cs2,cs2f,t_ice,t,ot,ts,mco2,sin_dec,sol_dist,filename=oname
print,'Model input file saved as: '+oname

return,oname
end
