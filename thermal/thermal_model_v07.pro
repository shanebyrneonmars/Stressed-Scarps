function thermal_model_v07,pname


;Implemented
;-----------
; Variable composition with depth
; Increasing layer thickness with depth
; Adjust # timesteps for integer number of days
; Surface temp. from linearized approach
; Automatic choice of dz0 and number of layers and timestep
; Variably implicit numerical scheme
; Option to generate and use flat terrain temperatures and albedos for upwelling scattered-flux calculation
; Did some minor rearrangement to enable running for airless bodies (frost and atmospheric radiation are now switched off if body_name NE MARS)
; Added alternate spice code for insolation calculation for non-Mars bodies
; Added option to use simple approximation of same e, A, T for surrounding flat terrain

; 12/31/2012 Squashed bug where icy soil thermophysical parameters didn't get set
; 12/31/2012 Squashed bug where seasonal frost of zero was not being recorded properly
; 12/31/2012 Set everything I could find that mattered to double precision
; 05/28/2013 Complete rewrite where model setup code is moved to new program
; 05/29/2013 Bug fix - heat flux should be Q*cos(slope) and not just Q - should probably add flux of q*dz*sin(slope) to each layer too...
; 06/03/2013 Changed temperature reset to account for subsurface heatflow
; 06/13/2013 Bug fix - Exported temperatures and rflat fluxes were one timestep off where they should have been
; 06/14/2024 Added cos_i and az to exported results

;
; Restore pname
;

restore,pname

;
; Needed constants
;

d2r    = !dpi/180.d0                  ; Convert degrees to radians
r2d    = 180.d0/!dpi                  ; Convert radians to degrees
sig    = 5.6705119d-8                 ; Stefan Boltzman constant
lfrost = 5.899d5                      ; Delta H of CO2 frost (J/kg)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Numerical Scheme          ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

f  = sp.numeric_scheme
q  = sp.heat_flux * cos(sp.slope*d2r)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Terms for conduction through layers when k and dz vary with depth ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

kdown         = f            * (dt/(rho*cap*dz)) * (2.0*kon*shift(kon,-1)) / (shift(kon,-1)*dz + kon*shift(dz,-1))
kdown(nl-1)   = 0.0
kup           = f            * (dt/(rho*cap*dz)) * (2.0*kon*shift(kon,+1)) / (shift(kon,+1)*dz + kon*shift(dz,+1))
kup(0)        = 0.0
keven         = 1.0 - kup - kdown

kdownim       = -1.0*(1.0-f) * (dt/(rho*cap*dz)) * (2.0*kon*shift(kon,-1)) / (shift(kon,-1)*dz + kon*shift(dz,-1))
kdownim(nl-1) = 0.0
kupim         = -1.0*(1.0-f) * (dt/(rho*cap*dz)) * (2.0*kon*shift(kon,+1)) / (shift(kon,+1)*dz + kon*shift(dz,+1))
kupim(0)      = 0.0
kevenim       = 1.0 - kupim - kdownim


kupzero     = f       * (dt/(rho(0)*cap(0)*dz(0))) * (kon(0)/dz(0)) * 2.0         ; Factor of two because conduction exchanges with surface temp only half a layer away
kupzeroim   = (1.0-f) * (dt/(rho(0)*cap(0)*dz(0))) * (kon(0)/dz(0)) * 2.0         ; NOTE: there's no minus sign here!!!!!!

boundary        = dblarr(nl)
boundary(nl-1L) = (dt/(rho(nl-1L)*cap(nl-1L)*dz(nl-1L))) * Q                      ; Bottom boundary doesn't change with time (constant flux) so set it up here


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Terms for efficient calculation ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

tkoz      = 2.0*kon(0)/dz(0)

finkup    = kupzero    * (sp.emissivity*lwf  +  (1.0-sp.albedo)*(swf+swf2)) / tkoz
finkupim  = kupzeroim  * (sp.emissivity*lwf  +  (1.0-sp.albedo)*(swf+swf2)) / tkoz
finkupim  = shift(finkupim,-1)
finkupt   = finkup + finkupim

finkupima = (sp.emissivity*lwf  +  (1.0-sp.albedo)*swf) / tkoz
finkupim  = shift(finkupim,-1)

tes       = cs2 * 3.0*sp.emissivity*sig / tkoz   ; cs2 is cos(slope/2)^2 or 1.0 depending on whether the simple flux rule is being used or not ; no longer true 12/31/2012
teskup    = kupzero   * tes
teskupim  = kupzeroim * tes
teskupt   = teskup + teskupim

bes       = cs2 * 4.0*sp.emissivity*sig / tkoz
beskup    = kupzero   * bes
beskupim  = kupzeroim * bes

t_iceim   = shift(t_ice,-1)
frost_a   = kupzero*t_ice + kupzeroim*t_iceim 

dM_ice2      = (    f)*dt*(2.0*kon(0)/dz(0))/lfrost
dM_ice2_im   = (1.0-f)*dt*(2.0*kon(0)/dz(0))/lfrost
dM_ice       = dt*(1.0-sp.frost_albedo)*(swf+swf2f)/lfrost + dt*sp.frost_emissivity*(lwf)/lfrost  - cs2f*dt*sp.frost_emissivity*(sig*t_ice^4) /lfrost
dM_ice_t     = f*dM_ice + (1.0-f)*shift(dM_ice,-1)  - dM_ice2*t_ice - dM_ice2_im*t_iceim

dT0_endice   = 4.00 * lfrost/(dz(0)*rho(0)*cap(0)) 
dT0_startice = 0.25 * (dz(0)*rho(0)*cap(0))/lfrost

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Arrays for saving results
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

sav_allt = fltarr(nt,nl)             ; record Temperatures
mc       = dblarr(nt)                ; record co2 frost mass
itt      = dblarr(nt,n_elements(ci)) ; record temperatures at compositional interfaces (including the surface)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Fire up the loops!              ;       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for yrs = 1L,sp.run_cycles do begin
print,"Starting cycle "+strcompress(yrs,/re)
If (sp.flat_calculate NE '') then rflat = transpose([[sp.frost_albedo*swf], [sp.frost_emissivity*sig*t_ice^4]]) ; Pre-populate the flat terrain fluxes with frost covered values, replace these during frost free periods

for i=0L,nt-1L do begin

 mc(i)         = mco2
 itt(i,*)      = vva*t(ci) + vvb*t(ci-1L)  
 itt(i,0)      = ts 
 sav_allt(i,*) = t
 ot            = t

if (mco2 EQ 0.0) then begin                                        ; Non-frosted case
  
 ts3 = ts^3
 ts4 = ts^4
 bts = 1.0/(1.0 + bes*ts3)
 aa  = (finkupima(i) + tes * ts4 ) * bts
 
 b   = beskup   * ts3 * bts
 bim = beskupim * ts3 * bts
 at  = (finkupt(i) + teskupt * ts4 ) * bts


 boundary(0) = at                                                 ; Top boundary condition is conduction to/from the surface (bottom boundary never changes and is set earlier)
 keven(0)    = 1.0 - kdown(0) - b                                 ; Update top right matrix element for explicit calculation
 kevenim(0)  = 1.0 - kdownim(0) + bim                             ; Update top right matrix element for implicit calculation
 t = kup*shift(ot,1) + kdown*shift(ot,-1) + keven*ot + boundary   ; Explicit solution: Handles all conduction and the boundary conditions
 t = trisol(kupim,kevenim,kdownim,t)                              ; Implicit solution, Uses IDL's built in implementation of Gaussian elimination
 ts = aa + bts*t(0)

  if (ts LT t_iceim(i)) then begin                                ; Did frost form?
   mco2 = dT0_startice*(t_iceim(i)-ts)                            ; If yes then make mco2 non-zero
   ts   = t_iceim(i)                                              ; ...and fix ts = t_ice
  endif
 
 If (sp.flat_calculate NE '') then  rflat(*,i) = [sp.albedo*swf(i), sp.emissivity*sig*ts^4] ; If we're calculating the flat terrain behavior then overwrite upwelling fluxes from frost-free terrain

endif else begin                                                  ; Frosted case
 
 bts = 0.0
 aa  = t_iceim(i)
 
 b   = kupzero
 bim = kupzeroim
 at  = frost_a(i)

 boundary(0) = at                                                 ; Top boundary condition is conduction to/from the surface (bottom boundary never changes and is set earlier)
 keven(0)    = 1.0 - kdown(0) - b                                 ; Update top right matrix element for explicit calculation
 kevenim(0)  = 1.0 - kdownim(0) + bim                             ; Update top right matrix element for implicit calculation
 t = kup*shift(ot,1) + kdown*shift(ot,-1) + keven*ot + boundary   ; Explicit solution: Handles all conduction and the boundary conditions
 t = trisol(kupim,kevenim,kdownim,t)                              ; Implicit solution, Uses IDL's built in implementation of Gaussian elimination
 ts = aa + bts*t(0)

 mco2 = mco2 - dM_ice_t(i) - dM_ice2*ot(0) - dM_ice2_im*t(0)      ; Update Frost Mass
 if (mco2 LT 0) then begin                                        ; Did the surface defrost?
  ts   = ts - mco2*dT0_endice                                     ; If yes then use the left-over energy in top layer...
  mco2 = 0.0d0                                                    ; ...and zero the frost counter
 endif 
 
endelse

endfor                                                            ; End timestep loop

if (yrs EQ sp.reset_cycle) then begin                             ; Reset subsurface temperatures
 t(*) = mean(itt(*,0))
 t(0)  = t(0) + 0.5*dz(0)*sp.heat_flux*cos(sp.slope*d2r)/kon(0)
 for i=1,nl-1L do t(i) = t(i) + (0.5*dz(i)/kon(i) + total(dz(0:i-1)/kon(0:i-1)) )*sp.heat_flux* cos(sp.slope*d2r)

 print,'Reset mean to '+string(mean(itt(*,0)))
endif

endfor                                                            ; End yearly loop



If (sp.flat_calculate NE '') then begin
 timflat = tim
 oname = sp.flat_calculate
 rflat(1,*) = shift(rflat(1,*),1)                ; Realized that I'm putting future temps and frost masses into an array that's later used in the present - this shift corrects that
 save,timflat,rflat,filename=oname
 print,'Flat terrain input file saved to: '+oname
 
 if (sp.run_number NE '') then begin 
  oname = sp.run_number+'_results.sav'
  save,ls,ltst,mc,sav_allt,itt,ci,z,dz,sp,cos_i,az,filename=oname
  print,'Model results saved to: '+oname
 endif
endif

If (sp.flat_calculate EQ '') then begin
 oname = sp.run_number+'_results.sav'
 save,ls,ltst,mc,sav_allt,itt,ci,z,dz,sp,cos_i,az,filename=oname
 print,'Model results saved to: '+oname
endif

return,oname

end
