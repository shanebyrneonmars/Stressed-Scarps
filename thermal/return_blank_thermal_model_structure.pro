function return_blank_thermal_model_structure

;
; Create a blank template structure that will contain the parameters
; for thermal model runs
;
; 5/25/2013 Initially created - SB
;

sp = CREATE_STRUCT(   $

'body_name', '',      $
'orbit',      [0.d0,0.d0,0.d0,0.d0,0.d0], $
'latitude', 0.d0,     $
'longitude', 0.d0,    $
'slope', 0.d0,        $
'aspect', 0.d0,       $
'elevation', 0.d0,    $
'flat_flux_file', '', $
'flat_calculate', '', $

'albedo',     0.d0,   $
'emissivity', 0.d0,   $
'heat_flux',  0.d0,   $ 

'layer_top',      [0.d0,0.d0,0.d0,0.d0], $
'layer_inertia',  [0.d0,0.d0,0.d0,0.d0], $
'layer_density',  [0.d0,0.d0,0.d0,0.d0], $
'layer_capacity', [0.d0,0.d0,0.d0,0.d0], $

'daily_layers',     0.d0, $
'annual_layers',    0.d0, $
'layer_growth',     0.d0, $

'dt',                0.d0, $
'start_time',        0.d0, $ 
'run_cycles',        0.d0, $
'reset_cycle',       0.d0, $
'initial_temperature',0.d0,$ 
'numeric_scheme',    0.d0, $

'atm_ir',           0.d0, $
'atm_vis',          0.d0, $
'frost_switch',     0.d0, $
'frost_albedo',     0.d0, $
'frost_emissivity', 0.d0, $

'run_number',       '', $
name='surface_properties')

;save,sp,filename='blank_thermal_input.sav'
return,sp
end
