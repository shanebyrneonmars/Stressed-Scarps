function stress_from_temp_results, fname, d = d
  compile_opt idl2

  if (n_elements(d) eq 0) then begin
    print, 'Assuming a grain size of 1mm'
    d = 1d-3
  endif

  restore, fname

  sav_allt = double(sav_allt)
  dT_dt = -0.5d0 * ((shift(sav_allt, 1, 0) - sav_allt) / sp.dt + (sav_allt - shift(sav_allt, -1, 0)) / sp.dt)
  t0 = rebin(transpose(mean(sav_allt, dimension = 1)), n_elements(ls), n_elements(z))
  lss = ((2.0 * ls[0] - ls[1] - ls[-1]) * dindgen(n_elements(ls)) / (n_elements(ls) - 1.d0) + ls + 360.d0) mod 360.d0

  E = 2.339d10 - 6.48d7 * sav_allt
  dE_dT = -6.48d7

  alp = 2.47d-7 * sav_allt - 1.17d-5
  dalp_dT = 2.47d-7

  mu = 0.4d0 ; Poissons Ratio
  r = 8.31d0 ; Universal Gas Constant
  vm = 1.97d-5 ; Molar volume of Ice I
  del = 9.04d-10 ; Width of grain boundary (taken to be 2xBurger's vector)
  dv0 = 9.10d-4 ; Pre-exponential constant for volume diffusion in ice
  db0 = 5.80d-4 ; Pre-exponential constant for grain-boundary diffusion in ice

  p = [-2.d0, -3.d0, 0.0d0, -1.4d0, 0.0d0] ; Grain size dependance exponents for five deformation mechanisms
  n = [1.d0, 1.d0, 2.4d0, 1.8d0, 4.0d0] ; Stress exponents for five deformation mechanisms
  q = [59.4, 49, 60, 49, 60] * 1d3 ; Activation energies for five deformation mechanisms
  a = [(42.d0 * vm / (3.d0 * r)) * dv0 * 1d6, (42.d0 * vm / (3.d0 * r)) * !pi * del * db0 * 1d6, 5.5d7, 3.9d-3, 4.0d5] * (1d-6) ^ n ; Strain rate prefactors including correction for MPa to Pa
  a[2 : 4] = a[2 : 4] / (2.d0 * 3.d0 ^ (-(n[2 : 4] + 1.d0) / 2.d0)) ; Correct 'A' measured in uniaxial compression to true 'A'

  e_diffv = 3.d0 ^ (-(n[0] + 1.d0) / 2.d0) * a[0] * d ^ p[0] * exp(-q[0] / (r * sav_allt)) / sav_allt ; Calculate strain rate coefficients that are temperature dependant
  e_diffb = 3.d0 ^ (-(n[1] + 1.d0) / 2.d0) * a[1] * d ^ p[1] * exp(-q[1] / (r * sav_allt)) / sav_allt
  e_basal = 3.d0 ^ (-(n[2] + 1.d0) / 2.d0) * a[2] * d ^ p[2] * exp(-q[2] / (r * sav_allt))
  e_gbs = 3.d0 ^ (-(n[3] + 1.d0) / 2.d0) * a[3] * d ^ p[3] * exp(-q[3] / (r * sav_allt))
  e_disl = 3.d0 ^ (-(n[4] + 1.d0) / 2.d0) * a[4] * d ^ p[4] * exp(-q[4] / (r * sav_allt))

  e_diffv = E / (1.d0 - mu) * e_diffv
  e_diffb = E / (1.d0 - mu) * e_diffb
  e_basal = E / (1.d0 - mu) * e_basal
  e_gbs = E / (1.d0 - mu) * e_gbs
  e_disl = E / (1.d0 - mu) * e_disl

  aa = 1.d0 / E * dE_dT * dT_dt
  bb = -E / (1.d0 - mu) * (dalp_dT * (sav_allt - t0) + alp) * dT_dt
  nm1 = n - 1.d0 ; Used to raise the deviatoric stress to the power n-1

  sig = rebin([1d0], (size(sav_allt))[1], (size(sav_allt))[2]) ; Define arrays to store stress results
  d_sig_dt = rebin([1d0], (size(sav_allt))[1], (size(sav_allt))[2]) ; Avoid populating them with zeros for numerical niceness

  for jj = 0l, 4l do begin ; Initial Stress state is unknown so run for a few years to close the loop
    for ii = 0l, (size(sav_allt))[1] - 1l do begin
      ss = sig[ii - 1l, *]
      ass = abs(ss)
      strainrate_pf = e_diffv[ii - 1l, *] * ass ^ nm1[0] + e_diffb[ii - 1l, *] * ass ^ nm1[1] + 1.d0 / (1.d0 / (e_basal[ii - 1l, *] * ass ^ nm1[2]) + 1.d0 / (e_gbs[ii - 1l, *] * ass ^ nm1[3])) + e_disl[ii - 1l, *] * ass ^ nm1[4]

      d_sig_dt[ii - 1l, *] = aa[ii - 1l, *] * ss + bb[ii - 1l, *] - strainrate_pf * ss ; 0 = e_thermal + e_elastic + e_viscous (Mellon's answer would have used -bb)
      sig[ii, *] = ss + d_sig_dt[ii - 1l, *] * sp.dt
    endfor
    print, 'Stress, year: ' + strcompress(jj, /re)
  endfor

  aa2 = (2.d0 * mu / E ^ 2) * dE_dT * dT_dt * sig
  bb2 = (dalp_dT * (sav_allt - t0) + alp) * dT_dt
  cc2 = -(2.d0 * mu / E) * d_sig_dt
  ass = abs(sig)
  strainrate_pf = e_diffv * ass ^ nm1[0] + e_diffb * ass ^ nm1[1] + 1.d0 / (1.d0 / (e_basal * ass ^ nm1[2]) + 1.d0 / (e_gbs * ass ^ nm1[3])) + e_disl * ass ^ nm1[4]
  dd2 = -2.d0 * strainrate_pf * sig / (E / (1.d0 - mu))

  e_z = aa2 + bb2 + cc2 + dd2 ; Vertical Strain Rate   (Mellon's answer would have used -bb2)
  d_z = e_z * sp.dt * rebin(transpose(dz), (size(e_z))[1], (size(e_z))[2]) ; Vertical displacement increments
  d_zz = total(d_z, 1, /cumulative) ; Vertical displacement cumulative
  ; dz0  = total(d_zz(*,ci(1):*),2)                                           ; Sum of vertical displacements over all layers - this is the surface displacement

  if (n_elements(ci) eq 1) then dz0 = total(d_zz[*, ci[0] : *], 2) ; Sum of vertical displacements over all layers - this is the surface displacement
  if (n_elements(ci) gt 1) then dz0 = total(d_zz[*, ci[1] : *], 2) ; Sum of vertical displacements over all layers - this is the surface displacement

  ffname = 'stess_' + string(round(d * 1d6), format = '(I04)') + '_' + fname

  save, ci, dz, itt, ls, ltst, mc, sav_allt, sp, z, lss, sig, dz0, filename = ffname

  return, ffname
end