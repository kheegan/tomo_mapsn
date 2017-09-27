function qso_spaceden_z, z, gmag_bright, gmag_faint, $
                         M_lim_faint=M_lim_faint, M_lim_bright= $
                         M_lim_bright
; Integrate over observed g-magnitude to give QSO space density at
; fixed z
;
; Use lookup table to increase speed
common qsolumfunc, M_grid, phi_grid

; First, calculate luminosity distance
d_l = lumdist(z, Omega_m = 0.267d, k=0.d,H0 = 71.d,/silent)

; Compute K-correction using Croom et al 2009 approximation
Kcorr = -2.5d * (-0.5d) * alog10(1.+z) - 0.59640d

M_lim_faint = gmag_faint - 5.d*alog10(d_l * 1.d6/10.d) - Kcorr
M_lim_bright  = gmag_bright  - 5.d*alog10(d_l * 1.d6/10.d) - Kcorr

if not keyword_set(M_grid) then begin
   n_points = round((29. - 21.) / 0.002d)+1
   
   M_grid = -29. + 0.002d * dindgen(n_points)
   
   
   phi_grid = dblarr(n_points)
   
   for ii=0, n_points - 1 do begin
      
      phi_grid[ii] = qso_lumfunc_bossmmt((M_grid[ii])[0], z)
      
   endfor
endif

Mcut = where(M_grid GE M_lim_bright AND M_grid LT M_lim_faint)
M_grid_out = M_grid[Mcut]
phi_grid_out = phi_grid[Mcut]

return, int_tabulated(M_grid_out, phi_grid_out)
end
