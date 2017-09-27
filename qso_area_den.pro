function qso_area_den, z1, z2, gmag=gmag
;; Evaluate area density of quasars (per sq deg) at redshifts between
;; z1 and z2, and brighter than some g-magnitude limit
;;
;; KG Lee (4.12.2014): gmag can now be entered as the faint-end limit
;; instead of as 2-element vector. Assume g=17 as bright limit.

delta_z = 0.005
gmag_0 = 17.

n_z = round((z2 - z1)/delta_z) + 1
zgrid = z1 + delta_z * findgen(n_z)

;; c/H_0/sqrt(Omega_matter) ... multiplying dz with this gives comoving
;; radial distance
dw2dz = 3.e5 / 71. ;/ sqrt(0.267) 

integrand = fltarr(n_z)
Mpc_per_deg = fltarr(n_z)

for ii = 0, n_z - 1 do begin

   arcsec_per_Mpc = zang(1000.,zgrid[ii], H0=71.,Omega_m = 0.267,k=0,/silent)
   
   ;; (1+z) is to convert from physical to comoving
   Mpc_per_deg[ii] = (1.+zgrid[ii]) * 3600./arcsec_per_Mpc 

   integrand[ii] = dw2dz * qso_spaceden_z(zgrid[ii], gmag_0, $
                                         gmag) * $
                   Mpc_per_deg[ii]^2 / $
                   sqrt(0.267*(1.+zgrid[ii])^3 +0.733)

endfor
;stop

lowzcut= where(zgrid LE 2.1999)
hizcut = where(zgrid GT 2.2)

if lowzcut NE [-1] then $
   integral1 = int_tabulated(zgrid[lowzcut], integrand[lowzcut]) else $
      integral1 = 0.
if hizcut NE [-1] then $
   integral2 = int_tabulated(zgrid[hizcut], integrand[hizcut]) else $
      integral2 = 0

return, integral1 + integral2
end
