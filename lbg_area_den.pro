function lbg_area_den, z1, z2, gmag=gmag
;; Evaluate area density of LBG (per sq deg) at redshifts between
;; z1 and z2, and brighter than some g-magnitude limit
;
;  KG Lee 11/08/2013: Added magnitude offset to account for the fact
;  that Reddy et al 2008 uses m_{GR} not g. We assume m_{GR}-g=-0.2

delta_z = 0.002

n_z = round((z2 - z1)/delta_z) + 1
zgrid = z1 + delta_z * findgen(n_z)

;; c/H_0/sqrt(Omega_matter) ... multiplying dz with this gives comoving
;; radial distance
dw2dz = 3.e5 / 71. ;/ sqrt(0.267) 

integrand = fltarr(n_z)
Mpc_per_deg = fltarr(n_z)

for ii = 0, n_z - 1 do begin

   arcsec_per_Mpc = zang(1000.,(zgrid[ii])[0], H0=71.,Omega_m = 0.267, $
                         k=0,/silent)
   
   ;; (1+z) is to convert from physical to comoving
   Mpc_per_deg[ii] = (1.+zgrid[ii]) * 3600./arcsec_per_Mpc 

   integrand[ii] = dw2dz * reddy_lumfunc(zgrid[ii], mag=gmag-0.2) * $
                   Mpc_per_deg[ii]^2 / $
                   sqrt(0.267*(1.+zgrid[ii])^3 +0.733)

endfor
;stop

integral = int_tabulated(zgrid, integrand)

return, integral
end
