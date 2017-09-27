;+
; NAME:
;    PK1D_FOREST
;
;
; PURPOSE:
;    Compute theoretical 1D Lyman-alpha forest power spectrum, based
;    on a biased linear-theory model. k is units of s km^-1
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:
;
;
;
; OPTIONAL INPUTS:
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;    pk3d_forest()
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-
function pk1d_forest, k_par, z=z1, bias=bias1, beta=beta1, k_D=k_D

if not keyword_set(z1) then z =2.25 else z=z1

ngrid = 71
kmax = 0.6
dk = (alog10(kmax) + 3.) / float(ngrid-1) 
kperp = 10.d^(-4. + dindgen(ngrid) * dk) 

if n_elements(k_par) EQ 1 then begin

   
   k_mag = sqrt(kperp^2 + k_par^2)
   mu    = k_par/k_mag
   
   k_mag_mpc = k_mag * (94.d * sqrt((1.+z)/3.25))
   
   pk3d = pk3d_forest(k_mag_mpc, mu, z=z, bias=bias1, beta=beta1, $
                      k_D=k_D)
   
   integrand =  kperp * (pk3d * (94.d * sqrt((1.+z)/3.25))^3)

   pk1d = int_tabulated(kperp, integrand)/2.d/ !dpi
   
endif else begin
   npar = n_elements(k_par)
   pk1d = dblarr(npar)

   for ipar=0, npar-1 do begin
      k_par_tmp = (k_par[ipar])[0]
      
      k_mag = sqrt(kperp^2 + k_par_tmp^2)
      k_mag_mpc = k_mag * (94.d * sqrt((1.+z)/3.25))
      mu = k_par_tmp/ k_mag

      pk3d = pk3d_forest(k_mag_mpc, mu, z=z, bias=bias1, beta=beta1 $
                         , k_D=k_D)
      
      integrand =  kperp * (pk3d * (94.d * sqrt((1.+z)/3.25))^3)
      
      pk1d[ipar] = int_tabulated(kperp, integrand)/2.d/ !dpi
      
   endfor
   
endelse

return, pk1d
end
