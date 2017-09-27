function qso_lumfunc_bossmmt, Mg, z
;+
;
; Return best-fit model QSO luminosity function (in units of
; /Mpc^3/gmag), for an input Mg and z
;
; Usage: Phi = qsolumfunc_bossmmt(Mg,z)
;-

Mg_zp = -26.36
zp = 2.2
log_Phi_star = -5.89
if z LT 2.2 then begin
   alpha = -3.5
   beta = -1.43
   k1 = 0.03
   k2 = -0.34
endif else begin
   alpha = -3.19
   beta = -1.17
   k1 = -0.35
   k2 = -0.02
endelse

Mg_z = Mg_zp - 2.5*(k1 * (z - zp) + k2 * (z - zp)^2)

exp1 = 0.4 * (alpha + 1)*(Mg - Mg_z)
exp2 = 0.4 * (beta + 1) *(Mg - Mg_z)

return, 10.^log_Phi_star / (10.^exp1 + 10.^exp2)
end
