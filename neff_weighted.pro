;+
; NAME:
;    NEFF_WEIGHTED
;
; PURPOSE:
;    Calculate n_eff(k_parallel), effective area density of a LyaF
;    survey, Eq 13 in McQuinn & White 2011. In units of h^2 Mpc^-2 
;
;    k_par may be a vector but not t_exp or mag_lim
;
; CALLING SEQUENCE:
;    neff = neff_weighted(mag_lim=mag_lim, k_par=k_par, t_exp=t_exp,
;                         nlos_total=nlos_total,bias=bias,beta=beta)
;
; INPUTS:
;    mag_lim   - Limiting magnitude of survey. Implicitly sets number
;                of sightlines
;    k_par     - LOS k-mode in s km^-1. Optional.
;    t_exp     - Exposure time, based on VLT-VIMOS
;    bias
;    beta
;
; OPTIONAL INPUTS:
;    mag_50    - Below this magnitude, take n_los->0.5n_los. mag_50 <
;                mag_lim  
;    nlos_in   - Alternative input: vector of sightlines/deg^2 (associated
;                with snr_in vector)
;    snr_in    - Vector of pixel snr (per angstrom)
;
;
; OPTIONAL OUTPUTS:
;    nlos_total- Areal density of sightlines per sq deg
;
; COMMENTS:
;    Current version has the spectral noise power slightly wrong (not
;    taking into account mean LyaF absorption)
;
;-
function P_noise, snr_ang, z=z1

return, 0.8/ snr_ang^2 * ((1.+z1)/4.)^(-1.5) ;/exp(-taueff_evo(z1))^2 
end

function neff_weighted, t_exp, mag_lim=mag_lim, k_par=k_par, z=z1, $
                        nlos_total=nlos_total, bias=bias,beta=beta, $
                        mag_50=mag_50, nlos_in=nlos_in, snr_in=snr_in

common neff_block, snr, nlos, mag,t_exp_block, mag_block

if not keyword_set(z1) then z = 2.25 else z = z1
if not keyword_set(t_exp) and not keyword_set(nlos_in) then t_exp = 4.
if not keyword_set(mag_lim) and not keyword_set(nlos_in) then mag_lim = 1.25*alog10(texp / 16.)+ 24.75
if not keyword_set(k_par) then pk_los = 2.8 else $ ; P_F(k) at k=190 s/km 
   pk_los = pk1d_forest(k_par,z=z,bias=bias, beta=beta)

pk_los = pk_los / (94.* sqrt((1.+z)/3.25))

if not keyword_set(t_exp_block) then t_exp_block = -1
if not keyword_set(mag_block) then mag_block = -1

if keyword_set(nlos_in) AND keyword_set(snr_in) then begin
   nlos = nlos_in
   snr  = snr_in
endif else begin
   if (t_exp_block NE t_exp) OR (mag_lim NE mag_block) then begin
      gen_snrlist_magcut, 2.35, 2.7, mag_lim, 1., snr=snr, nlos_out=nlos, $
                          dv_pix=76.,/silent, t_exp=t_exp,mag=mag ;, create_block=mag_lim
      t_exp_block = t_exp
      mag_block = mag_lim
   endif
endelse

if keyword_set(mag_50) then nlos[where(mag GE mag_50)] = 0.5*nlos[where(mag GE mag_50)]

nlos_total = total(nlos)

nu_n = dblarr(n_elements(pk_los), n_elements(snr))

npk = n_elements(pk_los)
nnoise = n_elements(snr)
if size(pk_los,/dim) EQ 0 then pk_los = [pk_los]
pk_los = rebin(pk_los,npk,nnoise)
Pnoise = transpose(rebin(P_noise(snr,z=z),nnoise,npk))

nu_n = pk_los / (pk_los + Pnoise)

neff = nu_n # nlos 

area_1deg = (comdis(z, 0.3, 0.7)*2998./!RADEG)^2

return, neff/area_1deg
end
