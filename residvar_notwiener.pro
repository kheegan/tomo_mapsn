function corrfunc_k, k, mu, Lperp=Lperp, Lpar=Lpar, z=z
;; This is the Fourier transform of the C_MD and C_DD 
;;
;; Does NOT include the sigma^2 prefactor in from of C_MD

k_par = k * mu
k_perp = sqrt(replicate(k,n_elements(k_par))^2 - k_par^2)


prefactor = sqrt(!dpi * double(Lpar)^2) * sqrt(!dpi * double(Lperp)^2)

gauss_perp = exp(- !dpi^2 * k_perp^2 * Lperp^2)
gauss_par  = exp(- !dpi^2 * k_par^2  * Lpar^2)

return, prefactor*gauss_perp*gauss_par
end

function resid_k, k, mu, z=z1, t_exp=t_exp, mag_lim=mag_lim, $
                  bias=bias, beta=beta,Lperp=Lperp, Lpar=Lpar, $
                  sigf_sq = sigf_sq, mag_50=mag_50
; Expression for residual of Wiener filter, in k-space
;
; k is in units of h^-1 Mpc, but k_par is converted to s km^-1 for
; input to pk1d_forest/
;
; mag_50 is the magnitude above which nlos gets cut in half

if not keyword_set(t_exp) then t_exp = 4.
if not keyword_set(z1) then z=2.25 else z=z1
if not keyword_set(mag_lim) then $
   mag_lim = 1.25*alog10(t_exp/16.) + 24.75

k_par = (mu * k)/(94.d * sqrt(1.d + z) / sqrt(3.25d))

pk3d = pk3d_forest(k, mu, z=z, bias=bias, beta=beta)
neff = neff_weighted(abs(t_exp), mag_lim=mag_lim, z=z, k_par=k_par, $
                     bias=bias, beta=beta,nlos_total=nlos_total, $
                    mag_50=mag_50)
if t_exp LT 0 then neff = nlos_total /(comdis(z, 0.3, 0.7)*2998./!RADEG)^2
pknoise = pk1d_forest(k_par,z=z) / (94.d * sqrt(1.d + z) / sqrt(3.25d)) $
          /  neff

; This is the FT of the Gaussian correlation function used in the
; reconstruction 
Corr_k = corrfunc_k(k, mu, Lperp=Lperp, Lpar=Lpar, z=z) * $
         sigf_sq

term1 = pknoise^2 * pk3d / (pk3d + pknoise)^2
term2 = pk3d^2 * pknoise / (pk3d + pknoise)^2


return, term1+term2
end

function residvar_notwiener, R, z=z1, t=t_exp, mag_lim=mag_lim, bias=bias, $
                          beta=beta, Lperp=Lperp, Lpar=Lpar, mag_50=mag_50
;R = 2. ; Mpc/h
if not keyword_set(z1) then z = 2.25 else z=z1

sigf_sq = (reform(sigmaf_3d_theory(R, bias=bias, beta=beta,z=z)))[0]

if not keyword_set(Lpar) then Lpar = 0.17d
if not keyword_set(Lperp) then Lperp = 2.95d * (t_exp/9.d)^(-0.625)

ngrid_k = 21
ngrid_mu = 17

kgrid  = 0.25*dindgen(ngrid_k)+0.01
mugrid = dindgen(ngrid_mu) /(double(ngrid_mu)-0.992)+0.0001

integrand_k = dblarr(ngrid_k)

for ik=0, ngrid_k-1 do begin

   k_tmp = (kgrid[ik])[0]

   integrand_mu = resid_k(k_tmp, mugrid, t_exp=t_exp, mag_lim=mag_lim, $
                          bias=bias, beta=beta, Lperp=Lperp, Lpar=Lpar, $
                          sigf_sq=sigf_sq,mag_50=mag_50)*  k_tmp^2 * $
                  (exp(-0.5d *(k_tmp * R)^2))^2
                  ;(3. * (sin(k_tmp*R) - k_tmp*R* cos(k_tmp*R)) / $
                  ; (k_tmp*R)^3)^2

   ;print, integrand_mu
   
   ;stop
   integrand_k[ik] = int_tabulated(mugrid, integrand_mu, /double)

   ;print, integrand_k[ik]

endfor
sigmasq_R = int_tabulated(kgrid, integrand_k, /double) / 2./!dpi

return, sigmasq_R
end
