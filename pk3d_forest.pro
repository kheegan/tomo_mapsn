;+
; NAME:
;    PK3D_FOREST
;
; PURPOSE:
;    Return the 3D Lya forest flux power spectrum, <delta_f(k)
;    deltaf*(k)>. k is in units of h Mpc^-1
;
;    Uses a simple analytic form from McQuinn & White 2010, that gives
;    the forest power as a biased version of the underlying matter
;    P(k):
;
;    P_F,3d(k) = b^2 (1+beta * mu^2)^2 * P_dm(k) * exp(-k_par^2/k_D^2)
;
; CATEGORY:
;    Function
;
; CALLING SEQUENCE:
;    pk3s = pk3d_forest(k, mu, bias=bias1, beta=beta1, z=z1, k_D=k_D1)
;
; INPUTS:
;    k        - Wavenumber (in h Mpc^{-1}
;    mu       - Angle w.r.t. line-of-sight, mu = k_par/k
;
; OPTIONAL INPUTS:
;    bias     - Lya forest bias
;    beta     - Redshift distortion parameter
;    z        - redshift
;    k_D      - Wavenumber of forest thermal cutoff. Must be in h Mpc^-1
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
;   T_bbks(), pk_linear()
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;    KG Lee 16/08/2013  - Original code
;    KG Lee 27/09/2017  - Folded BBKS functions

function T_bbks, k, Omega_M = Omega_M, Omega_b=Omega_b, h=h
;; BBKS transfer function
if not keyword_set(Omega_M) then Omega_M = 0.26d
if not keyword_set(Omega_b) then Omega_b = 0.04d
if not keyword_set(h)  then h = 0.7d

;; From Sugiyama 1995. Gamma = Omega_M * h in the absence of baryons
Gamma = Omega_M * h * exp(-Omega_b * (1.d + ( sqrt(2.d*h)/Omega_M)))

q = k / Gamma

prefac = alog(1.d + 2.34d*q) / 2.34d  /q

poly = 1.d + 3.89d*q + (16.2d * q)^2 + (5.47d * q)^3 + (6.71d * q)^4 

return, prefac * poly^(-0.25)
end

function pk_linear, k, z
;; This is the BBKS P(k) linear
if not keyword_set(z) then z = 2.3

ns = 0.965d
c_over_H0 = 2997.92458d
;; This comes from Delta_R^2 = 2.41e-9 from WMAP9
delta_H = 5.06d-5  

prefactor = 2.d * !dpi^2  * delta_H^2 * c_over_H0^3

pk = prefactor * (k * c_over_H0)^ns * (T_bbks(k))^2 * $
     (growthfac(1./(1.+z)) / growthfac(1.))^2

return, pk
end


;-
function pk3d_forest, k, mu, bias=bias1, beta=beta1, z=z1, k_D=k_D1

if keyword_set(z1) then z = z1 else z = 2.25d
;; These bias and beta are the central values from Slosar et al 2011
;; at z=2.25
;; For bias, we use the central value of Slosar but interpolate the
;; redshift dependence from Table 1 in McQuinn & White 2010
if keyword_set(bias1) then bias = bias1 else begin
   bias_tab = [0.12, 0.18, 0.27, 0.37, 0.55]
   z_tab = [2.0, 2.5, 3., 3.5, 4.]
   bias = 0.2d/0.15d * interpol(bias_tab, z_tab, z) 
endelse
if keyword_set(beta1) then beta = beta1 else beta = 0.8

;; This value is obtained from McQuinn & White 2010. k_D=0.08 s
;; km^{-1}. The redshift-dependent factor converts to h Mpc^-1 
if keyword_set(k_D1) then k_D = k_D1 else $
   k_D = 0.08d * 94.d * sqrt(1.d + z) / sqrt(3.25d) 


k_par = mu * k

pkf = bias^2 * (1.d + beta*mu^2)^2 * pk_linear(k, z) * $
      exp(-(k_par/k_D)^2)

return, pkf
end
