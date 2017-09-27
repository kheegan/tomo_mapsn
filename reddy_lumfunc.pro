FUNCTION REDDY_LUMFUNC, z, lum_ratio = lum_ratio, mag = mag, absM = M_lum $
                        , mstar_mag = mag_star

;; extent of validity including space-based data is: 
;; M_lum  = -18.2d which corresponds to lum_ratio = L/L_* = 0.088, 
;; mag = 27.3, and phi = 6.9e-3 Mpc^-3

;; Reddy used this cosmology and units of Mpc^-3. 
LIT_H = 0.70D
OMEGA0 = 0.30D
LAMBDA0 = 0.70D
W = -1.0D

;; Lum function of LBGs at z ~ 1.9-3.4 using ground + HDF data, from Table 7
;; of (Reddy et al. 2008, ApJS, 175, 48)
IF z GT 2.7 AND z LE 3.4 THEN BEGIN
   alpha = -1.57d ;; Lum function of LBGs at z ~ 3 using HDF data (Reddy et al.)
   M_star = -20.84d
   phi_star = 1.66d-3 ;; units h0.7^3 Mpc^(-3) mag^-1
ENDIF ELSE IF z GT 1.9 AND z LE 2.7 THEN BEGIN
   alpha = -1.60d
   M_star = -20.60d
   phi_star = 3.31d-3
ENDIF ELSE BEGIN 
   splog, 'Redshift outside range of validity for Reddy lumfunc'
   splog, 'PROCEED WITH CAUTION DUDE. We are blindly extrapolating!!'
   IF z GE 3.4 THEN BEGIN
      alpha = -1.57d ;; Lum function of LBGs at z ~ 3 using HDF data (Reddy et al.)
      M_star = -20.84d
      phi_star = 1.66d-3 ;; units h0.7^3 Mpc^(-3) mag^-1
   ENDIF ELSE IF z LE 1.9 THEN BEGIN
      alpha = -1.60d
      M_star = -20.60d
      phi_star = 3.31d-3
   ENDIF
ENDELSE

anow = 1.0d/(1.0d + z)
DM = distance_modulus(anow, lit_h, omega0, lambda0, w)
mag_star = M_star + DM + 2.5d*alog10(anow) 

;; apparent magnitude corresponding M_star
IF KEYWORD_SET(LUM_RATIO) AND KEYWORD_SET(mag) THEN $
   message, 'Cant specifiy both' $
ELSE IF KEYWORD_SET(LUM_RATIO) THEN BEGIN
   M_lum = M_star - 2.5d*alog10(lum_ratio)
   mag = M_lum + DM + 2.5d*alog10(anow)
ENDIF ELSE IF KEYWORD_SET(mag) THEN BEGIN
   M_lum = mag - DM - 2.5d*alog10(anow)
   lum_ratio = 10.0d^(0.4d*(M_star - M_lum))   
ENDIF

   igam_func = gamma(alpha+1.0d)*(1.0d - igamma(alpha+1.0d, lum_ratio)) 
   n_faint = phi_star*igam_func
   
   RETURN, n_faint
END
