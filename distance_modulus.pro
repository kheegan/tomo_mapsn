FUNCTION distance_modulus, a, h, omega0, lambda0, w

a = double(a)
lit_h = double(h)
omega0 = double(omega0)
lambda0 = double(lambda0)
w = double(w)

HORIZON_DM = 42.38410353        ; this is 5*log10(HORIZON/10pc)

D = dofa(a, omega0, lambda0, w)
; luminosity distance in Mpc (converted from Mpc/h)
D_L_Mpc = D/a/lit_h

DM = HORIZON_DM + 5.0*alog10(D_L_Mpc)
  


RETURN, DM
END







