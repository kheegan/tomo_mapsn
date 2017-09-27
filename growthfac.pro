;+
; Returns cosmological growth factor D1(a). See Eq 7.77 in Dodelson
;
; Usage:
;   D_1 = growthfac(a, Omega_M = Omega_M, Omega_Lambda=Omega_Lambda,
;                   h=h)
; Input parameters:
;     a            - 1/(1+z)
;     Omega_M     - Omega_baryons+Omega_darkmatter
;     Omega_Lambda 
;     h 
;
; On initial call, sets up a lookup table and subsequently
; interpolates from this table. 
;
;-
function E_a, a, Omega_M, Omega_Lambda
;; Returns (H(a)/H_0)^2

return, Omega_M / a^3 + Omega_Lambda
end

function growthfac, a1, Omega_M = Omega_M1, Omega_Lambda=Omega_Lambda1,$
                   h=h1

common d1_block, a_grow, d1_grow, Omega_M_bl, Omega_Lambda_bl, h_bl

if keyword_set(Omega_M1) then Omega_M=Omega_M1 else $
   if keyword_set(Omega_M_bl) then Omega_M = Omega_M_bl
if keyword_set(Omega_Lambda1) then Omega_Lambda = Omega_Lambda1 else $
   if keyword_set(Omega_Lambda_bl) then Omega_Lambda = Omega_Lambda_bl
if keyword_set(h1) then h = h1 else $
   if keyword_set(h_bl) then h = h_bl

if not keyword_set(Omega_M) AND not keyword_set(Omega_M_bl) then $
   Omega_M = 0.3d
if not keyword_set(Omega_Lambda) AND not keyword_set(Omega_Lambda_bl) $
then Omega_Lambda = 0.7d
if not keyword_set(h) AND not keyword_set(h_bl) then h= 0.7d

if not keyword_set(Omega_M_bl) then Omega_M_bl = 0.
if not keyword_set(Omega_Lambda_bl) then Omega_Lambda_bl=0.
if not keyword_set(h_bl) then h_bl = 0.

;; Generate lookup table if not already existent, or if cosmo
;; parameters have changed
if not keyword_set(d1_grow) or Omega_M NE Omega_M_bl OR $
   Omega_Lambda NE Omega_Lambda_bl OR h NE h_bl then begin
   
   Omega_Lambda_bl = Omega_Lambda
   Omega_M_bl     = Omega_M

   delta_a = 0.02d

   ngrid = 48
   agrow = reverse(1.d -  delta_a * dindgen(ngrid))
   dl_grow = dblarr(ngrid)

   ;; This is to set up the integrand for the first bin
   agrid = 0.001d + 0.001*dindgen(20) 
   E_a_grid = E_a(agrid,Omega_M, Omega_Lambda )

   for igrid = 0, ngrid-1 do begin
      
      atmp = (agrow[igrid])[0]
      agrid = [agrid, atmp]
      E_a_grid = [E_a_grid, E_a(atmp, Omega_M, Omega_Lambda)]

      integrand = 1.d/ (agrid * sqrt(E_a_grid))^3
      prefactor = 2.5d * Omega_M * sqrt(E_a(atmp, Omega_M, Omega_Lambda))

      dl_grow[igrid] = prefactor * int_tabulated(agrid, integrand, /double)

   endfor

endif

return, interpol(dl_grow, agrow, a1)
end
