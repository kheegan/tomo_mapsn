pro gen_snrlist_magcut, zlow, zhi, maglim, area, outfile=outfile, $
                        dv_pix=dv_pix, t_exp=t_exp, silent=silent, $
                        nlos_out=nlos_vec_out, mag=mag_vec_out, snr=snr_vec, $
                        nlos_lbg = nlos_vec_lbg_out, create_block=create_block
;+ 
; Generate list of spectral SNR and number of LOS, assuming some
; redshift/magnitude cuts, and area. Scaled to the 11.64 sq deg area
; of Martin's 250 Mpc/h box. SNR assume VIMOS HR-Blue with 4
; hrs exposure, 1" seeing, airmass 1.3, 3 nights from new
; moon (estimated scaling with t_exp and B)
;
; Usage:
;   gen_snrlist_magcut, zlow, zhi, maglim, area, outfile=outfile, 
;       deltav_pix=deltav_pix, t_exp=t_exp
;  
; Inputs:
;   zlow, zhi   - Lower and upper redshift limits, respectively
;   maglim      - B-magnitude cut
;   area        - Area coverage on the sky, in sq degrees
;   outfile     - Output filename
;   dv_pix      - Delta(v) of output SNR. Default corresponds to
;                 ~11 km/s, pixel width of 250 Mpc/h sims 
;   t_exp       - Telescope exposure time, in hours. Default to 4 hrs 
;
; Keyword options:
;   silent      - When set, suppresses verbosity. Pointless unless
;                 outfile is set. 
;   mag         - Outputs the vector of source magnitudes
;   snr         - Outputs the vector of source SNR
;   nlos_out        - Outputs vector of sightline densities
;   nlos_lbg    - Outputs vector LBG densities
;   create_block- Force evaluation of common block. If set to the
;                 limiting magnitude (i.e. a float) then this is used
;                 as upper limit to magnitude grid
;   
;
;  Uses Reddy et al 2009 LBG luminosity function, and
;  Palanque-Delabrouille 2013's quasar luminosity functions.  
;
;  K.G. Lee  1/9/2013 - Fixed common block to not rerun when maglim is
;                       changed. Results may change compared to
;                       previous iterations because magnitude bins are
;                       now fixed.
;-

  common snrmag_block, maglim_bl, area_bl, mag_vec, nlos_vec, zlow_bl, zhi_bl, $
     nlos_vec_lbg


if not keyword_set(dv_pix) then dv_pix = 23500./2048.
dv_vimos = 52. ; velocity width per VIMOS HR-Blue pixel. Necessary because we're using the VIMOS ETC for SNR evaluation.
pixratio = dv_vimos / dv_pix

if not keyword_set(area) then area = 11.64
if not keyword_set(t_exp) then t_exp=4.


if not keyword_set(mag_vec) then create_block=1 else $
   if area_bl NE area OR zlow NE zlow_bl $
   OR zhi NE zhi_bl then create_block=1

if keyword_set(create_block) then begin
;; Generate magnitude vector in Delta(mag) = 0.5 down to 18th. 
   if size(create_block,/type) EQ 4 OR size(create_block,/type) EQ 5 then $
      mag_upper = create_block else mag_upper = 25.1
   nbins_mag = ceil((mag_upper - 19.)/0.1) + 2 
   magbins= reverse(mag_upper - 0.1 * findgen(nbins_mag))
   nlos_vec = lonarr(nbins_mag-1)
   nlos_vec_lbg = nlos_vec
   snr_vec  = fltarr(nbins_mag-1)
   mag_vec = fltarr(nbins_mag-1) 
   
   nlos_total = 0
   
   nlos_gal_cumul = lbg_area_den(zlow, zhi, gmag=(magbins[0])[0])
   

nlos_gal_cumul = lbg_area_den(zlow, zhi, gmag=(magbins[0])[0])

;; Loop over mag bins
   for ii = 0, nbins_mag -2  do begin
      mag_low = (magbins[ii])[0]
      mag_hi  = (magbins[ii+1])[0]
      
      nlos_gal_cumul_tmp = lbg_area_den(zlow, zhi, gmag= mag_hi)
      
      nlos_gal = round((nlos_gal_cumul_tmp - nlos_gal_cumul)*area)
      nlos_gal_cumul = nlos_gal_cumul_tmp
      
      nlos_qso = round(qso_area_den(zlow, zhi, gmag=[mag_low, mag_hi]) * area)
      
      nlos_vec[ii] = nlos_gal+nlos_qso
      nlos_vec_lbg[ii] = nlos_gal
      mag_vec[ii] = mag_hi
      
      nlos_total += nlos_gal + nlos_qso
      
   endfor

   zhi_bl = zhi
   zlow_bl = zlow
   area_bl = area

endif

cutmag = where(mag_vec LE maglim)

snr_vec = 10.^(0.5*(alog10(t_exp[cutmag]/4.) - 0.8*(mag_vec[cutmag] - 25.3)))/ $
          sqrt(pixratio)
mag_vec_out = mag_vec[cutmag]
nlos_vec_out = nlos_vec[cutmag]
nlos_vec_lbg_out = nlos_vec_lbg[cutmag]

if not keyword_set(silent) then $
   print, [transpose(mag_vec[cutmag]), transpose(nlos_vec[cutmag]), $
           transpose(snr_vec)]

if keyword_set(outfile) then writecol, outfile, $
                                       snr_vec, nlos_vec[cutmag], mag_vec[cutmag], $
                                       fmt='(f7.2,2x,i6,2x,f5.2)'

if not keyword_set(silent) then begin
   print, 'Total sightlines within the area:', fix(total(nlos_vec[cutmag]))
   print, 'Area density = ', float(total(nlos_vec[cutmag]))/area, ' per sq deg'
endif

end
