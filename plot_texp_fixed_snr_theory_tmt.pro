;; Reproduce scaling seen eps_scale.ps, using theory residuals

;;; Same as plot_texp_fixed_snr_theory.pro, but with relabeled with
;;; exposure times divided by 13.4 to simulate TMT


readcol, 'eps_scale_factor.txt', nlos, texp, scale, scale2, scale3

eps1 = scale * 67./sqrt(nlos)
eps2 = scale2 * 67./sqrt(nlos)
eps3 = scale3 * 67./sqrt(nlos)

texp_a1 = 20. * eps1^(-1.6)

set_plot, 'ps'
device, file='texp_fixed_snr_theory_tmt.ps', /encap, /color, xsize=7.3,ysize=5., $
        /inches

!p.font=0
texp = texp
plot, eps1, texp/13.4,  xran=[0.5, 5.5], yran=[0, 4.], ysty=9,$
      charsize=1.4, charthick=3,  ytit=textoidl('Estimated WFOS Exposure Time (hrs)'),  $
      xtit=textoidl('3D Tomographic Resolution (h^{-1} Mpc)'), $
       posit= [0.11, 0.12, 0.865, 0.95], $
      /norm, thick=3, psym=5,xsty=1

oplot, eps2, texp/13.4, color=djs_icolor('red'), psym=4
oplot, eps3, texp/13.4, color=djs_icolor('blue'), psym=6

nlos_tab =[112, 225, 359, 504, 657, 971, 1289,  1604, 1914, 2218, 2805, 3367]
texp_tab = [2.,  3,   4,   5. , 6.,  8.,  10.,   12.,  14.,  16,   20., 24]/13.4

nlos_tickv = [500, 1000, 1500, 2000, 2500, 3000,3500,4000, 5000, 6000, 7000]
texp_tickvals = interpol(texp_tab, nlos_tab, nlos_tickv)

axis, yaxis=1, ytickv = texp_tickvals, ytickname = strtrim(nlos_tickv,2), $
  ystyle=1, charsize=1.3, ytit=textoidl('n_{los} (per sq deg)'), $
  ythick=1, charthick=3, yticks = 10

;oplot, eps1, texp_a1, linesty=2, thick=3

legend, [textoidl('S/N_{map} = 2.0'),textoidl('S/N_{map} = 2.5'), $
         textoidl('S/N_{map} = 3.0')], linesty=0, $
        color=[djs_icolor('red'),djs_icolor('black'), djs_icolor('blue')], $
        charsize=1.2, charthick=2, thick=2, /top,/right_legend, $
        delimiter=' '

 ;;--------------------------------------------------------------------

beta = 0.5
bias = 0.207

nlos_tab = [112, 225, 359, 504, 657, 971, 1289,  1604, 1914, 2218, 2805];, 3367]
texp_tab = [2.,  3,   4,   5. , 6.,  8.,  10.,   12.,  14.,  16,   20.];, 24]

n_texp = n_elements(texp_tab)

epsilon = [0.8, 1.0, 1.2, 1.5,1.7, 2.0,2.2, 2.5, 3.,3.5, 4.,4.5,5.,6.]
;epsilon = [1.,1.6, 2.5, 4., 5.5]
n_epsilon = n_elements(epsilon)

snrgrid = fltarr(n_epsilon, n_texp)

;; Generate grid to interpolate over
for ieps = 0, n_epsilon - 1 do begin
   eps_tmp = (epsilon[ieps])[0]
   print, 'eps = ', eps_tmp
   for it = 0, n_texp - 1 do begin
      texp_tmp = (texp_tab[it])[0]
      print, '     texp = ', texp_tmp
      snrtmp = sigmaf_3d_theory(eps_tmp, beta=beta, bias=bias,/gaussian)/ $
               residvar_notwiener(eps_tmp, t=texp_tmp, beta=beta, $
                               bias=bias, Lpar=0.2)
      snrgrid[ieps, it] = sqrt(snrtmp)
   endfor

endfor
print, 'Finished generating grid'

snr_out = [2.0, 2.5, 3.]
snrcolors = [djs_icolor('red'), djs_icolor('black'), djs_icolor('blue')]
n_snr = 3

for isnr=0, n_snr-1 do begin
   snr = (snr_out[isnr])
   texp_vec = fltarr(n_epsilon)
   for ieps = 0, n_epsilon-1 do begin
      snr_fixed_eps = reform(snrgrid[ieps,*])
      texp_vec[ieps] = interpol(texp_tab, snr_fixed_eps, snr)
   endfor
   
   oplot, epsilon, texp_vec/13.4, linesty=2, color=snrcolors[isnr], $
          thick=3
   
endfor

device, /close
set_plot, 'x'

end
