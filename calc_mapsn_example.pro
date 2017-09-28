;; Example script to illustrate how to run IGM tomography map S/N code
;; from Lee+2014. This basically implements Eq 18 in that paper
;;
;; To run, '.run calc_mapsn_example' at terminal


;; Desired map resolution
epsilon_3d = 1.5

;; Parameters for 3D Lya forest power spectrum 
beta = 0.5     
bias=0.207

;; This is an example input file lists the magnitudes, number
;; of sightlines (per sq deg) and corresponding spectra S/N (per
;; angstrom) given the exposure time. 
;;
;; Generated with the assumptions of Lee+2014, i.e. luminosity
;; function and S/N vs t_exp assumed there
print, 'Calculating S/N by reading in sightline S/N distribution'
readcol, 'nlos_sn_g25.0_texp8.0.txt', gmag, nlos_in, snr_in, /silent

print, 'Sightline density (per sq deg) = ', total(nlos_in)

sigmaf_3d = sigmaf_3d_theory(epsilon_3d, beta=beta, bias=bias, /gaussian)

;; Lperp should be roughly the sightline separation (in Mpc/h)
Lperp = !pi/180. * 2998. * comdis(2.3, 0.31, 0.69)/sqrt(total(nlos_in))

residvar = residvar_notwiener(epsilon_3d,beta=beta, bias=bias, Lpar=1., $
                                 Lperp=Lperp,nlos_in=nlos_in,snr_in=snr_in)

print, 'S/N_map = ', sqrt(sigmaf_3d/residvar)

;; One can also use the internal scripts to
;; generates a sightline S/N distribution given a
;; limiting magnitude and exposure time. This adopts the assumptions
;; in Lee+2014 (which are a bit dodgy) 
;;
print, 'Now, carry out the same using internal scripts'
residvar2 = residvar_notwiener(epsilon_3d, t=8., mag_lim=25, $
                               beta=beta, bias=bias, Lpar=1.)
print, 'S/N_map = ', sqrt(sigmaf_3d/residvar2)
stop
end
