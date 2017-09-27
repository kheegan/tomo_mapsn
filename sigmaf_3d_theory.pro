function tophat_filter, k, R

th = 3. * (sin(k*R) - k*R* cos(k*R)) / (k*R)^3
return, th
end

function gaussian_filter, k, R

return, exp(-0.5d *(k * R)^2)
end

function sigmaf_3d_theory, R1, beta=beta, bias=bias,z=z, gaussian=gaussian 

Rgrid = R1
nR = n_elements(Rgrid)
sigmasq_R = dblarr(nR)

ngrid_k = 26
ngrid_mu = 26

kgrid  = 10.d^(-2.3 + 0.1*dindgen(ngrid_k))
mugrid = dindgen(ngrid_mu) /(ngrid_mu-1)

for iR = 0, nR-1 do begin
   R = (Rgrid[iR])[0]
   integrand_k = dblarr(ngrid_k)
   for ik=0, ngrid_k-1 do begin
       
      k_tmp = (kgrid[ik])[0]

      if keyword_set(gaussian) then filter = gaussian_filter(k_tmp, R) $
      else filter = tophat_filter(k_tmp, R)
      
      integrand_mu =  pk3d_forest(k_tmp, mugrid,bias=bias, $
                                                     beta=beta,z=z,k_D=k_D) * $
                     k_tmp^2 * filter^2
      
      integrand_k[ik] = int_tabulated(mugrid, integrand_mu, /double)
      
   endfor
   
   sigmasq_R[iR] = int_tabulated(kgrid, integrand_k, /double) / 2./!dpi

endfor

return, sigmasq_R
end
