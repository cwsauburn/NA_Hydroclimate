
 random_setallseed(36454749,10494848)

 url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"

 dir = url+"CESM-LE/lnd/postprocessed/"

 dirout = "/glade/work/kumar34/hydroclimate_ncc/"

 dirIn = "/glade/work/kumar34/hydroclimate_ncc/"


 ;dir = "/Projects/SERDP_regimeshifts/data/CESM-LE/lnd/postprocessed/"
 ;dirout = dir+"regions_clean/"
 
 nyrp = 30
 yrstr = (/1940,1970,2000,2030,2060/)
 yrend = yrstr + nyrp - 1

 ;region = (/"Southwest","Southeast","SouthernUSPlains","NorthernUSPlains","CanadianPlains"/)
  ;;region = (/"Southwest","Southeast", "CanadianPlains_rev1"/)
    region = (/"CanadianPlains_rev1"/)
    nrgs = dimsizes(region)

  ;latS = (/30.,25.,28.,43.,53./)
  ;latN = (/40.,35.,38.,49.,61./)
  ;lonW = (/239.,265.,255.,254.,242./)
  ;lonE = (/252.,282.,263.,264.,257./)

  ;latS = (/30.,25.,49./)
  ;latN = (/40.,35.,60./)
  ;lonW = (/239.,265.,248./)
  ;lonE = (/252.,282.,263./)

   latS = (/49./)
   latN = (/60./)
   lonW = (/248./)
   lonE = (/263./)




 do r = 0,nrgs-1
   print_table([/region(r)+"*****************************************************"/],"%s")

   filereg = dirIn+"H2OSOI_CESM-LE_1940-2089_30yr_MemoryENSOregressions_"+region(r)+".nc"
   fr = addfile(filereg,"r")
 
   ensi = ispan(1,40,1)
   nEns = dimsizes(ensi)
 
   periods = fr->periods
   nps = dimsizes(periods)
   ntim = nyrp - 1 ; -1 because by smoothing enso we lose one year
 
   soilm := fr->lvar
;   grabAnnual = ispan(2,nyrp*12-1,12) ; first Oct-Sep annual mean is in record "2"
   grabAnnual = ispan(14,nyrp*12-1,12) ; but because of enso missing data we skip the first year
   soilm := soilm(periods|:,ens|:,time|grabAnnual)

   enso := fr->enso
   enso := enso(periods|:,ens|:,months|grabAnnual)
   
   alfa = fr->coefs(2,:,0,:) ; soil memory coefficients for Memory+ENSO model
   beta = fr->coefs(2,:,1,:) ; enso coefficients for Memory+ENSO model
   
   errors := fr->error(2,:,:,:)
   errors := errors(periods|:,ens|:,time|grabAnnual)

   ; get the mean and stddev of the regression errors to generate a random white noise component
   erravg = dim_avg_n(errors,(/1,2/))
   errstdv = dim_stddev_n(errors,(/1,2/))

   ; check distribution of errors for normality
   ; 1. first generate a random normal distribution with same avg, stddev, and sample size
   ; 2. apply KS test to check if distribituion of errors is statistically different from the normal distribution
   do p = 0,nps-1
     p_errors = ndtooned(errors(p,:,:))
     randnorm = random_normal(erravg(p),errstdv(p),dimsizes(p_errors))
     kstest = kolsm2_n(p_errors,randnorm,0)
     if(kstest.gt.0.05) then
       print_table([/kstest,": errors are normally distributed"/],"%4.2f %s")
     else
       print_table([/kstest,": errors are NOT normally distributed"/],"%4.2f %s")
     end if
   end do
   
   ; generate nit*nEns time series of synthetic soilm from Memory+ENSO model using random initial conditions and error
   ; (nEns here represents the random sampling (with replacement) of the 40 ensemble members)
   nit = 1000
   ens_sample = new((/nit,nEns/),"integer")
   do b = 0,nit-1
     sample = toint(random_uniform(0,nEns,nEns))
     ens_sample(b,:) = sample
     delete(sample)
   end do

   soilm_synt := new((/nps,nit,nEns,ntim/),"float",soilm@_FillValue)
   do p = 0,nps-1
     randerr := random_normal(erravg(p),errstdv(p),(/nit,nEns,ntim/))
     randSM0 := random_normal(0,1,(/nit,nEns/))
     do i = 0,nit-1
       do e = 0,nEns-1
         em = ens_sample(i,e)
         sm0 = randSM0(i,e)
         do t = 0,ntim-1
           soilm_synt(p,i,e,t) = alfa(p,em)*sm0 + beta(p,em)*enso(p,em,t) + randerr(i,e,t)
           sm0 = (/soilm_synt(p,i,e,t)/)
         end do
       end do
     end do
   end do

   ; read in original soil moisture anomalies to get stddev that was used to standardize before regression 
   filesm = dir+"H2OSOI_CESM-LE_40em_192001-210012_12mormAnomsEnsMean.nc"
   fl = addfile(filesm,"r")
   date := fl->date
   sti = ind(date.eq.194001)
   endi = ind(date.eq.208912)
   date := date(sti:endi)
   smvar := fl->H2OSOI(sti:endi,:,{latS(r):latN(r)},{lonW(r):lonE(r)})
   smtsr = dim_avg_n_Wrap(smvar,(/2,3/))
   smstdv := dim_stddev_n_Wrap(smtsr,0)

   ; use sm stddev to unstandardize soil moisture anomalies (both original and synthetic)
   soilm = soilm*conform_dims(dimsizes(soilm),smstdv,1)
   smstdv_rsmp = new((/nit,nEns/),"float")
   do i = 0,nit-1
     ems = ens_sample(i,:)
     smstdv_rsmp(i,:) = smstdv(ems)
   end do
   soilm_synt = soilm_synt*conform_dims(dimsizes(soilm_synt),smstdv_rsmp,(/1,2/))

   ; read in ensemble mean+trend from original soil moisture file and reshape for # of periods
   smem := fl->ensmean(sti:endi,{latS(r):latN(r)},{lonW(r):lonE(r)})
   smemtsr = dim_avg_n_Wrap(smem,(/1,2/))
   smem2d := onedtond(smemtsr,(/nps,nyrp*12/))
   smem2d := smem2d(:,grabAnnual)

   ; use ensemble mean+trend to recover soil moisture monthly values
   soilm_w_trend = soilm + conform_dims(dimsizes(soilm),smem2d,(/0,2/))
   soilm_synt_w_trend = soilm_synt + conform_dims(dimsizes(soilm_synt),smem2d,(/0,3/))

   ; compute ensemble mean for the 1970-1999 period
   soilm_temp = soilm_w_trend(1,:,:)
   soilm_30yrmn = avg(soilm_temp)
   soilm_synt_temp = soilm_synt_w_trend(1,:,:,:)
   soilm_synt_30yrmn = avg(soilm_synt_temp)

   ; compute anomalies from ensemble mean
   soilm_trendanoms = soilm_w_trend - soilm_30yrmn
   soilm_synt_trendanoms = soilm_synt_w_trend - soilm_synt_30yrmn

   copy_VarCoords(soilm,soilm_trendanoms)
   soilm_trendanoms@long_name = "CESM Oct-Nov SOILM"
   copy_VarCoords(soilm,soilm_synt_trendanoms(:,0,:,:))
   soilm_synt_trendanoms!1 = "iter"
   soilm_synt_trendanoms&iter = ispan(1,1000,1)
   soilm_synt_trendanoms@long_name = "CESM synthetic Oct-Nov SOILM"
   
   ; save to netcdf
   outfile = dirout+"H2OSOI_CESM-LE_1940-2089_30yr_TimeSeries_"+region(r)+".nc"
   system("/bin/rm -f "+outfile)
   ncdf = addfile(outfile,"c")
   ncdf->cesm_soilm = soilm_trendanoms
   ncdf->cesm_synt_soilm = soilm_synt_trendanoms

 end do
 
