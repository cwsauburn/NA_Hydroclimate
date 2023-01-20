;***************************************************************************************

 url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"


 dir = url+"GFDL-CM3-LE/"

; specify start and end dates for each 20-yr period
 nyrp = 30
 s0 = 194001
 wstep = nyrp*100
 nps = 5
 strdates = new((/nps/),"integer")
 do p = 0,nps-1
   strdates(p) = s0
   s0 = s0 + wstep
 end do
 enddates = strdates + wstep - 89
 ntimp = nyrp*12
 timeper = ispan(1,ntimp,1)
 yrstr = (/1940,1970,2000,2030,2060/)
 yrend = yrstr + nyrp - 1
 periods = (/19401969,19701999,20002029,20302059,20602089/)

;***************************************************************************************

; reading the ocean EOF timeseries
 ocnvar = "enso"
 fileocn = dir+"SST/GFDL-CM3-LE_sst_pc1_1940-2089_30yrperiods_20em.nc"
 fo = addfile(fileocn,"r")
 pacpc1 = fo->pacpc1_12mrm
 stdv = dim_stddev_n_Wrap(pacpc1,(/0,1/))
 pacpc1 = pacpc1 / conform_dims(dimsizes(pacpc1),stdv,2)
 delete(stdv)
 pacpc1!2 = "ens"

;***************************************************************************************

; reading the land variable at n and n-1yr timesteps, subsetting over region
 varname = "SOILM"
 
 filelnd = dir+varname+"/"+varname+"_GFDL-CM3-LE_20em_192001-210012_12mormAnomsEnsMean.nc"
 fl = addfile(filelnd,"r")
 
 date = fl->date
 sti = ind(date.eq.strdates(0))
 endi = ind(date.eq.enddates(nps-1))
 date := date(sti:endi)
 
 ensi = ispan(0,19,1)
 nEns = dimsizes(ensi)
 
;***************************
; region = "Southwest"
; latS = 30.
; latN = 40.
; lonW = 239.
; lonE = 252.
;***************************
 region = "Southeast"
 latS = 25.
 latN = 35.
 lonW = 265.
 lonE = 282.
;***************************
; region = "SouthernUSPlains"
; latS = 28.
; latN = 38. ;45.
; lonW = 255.
; lonE = 263.
;***************************
; region = "NorthernUSPlains"
; latS = 43.
; latN = 49.
; lonW = 254.
; lonE = 264.
;***************************
; region = "CanadianPlains_rev1"
; latS = 49.
; latN = 60.
; lonW = 248.
; lonE = 263.
;***************************

 lat = fl->lat({latS:latN})
 nlat = dimsizes(lat)
 lon = fl->lon({lonW:lonE})
 nlon = dimsizes(lon)
 
 lvar4d = fl->$varname$(sti:endi,:,{latS:latN},{lonW:lonE})
 lvarpy4d = fl->$varname$(sti-12:endi-12,:,{latS:latN},{lonW:lonE})

;***************************************************************************************

; compute spatially averaged time series
 lvar40e = dim_avg_n_Wrap(lvar4d,(/2,3/))
 lvarpy40e = dim_avg_n_Wrap(lvarpy4d,(/2,3/))

; standardize both variables
 lvstdv = dim_stddev_n_Wrap(lvar40e,0)
 lvar40e = lvar40e / conform_dims(dimsizes(lvar40e),lvstdv,1)

 lvpystdv = dim_stddev_n_Wrap(lvarpy40e,0)
 lvarpy40e = lvarpy40e / conform_dims(dimsizes(lvarpy40e),lvpystdv,1)

; reshape into 30yr chunks
 lvar = new((/nps,ntimp,nEns/),"float",lvar40e@_FillValue)
 lvarpy = lvar
 do p = 0,nps-1
   sti = ind(date.eq.strdates(p))
   endi = ind(date.eq.enddates(p))
   lvar(p,:,:) = lvar40e(sti:endi,:)
   lvarpy(p,:,:) = lvarpy40e(sti:endi,:)
 end do

 lvar!0 = "periods"
 lvar&periods = periods
 lvar!1 = "time"
 lvar&time = timeper
 lvar&ens = pacpc1&ens
 copy_VarMeta(lvar,lvarpy)
 lvarpy@long_name = "previous year "+lvar@long_name

 delete([/lvar4d,lvarpy4d,lvar40e,lvarpy40e,lvstdv,lvpystdv/])
 
;***************************************************************************************
 
; model lvar as a function of lvarpy and pacpc1, with crossvalidation

 regressions = (/"Memory","ENSO","Memory+ENSO"/)
 regs = ispan(1,3,1)
 regstring = oneDtostring(regressions)
 nregs = dimsizes(regressions)
  
 coefs = new((/nregs,nps,2,nEns/),"float",lvar@_FillValue)       ; regression coefs for Memory and ENSO
 lvarmod = new((/nregs,nps,ntimp,nEns/),"float",lvar@_FillValue) ; lvarmod = forecast of lvar
 
 error = lvarmod                                                 ; forecast error (lvarmod - lvar) 
 actskill = new((/nregs,nps/),"float",lvar@_FillValue)           ; correlation between lvarmod and lvar
 sig2noise = actskill                                            ; signal-to-noise ratio of lvarmod
 predskill = actskill                                            ; predicted skill (rho infinite)
 
 do p = 0,nps-1
   print_table([/tostring(yrstr(p))+"-"+tostring(yrend(p))/],"%s")
   do i = 0,nEns-1
     ensgrab = ind(ensi.ne.i)    ; leaves one out for validation
     lv2d = lvar(p,:,ensgrab)
     lv := ndtooned(lv2d)
     lvpy2d = lvarpy(p,:,ensgrab)
     lvpy := ndtooned(lvpy2d)
     pc2d = pacpc1(p,:,ensgrab)
     pc := ndtooned(pc2d)

     notmiss = ind(.not.(ismissing(pc)))
     pc := pc(notmiss)
     lvpy := lvpy(notmiss)
     lv := lv(notmiss)

     ; regress against Memory and forecast using remaining ensemble member
     b1 = escorc(lv,lvpy)
     coefs(0,p,0,i) = b1
     lvarmod(0,p,:,i) = b1*lvarpy(p,:,i)
     
     ; regress against ENSO and forecast using remaining ensemble member
     b2 = escorc(lv,pc)
     coefs(1,p,1,i) = b2
     lvarmod(1,p,:,i) = b2*pacpc1(p,:,i)
     
     ; regress against Memory+ENSO and forecast using remaining ensemble member
     X = (/lvpy,pc/)
     Xt = transpose(X)
     XtX = X # Xt
     Xty = X # lv
     B = inverse_matrix(XtX) # Xty
     coefs(2,p,:,i) = B
     lvarmod(2,p,:,i) = B(0)*lvarpy(p,:,i) + B(1)*pacpc1(p,:,i)

     delete([/ensgrab,lv2d,lv,lvpy2d,lvpy,pc2d,pc,b1,b2,X,Xt,XtX,Xty,B/])
   end do

   lv = ndtooned(lvar(p,:,:))

   do r = 0,nregs-1
     ; get error/noise
     error(r,p,:,:) = lvarmod(r,p,:,:) - lvar(p,:,:)
     error1d = ndtooned(error(r,p,:,:))
     lvmod = ndtooned(lvarmod(r,p,:,:))
     ; get actual skill, or the (anomaly) correlation between forecast and obs
     actskill(r,p) = escorc(lv,lvmod)
     ; get signal-to-noise ratio
     s2nsq = variance(lvmod) / variance(error1d)
     sig2noise(r,p) = sqrt(s2nsq)
     ; get predicted skill, or rho(infinity)
     predskill(r,p) = s2nsq / sqrt((s2nsq+1)*s2nsq)
     delete([/error1d,lvmod,s2nsq/])
   end do
   
   delete(lv)
 end do
 
 copy_VarCoords(lvar,lvarmod(0,:,:,:))
 lvarmod!0 = "regressions"
 lvarmod&regressions = regs
 lvarmod@regress = regstring
 
 copy_VarMeta(lvarmod(:,:,0,:),coefs(:,:,0,:))
 coefs!2 = "predictors"
 coefs&predictors = (/1,2/)
 
 copy_VarMeta(lvarmod,error)
 copy_VarCoords(lvarmod(:,:,0,0),actskill)
 actskill@regress = regstring
 copy_VarMeta(actskill,sig2noise)
 copy_VarMeta(actskill,predskill)
 
;***************************************************************************************
 dirout = "/glade/work/kumar34/hydroclimate_ncc/"

; save to netcdf
 print_table([/"saving regression outputs to netcdf..."/],"%s")
 outfile = dirout+varname+"_GFDL-CM3-LE_1940-2089_30yr_MemoryENSOregressions_"+region+".nc"
 system("/bin/rm -f "+outfile)
 ncdf = addfile(outfile,"c")

 ncdf->lvar = lvar
 ncdf->lvarpy = lvarpy
 ncdf->$ocnvar$ = pacpc1
 ncdf->coefs = coefs
 ncdf->lvarmod = lvarmod
 ncdf->error = error
 ncdf->actskill = actskill
 ncdf->sig2noise = sig2noise
 ncdf->predskill = predskill
