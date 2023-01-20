
  data1 = asciiread("./Tropical_IndoPacific_30S_30N_PC1_1950_2014_EE_Earth.txt", (/780, 2/), "double")
  ;data1 = asciiread("./Tropical_30S_30N_PC1_1950_2014_EE_Earth.txt", (/780, 2/), "double")
  data1@_FillValue = -9999.0


  enso = tofloat(data1(:, 1))
  print(enso(0:10))

  ;exit
 
  ;***********Reading the soil moisture data


  st1_h = (1950-1850)*12
  end1_h = (2015 - 1850) *12 -1 

  inDirLH = "/glade/p/cgd/tss/people/oleson/CLM_LAND_ONLY_RELEASE/CLM5/clm50_r270_1deg_GSWP3V1_iso_newpopd_hist/lnd/proc/tseries/month_1"

  ;fl = addfile(inDirLH+"/clm50_r270_1deg_GSWP3V1_iso_newpopd_hist.clm2.h0.H2OSOI.185001-201412.nc", "r")

  ;********************
   fm = addfile("/glade/work/kumar34/CESM2-LE/H2OSOI/b.e21.BHIST-SSP370.f09_g17.LE2-1001.001.clm2.h0.H2OSOI.185001-210012.nc", "r")

   lm = fm->landmask
   area = fm->area
   landfrac = fm->landfrac

   arwt = area
   arwt = area * landfrac

   lat = fm->lat
   lon = fm->lon


    slat = (/37.5/)
    nlat = (/42.75/)

    elon = (/360.0-91.25/)
    wlon = (/360.0-87.75/)

    dzsoi_3D = fm->DZSOI(0:8, {slat:nlat}, {elon:wlon})

    dzsoi = wgt_areaave2(dzsoi_3D, arwt({slat:nlat}, {elon:wlon}), 0)

    print(sum(dzsoi))
   
    dzsoi(8) = 1.0 - sum(dzsoi(0:7))

    print(dzsoi)

    print(sum(dzsoi))

   ;********************************************
   ;*****LAT LON FOR NORTH AMERICA*************


    latS = 15.
    latN = 72.
    lonW = 192.
    lonE = 308.


 ;***************************
    fl = addfile(inDirLH+"/clm50_r270_1deg_GSWP3V1_iso_newpopd_hist.clm2.h0.H2OSOI.185001-201412.nc", "r")

    H2OSOI_hist_temp = fl->H2OSOI(st1_h:end1_h, 0:8, {latS:latN}, {lonW:lonE})

    H2OSOI_hist = dim_avg_wgt_n_Wrap(H2OSOI_hist_temp, dzsoi, 1, 1)

    H2OSOI_hist_climo = clmMonTLL(H2OSOI_hist(12:371, :, :))

    H2OSOI_ANOM = calcMonAnomTLL(H2OSOI_hist, H2OSOI_hist_climo)

    printVarSummary(H2OSOI_ANOM)

  ;*************************************
  ;*******Caculate 12 month running mean and deseasonalize the data****************
   
  var1 = H2OSOI_ANOM

  var1 = runave_n_Wrap(var1, 12, 0, 0)


  delete(H2OSOI_ANOM)

  nyr = floattointeger(dimsizes(var1&time)/12)      

  var1_dtrend = var1

  temp1 = var1(0:nyr-1, :, :)

  do i = 0, 11
     do j = 0, nyr-1
        temp1(j, :, :) = var1(j*12+i, :, :)
     end do

     temp1_dtrend = temp1
     temp1_dtrend = dtrend_msg_n(ispan(1, nyr, 1), temp1, False, False, 0)
     do j = 0, nyr-1
        var1_dtrend(j*12+i, :, :) = (/temp1_dtrend(j, :, :)/)
     end do

  end do


;**********************

;**Use  deterned data*********************
 
   var1 = var1_dtrend
  
 ;**********************************************
  enso_std_num = dim_stddev_n_Wrap(enso, 0)

  enso_std = enso
  enso_std = enso / enso_std_num 

 ;************************************************

 ;standardize the soil moisture

 sm_std = dim_stddev_n_Wrap(var1, 0)

 var1 = (/(var1 /  conform_dims(dimsizes(var1), sm_std, (/1, 2/)))/)

 
 ;*****************************************************


  var1py = var1
  var1py = var1@_FillValue
  var1py(12:nyr*12-1, :, :) = var1(0:(nyr-1)*12-1, :, :)     ; one year prior soil moisture

  lvargrid = dim_sum_n_Wrap(var1, 0)
  
 ;*********************************************************

  var1ModM = var1
  var1ModM = var1@_FillValue

  var1ModE = var1ModM

  var1ModME = var1ModM
;exit

 ;***************************************************

  lat := var1&lat
  nlat := dimsizes(lat)
 
  lon := var1&lon
  nlon := dimsizes(lon)

  time := var1&time
  ntime = dimsizes(time) 

  timei = ispan(0, dimsizes(time)-1, 1) 
  timei@_FillValue = -99

  regressions = (/"Memory","ENSO","Memory+ENSO"/)
  regs = ispan(1,3,1)
  regstring = oneDtostring(regressions)
  nregs = dimsizes(regressions)


  nps = 13     ; five 13 years period - 65 years
     
 
  coefs = conform_dims((/nregs,nps, 2,nlat,nlon/), var1(0, :, :), (/3, 4/)) 

  coefs = coefs@_FillValue

  actskill = coefs(:, 0, 0, :, :)
  sig2noise = actskill
  predskill = actskill

  do j = 0, nlat -1
     do i = 0, nlon -1
         if(.not.ismissing(lvargrid(j,i))) then ; only land points
                do p = 0, nps-1
                  st1 = toint(p*5*12) 
                  end1 = st1 + 5*12 - 1      ; leave the 13 years data for verfication
                  time1 = time
                  ;print(avg(time1))
                  time1(st1:end1) = time@_FillValue
                  ;print(avg(time1))
                  timegrab = ind(.not.ismissing(time1))                 
 
                  lv := var1(timegrab, j, i)
                  lvpy := var1py(timegrab, j, i)
                  pc := enso_std(timegrab)                  
                  
                  ;exit

                  notmiss1 := ind(.not.ismissing(lv))
                  lv := lv(notmiss1)
                  lvpy := lvpy(notmiss1)
                  pc := pc(notmiss1)                  


                  notmiss2 := ind(.not.ismissing(lvpy))
                  lv := lv(notmiss2)
                  lvpy := lvpy(notmiss2)
                  pc := pc(notmiss2)
  
                  notmiss3 := ind(.not.ismissing(pc))
                  lv := lv(notmiss3)
                  lvpy := lvpy(notmiss3)
                  pc := pc(notmiss3)

                  b1 = escorc(lv,lvpy)
                  coefs(0, p, 0, j, i) = b1
                  var1ModM(st1:end1, j, i) = b1 * var1py(st1:end1, j, i)                  



                 b2 = escorc(lv, pc)
                 coefs(1, p, 1, j, i) = b2
                 var1ModE(st1:end1, j, i) = b2 * enso_std(st1:end1) 


                 X = (/lvpy,pc/)
                 Xt = transpose(X)
                 XtX = X # Xt
                 Xty = X # lv
                 B = inverse_matrix(XtX) # Xty
                 coefs(2,p,:,j,i) = B
                 var1ModME(st1:end1, j, i) = B(0)* var1py(st1:end1, j, i) + B(1) * enso_std(st1:end1)

                 delete([/lv, lvpy, pc, b1, b2, X, Xt, XtX, Xty, B/])
                end do
            
            end if
          end do
     end do

 

   error1 = var1
   error1 = var1 - var1ModM
   error2 = var1
   error2 = var1 - var1ModE

   error3 = var1
   error3 = var1 - var1ModME  


   actskill(0, :, :) = (/escorc_n(var1, var1ModM, 0, 0)/) 
   actskill(1, :, :) = (/escorc_n(var1, var1ModE, 0, 0)/)
   actskill(2, :, :) = (/escorc_n(var1, var1ModME, 0, 0)/)

   
   s2nsq = dim_variance_n_Wrap(var1ModM, 0)
   s2nsq = s2nsq / dim_variance_n_Wrap(error1, 0) 

   sig2noise(0, :, :) = (/sqrt(s2nsq)/)
   onemat = sig2noise(0, :, :)
   onemat = 1.0
   predskill(0, :, :) = (/(s2nsq / sqrt((s2nsq + onemat) * s2nsq))/)

   s2nsq := dim_variance_n_Wrap(var1ModE, 0)
   s2nsq := s2nsq / dim_variance_n_Wrap(error2, 0)  

   sig2noise(1, :, :) = (/sqrt(s2nsq)/)
   predskill(1, :, :) = (/(s2nsq / sqrt((s2nsq + onemat) * s2nsq))/)


   s2nsq := dim_variance_n_Wrap(var1ModME, 0)
   s2nsq := s2nsq / dim_variance_n_Wrap(error3, 0)

   sig2noise(2, :, :) = (/sqrt(s2nsq)/)
   predskill(2, :, :) = (/(s2nsq / sqrt((s2nsq + onemat) * s2nsq))/)


   coefs_ens = dim_avg_n_Wrap(coefs, 1)
   coefs_ens!0 = "regress"
   coefs_ens&regress = regs
   coefs_ens@regressions = regstring

    coefs_ens!1 = "predictors"
    coefs_ens&predictors = (/1,2/)
    printVarSummary(coefs_ens)






   
    outDir = "/glade/work/kumar34/hydroclimate_ncc/"

    outfile = outDir+"NorthAmerica_obs_enso_h2osoi_regression_memoryENSO_5yrout.nc"
   
     system("/bin/rm -f "+outfile)
     ncdf = addfile(outfile,"c")
     ncdf->coefs = coefs_ens
     ncdf->actskill = actskill
     ncdf->sig2noise = sig2noise
     ncdf->predskill = predskill
     ncdf->varmodM = var1ModM
     ncdf->varmodE = var1ModE
     ncdf->var     = var1
     ncdf->sm_std = sm_std
     


 delete(ncdf)







   






