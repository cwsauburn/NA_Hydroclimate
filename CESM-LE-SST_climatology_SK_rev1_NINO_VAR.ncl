; Here I test how many years and/or ensemble members are needed to achieve a stable climatology,
; as measured by a high signal-to-total ratio (see Dirmeyer et al. 2013, eqs 1-4)

; For faster computations, analysis is restricted to North America
  latmin = 10
  latmax = 75
  lonmin = 360-170
  lonmax = 360-50
 
;***********************


 url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"

 dirgfdl = url+"CESM-LE/ocn/postprocessed/"


 ;dirgfdl = "/Projects/SERDP_regimeshifts/data/GFDL-CM3-LE/SOILM/"
 ;filegfdl = dirgfdl+"SOILM_GFDL-CM3-LE_20em_192001-210012_12moRunningMeans.nc"
 ;filegfdl = dirgfdl+"SOILM_GFDL-CM3-LE_20em_192001-210012_12mormAnomsEnsMean.nc"
 filegfdl = dirgfdl+"IndoPacific_SST_CESM-LE_40members_1920-2100_12moRunMeans.nc" 
 flg = addfile(filegfdl,"r")
 print(flg)


 ;exit


 var1_t = flg->SST
 
 sm_smon = var1_t(ens|:, time|:, nlat|:, nlon|:)
 
 delete(var1_t)
 
 
 sm_smon_em = dim_avg_n_Wrap(sm_smon, 0)

 printVarSummary(sm_smon_em)

 sm_smean_anom = sm_smon

 
 sm_smean_em_temp =  sm_smon_em(0:179, :, :)


 do i = 0, 11

    do j = 0, 179
       sm_smean_em_temp(j, :, :) =  sm_smon_em(j*12+i, :, :)
    end do
    
    temp12 = runave_n_Wrap(sm_smean_em_temp, 21, 0, 0)

    do k = 0, 39
    do j = 0, 179
       sm_smean_anom(k, j*12+i, :, :) = (/sm_smon(k, j*12+i, :, :)-temp12(j, :, :)/)  
    end do
    end do
  
 end do 
  
  sm_memory = sm_smon(:, 0:4, :, :)
 
  styr = (/1940, 1970, 2000, 2030, 2060/)

  nlag = 24


;  T1 = fspan(1, 24, 96)
;  T2 = fspan(1, 2, 12)


;  data_er = conform_dims((/dimsizes(T1), nlag+1, dimsizes(sm_memory&lat), dimsizes(sm_memory&lon)/), sm_memory(0:nlag, 0, :, :), (/1, 2, 3/))  

 ; data_er = data_er@_FillValue

;  do i = 0, dimsizes(T1)-1
;     x1 = ispan(0, 24, 1)
;     y1 = exp(-x1/T1(i))
;     do j = 0, dimsizes(y1)-1
;        data_er(i, j, :, :) = (/y1(j)/)
;     end do
;     delete(y1)     
; end do


;  data_er2 = data_er


  do i = 0, 4

      st1 = (styr(i) - 1920)*12 
      end1 = st1 + 360 -1

      temp11 = sm_smean_anom(:, st1:end1, :, :)

      temp11_std = dim_stddev_n_Wrap(temp11, 1)

      temp11_std = (/temp11_std*1000.0/)



      ;temp121 = esacr_n(temp11, 24, 1)

     ; mem_temp = temp121(0, 0, :, :)
      ;mem_temp = mem_temp@_FillValue

     ;do i_ens = 0, 39
      ;do i_ens = 0, 0

      ;  data_er2 = data_er@_FillValue
      ;  mem_temp = mem_temp@_FillValue

      ;  do i_div = 0, dimsizes(T1)-1

      ;     data_er2(i_div, :, :, :) = (/(data_er(i_div, :, :, :) - temp121(i_ens, :, :, :))^2/)
     
      ;   end do
    
     ;   data_er2_sum = dim_sum_n_Wrap(data_er2, 1)            

     ;   data_er2_sum_min = dim_min_n_Wrap(data_er2_sum, 0)

       
      ;   do i_div = 0, dimsizes(T1)-1          
      ;       mem_temp = where(data_er2_sum(i_div, :, :).eq.data_er2_sum_min, T1(i_div), mem_temp)
      ;   end do    
      

       sm_memory(:, i, :, :) = (/temp11_std/)

       

      ;end do


       delete(temp11)
       delete(temp11_std)
       ;delete(temp121)


  end do

 ;exit 

   dirout = "/glade/work/kumar34/hydroclimate_ncc/"


   outfile = dirout+"SST_CESM-LE_five30yr_period_sm_var.nc"
   system("/bin/rm -f "+outfile)
   ncdf1 = addfile(outfile,"c")
   ncdf1->sst_var = sm_memory 
   
 
 exit

 
 

 







  
;****************************************************************************

; compute sample climatologies within each ensemble member varying the number of years included (Xen in Eq 1)
 
 nYrClim = ispan(5,40,5)       ; number of years used to compute climatology
 nClims = dimsizes(nYrClim)    ; number of climatology cases

 maxNSample = (2080-1920+1)/5  ; number of climatology periods (samples) available in time series (1920-2080)

 print_table([/"computing individual ensemble member sample climatologies..."/],"%s")
 
 EMsamplClim = new((/nEns,nClims,maxNSample,12,dimsizes(lat),dimsizes(lon)/),"float",area@_FillValue)
 EMsamplClim!0 = "nEns"
 EMsamplClim&nEns = ens_ind
 EMsamplClim!1 = "nYrClim"
 EMsamplClim&nYrClim = nYrClim
 EMsamplClim!2 = "nSamples"
 EMsamplClim&nSamples = ispan(1,32,1)

 do i = 0, nEns-1              ; loop through ensemble members
   EnsMem = EMtsr(i,:,:,:)
   
   do c = 0,nClims-1           ; loop through cases of different climatology periods
     ny = nYrClim(c)
     stdate := ispan(19200201,20750201,ny*10000)
     endate := stdate + ny*10000 - 100
     nSmpl = dimsizes(stdate)
     if(endate(nSmpl-1).gt.20800101) then
       stdate := stdate(0:nSmpl-2)
       endate := endate(0:nSmpl-2)
       nSmpl = dimsizes(stdate)
     end if

     do s = 0,nSmpl-1          ; loop through each sample and compute the climatology 
       sti = ind(date.eq.stdate(s))
       endi = ind(date.eq.endate(s))
       SmpClim = clmMonTLL(EnsMem(sti:endi,:,:))
       EMsamplClim(i,c,s,:,:,:) = SmpClim
     end do

   end do
 end do

 delete([/EMtsr,sti,endi,stdate,endate,EnsMem,SmpClim/])

;****************************************************************************

; compute ensemble mean sample climatologies varying the number of ensemble members used (Xnbar in Eq 1)
 
 nEnsUsed = ispan(5,40,5)      ; number of ensemble members used for climatology
 nEnsU = dimsizes(nEnsUsed)    ; number of ensemble cases

 print_table([/"computing ensemble mean sample climatologies and variances..."/],"%s")

 EnsMeanSmpClim = new((/nEnsU,nClims,maxNSample,12,dimsizes(lat),dimsizes(lon)/),"float",area@_FillValue)
 EnsMeanSmpClim!0 = "nEnsUsed"
 EnsMeanSmpClim&nEnsUsed = nEnsUsed
 EnsMeanSmpClim!1 = "nYrClim"
 EnsMeanSmpClim&nYrClim = nYrClim
 EnsMeanSmpClim!2 = "nSamples"
 EnsMeanSmpClim&nSamples = ispan(1,32,1)

 EnsMeanVar = EnsMeanSmpClim(:,:,0,:,:,:)
 EnsMembVar = EnsMeanVar
 
 do e = 0, nEnsU-1             ; loop through cases of different number of ensemble members considered
   ne = nEnsUsed(e)
   EnsUsmpClim := EMsamplClim(0:ne-1,:,:,:,:,:)

   do c = 0,nClims-1           ; loop through cases of different climatology periods
     ny = nYrClim(c)
     nSmpl = (2080-1920+1)/ny
     EnsUclims := EnsUsmpClim(:,c,0:nSmpl-1,:,:,:)

     do s = 0,nSmpl-1          ; loop through each climatology sample and compute its ensemble mean
       EnsUclimSmpl := EnsUclims(:,s,:,:,:)
       EnsUSmplMean = dim_avg_n_Wrap(EnsUclimSmpl,0)
       EnsMeanSmpClim(e,c,s,:,:,:) = EnsUSmplMean
     end do

     EnsMeanVar(e,c,:,:,:) = dim_variance_n_Wrap(EnsMeanSmpClim(e,c,:,:,:,:),0) ; (Vs in Eq 3)
     EnsMembVar(e,c,:,:,:) = dim_variance_n_Wrap(EnsUclims,(/0,1/))             ; (Vn in Eq 4)
   end do
 end do

;****************************************************************************

; compute Signal-to-Total (S2T = Vs/(Vs+Vn))

 S2T = EnsMeanVar / (EnsMeanVar + EnsMembVar)
 copy_VarCoords(EnsMeanVar,S2T)
 
 delete([/EnsUsmpClim,EnsUclims,EnsUclimSmpl,EnsUSmplMean/])
 
;****************************************************************************

; select coordinates for a grid point and plot one graphic for each month

 print_table([/"preparing plot..."/],"%s")

 la = 40.
 lo = 360. - 105.

 ptS2T = S2T(:,:,:,{la},{lo})

 ; Plot Resources

 res = True
 res@gsnDraw              = False
 res@gsnFrame             = False
 res@gsnCenterStringFontHeightF = 0.04
 
 res@cnFillOn             = True
 res@cnFillMode           = "RasterFill"
 res@cnLinesOn            = False
 res@cnLevelSelectionMode = "ManualLevels"
 res@cnMinLevelValF       = 0.
 res@cnMaxLevelValF       = 0.6
 res@cnLevelSpacingF      = 0.02

 res@lbOrientation        = "Vertical"

 res@tiXAxisString        = "# of years in climatology"
 res@tiXAxisFontHeightF   = 0.035
 res@tiYAxisString        = "# of ensemble members"
 res@tiYAxisFontHeightF   = 0.035
 res@tmXBMinorOn          = False
 res@tmXBLabelFontHeightF = 0.03
 res@tmYLMinorOn          = False
 res@tmYLLabelFontHeightF = 0.03

 pres = True
 pres@gsnPanelLabelBar    = False
 pres@gsnPanelMainString  = "Signal-to-Total Ratio for Monthly Soil Moisture Climatologies"+ \
                            "~C~       between 1920 and 2080 at "+la+"N, "+toint(360.-lo)+"W"
 pres@gsnPanelMainFontHeightF = 0.015
 
; outDir = "/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/plots/"
 outDir = "/glade/u/home/cdd4f/SERDP2017/output/CESM-LE/"
 filplot = outDir+"SoilMoisture_Signal2Total_"+la+"N-"+toint(360.-lo)+"W"
 type = "png"
 
 wks = gsn_open_wks(type,filplot)
 gsn_define_colormap(wks,"MPL_gist_earth")
 gsn_reverse_colormap(wks)

 plot = new(12,graphic)
 months = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

 do m = 0,11
   res@gsnCenterString = months(m)
   plot(m) = gsn_csm_contour(wks,ptS2T(:,:,m),res)
 end do
 
 gsn_panel(wks,plot,(/4,3/),pres)
 system("convert -trim "+filplot+"."+type+" "+filplot+"."+type)

 
;****************************************************************************
;****************************************************************************
