; Here I test how many years and/or ensemble members are needed to achieve a stable climatology,
; as measured by a high signal-to-total ratio (see Dirmeyer et al. 2013, eqs 1-4)

; For faster computations, analysis is restricted to North America
 latmin = 10
 latmax = 75
 lonmin = 360-170
 lonmax = 360-50


;********************************************************************************

 inDir = "/glade/work/kumar34/GFDL-CM3LE/PET"


 flm = addfile(inDir+"/PET_pt_Amon_GFDL-CM3_historical_rcp85_r10i1p1_192001-210012.nc", "r")

 lat = flm->lat({latmin:latmax})
 lon = flm->lon({lonmin:lonmax})
 time = flm->time

 ntim = dimsizes(time)
 
  pet_pt = flm->PET

 
 print(ntim)

 flm1 = addfile("/glade/collections/cdg/data/CLIVAR_LE/gfdl_cm3_lens/Lmon/mrsos/mrsos_Lmon_GFDL-CM3_historical_rcp85_r10i1p1_192001-210012.nc", "r")
 lm = flm1->mrsos
 
  mrsos = lm(0, {latmin:latmax},{lonmin:lonmax})

 

  
;exit

;****************************************************************************

; create variable with all 40 members' historical+rcp85 monthly time series 
 
 nEns = 20
 ens_ind = new((/nEns/), "integer", -99)
 ens_ind(0:19) = ispan(1, 20, 1)
 
 EMtsr = new((/nEns,ntim,dimsizes(lat),dimsizes(lon)/),"float",pet_pt@_FillValue)

 inDirPET = "/glade/work/kumar34/GFDL-CM3LE/PET/"


 
 do i = 0, nEns-1
     
      ft1 = addfile(inDir+"/PET_pt_Amon_GFDL-CM3_historical_rcp85_r"+ens_ind(i)+"i1p1_192001-210012.nc", "r")
     
     
      var := ft1->PET(:, {latmin:latmax},{lonmin:lonmax})
      var = mask(var, .not.ismissing(mrsos), True)        
      EMtsr(i, :, :,:) = var
     

     delete(var)

     
 end do

 ;exit

 ;date = array_append_record(date1(sti1:endi1),date2(sti2:endi2),0)
 
 time1 = ft1->time
 
 EMtsr&time = time
 EMtsr!0 = "ens"
 EMtsr&ens = ens_ind

 ;delete([/var,sti1,sti2,endi1,endi2,time1,time2,date1,date2,ft1,ft2/])

 printVarSummary(EMtsr)
 
 EMtsr = runave_n_Wrap(EMtsr, 12, 0, 1)
 
 printVarSummary(EMtsr)

 
 sm_smon = EMtsr

 delete(EMtsr)


 sm_smon_em = dim_avg_n_Wrap(sm_smon, 0)

 printVarSummary(sm_smon_em)

 sm_smean_anom = sm_smon

 
 sm_smean_em_temp =  sm_smon_em(0:179, :, :)

   
  sm_memory = sm_smon(:, 0:4, :, :)
  sm_memory = sm_smon@_FillValue
 
  styr = (/1940, 1970, 2000, 2030, 2060/)

  nlag = 24


  do i = 0, 3

      st1 = (styr(i) - 1920)*12 
      end1 = st1 + 360 -1

      temp11 := sm_smon(:, st1:end1, :, :)

      temp11_std := dim_avg_n_Wrap(temp11, 1)

       sm_memory(:, i, :, :) = (/temp11_std/)

           
       delete(temp11)
       delete(temp11_std)
    

  end do
  
  
  
  do i = 4, 4

      st1 = (styr(i) - 1920)*12 
      end1 = st1 + 360 -1

      temp11 := sm_smon(:, st1:end1, :, :)

      temp11_std := dim_avg_n_Wrap(temp11, 1)

       sm_memory(:, i, :, :) = (/temp11_std/)
    
       delete(temp11)
       delete(temp11_std)
    

  end do

  
   dirout = "/glade/work/kumar34/hydroclimate_ncc/"


   outfile = dirout+"PET_GFDL-CM3-LE_five30yr_period_sm_mean.nc"
   system("/bin/rm -f "+outfile)
   ncdf1 = addfile(outfile,"c")
   ncdf1->pet_mean = sm_memory 
   
 
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
