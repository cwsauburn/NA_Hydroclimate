; Here I test how many years and/or ensemble members are needed to achieve a stable climatology,
; as measured by a high signal-to-total ratio (see Dirmeyer et al. 2013, eqs 1-4)

; For faster computations, analysis is restricted to North America
 latmin = 10
 latmax = 75
 lonmin = 360-170
 lonmax = 360-50

;****************************************************************************

; read in surface data file for grid variables
 
 fm = addfile("/glade/work/kumar34/prj_data/surf_data/CLM4.0CN.h0.fixed.nc", "r")
; fm = addfile("/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/datafiles/CLM4.0CN.h0.fixed.nc", "r")
 lm = fm->landmask({latmin:latmax},{lonmin:lonmax})
 area = fm->area({latmin:latmax},{lonmin:lonmax})
 landfrac = fm->landfrac({latmin:latmax},{lonmin:lonmax})
 area@_FillValue = 1.0*10^36
 landfrac@_FillValue = 1*0*10^36
 arwt = area
 arwt = area*landfrac/100
 lat = fm->lat({latmin:latmax})
 lon = fm->lon({lonmin:lonmax})

;****************************************************************************

; get time series dimensions
 
 inDir = "/glade/campaign/cesm/collections/cesmLE/CESM-CAM5-BGC-LE/lnd/proc/tseries/monthly/H2OSOI/"
; inDir = "/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/datafiles/H2OSOI/"
 run = (/"B20TRC5CNBDRD","BRCP85C5CNBDRD"/)

 file1 = systemfunc("ls "+inDir+"b.e11."+run(0)+".f09_g16.002.clm2.h0.H2OSOI.192001-200512.nc")
 ft1 = addfile(file1,"r")
 date1 = ft1->mcdate
 nt1 = dimsizes(date1)

 file2 = systemfunc("ls "+inDir+"b.e11."+run(1)+".f09_g16.105.clm2.h0.H2OSOI.200601-210012.nc")
 ft2 = addfile(file2,"r")
 date2 = ft2->mcdate
 nt2 = dimsizes(date2)
 nt2 = nt2 - 12

 ntim = nt1 + nt2

 print(ntim)

 
 
;exit

;****************************************************************************

; create variable with all 40 members' historical+rcp85 monthly time series 
 
 nEns = 40
 ens_ind = new((/nEns/), "integer", -99)
 ens_ind(0:34) = ispan(1, 35, 1)
 ens_ind(35:39) = ispan(101, 105, 1)

 ly_wt = (/0.0175, 0.0276, 0.0455, 0.0750, 0.1236, 0.2038, 0.3360/)

 print_table([/"compiling ensemble members..."/],"%s")

 EMtsr = new((/nEns,ntim,dimsizes(lat),dimsizes(lon)/),"float",area@_FillValue)
 
 do i = 0, nEns-1
 ;do i = 0, 0
   file1 = systemfunc("ls "+inDir+"b.e11."+run(0)+".f09_g16."+sprinti("%0.3i", ens_ind(i))+".clm2.h0.H2OSOI.*.nc")
   ft1 := addfile(file1,"r")
   date1 := ft1->mcdate
   sti1 = ind(date1.eq.19200201)
   endi1 = ind(date1.eq.20060101)
   var_t := ft1->H2OSOI(sti1:endi1,0:6,{latmin:latmax},{lonmin:lonmax})
   var := dim_avg_wgt_n_Wrap(var_t, ly_wt, 1, 1) 
   delete(var_t)
   EMtsr(i,0:nt1-1,:,:) = var
   ;print(sti1)
   ;print(endi1)

   file2 := systemfunc("ls "+inDir+"b.e11."+run(1)+".f09_g16."+sprinti("%0.3i", ens_ind(i))+".clm2.h0.H2OSOI.*.nc")
   ft2 := addfiles(file2,"r")
   date2 := ft2[:]->mcdate
   sti2 = ind(date2.eq.20060201)
   endi2 = ind(date2.eq.21000101)
   print(sti2)
   print(endi2)
   ;
   var_t := ft2[:]->H2OSOI(sti2:endi2,0:6,{latmin:latmax},{lonmin:lonmax})
   var := dim_avg_wgt_n_Wrap(var_t, ly_wt, 1, 1) 
   delete(var_t)
   EMtsr(i,nt1:ntim-1,:,:) = var
 end do

 date = array_append_record(date1(sti1:endi1),date2(sti2:endi2),0)
 
 time1 = ft1->time(sti1:endi1)
 time2 = ft2[:]->time(sti2:endi2)
 time2 = time_to_newtime(time2,time1@units)
 time = array_append_record(time1,time2,0)

 EMtsr&time = time
 EMtsr!0 = "ens"
 EMtsr&ens = ens_ind

 delete([/var,sti1,sti2,endi1,endi2,time1,time2,date1,date2,ft1,ft2/])

 printVarSummary(EMtsr)
 
 EMtsr = runave_n_Wrap(EMtsr, 12, 0, 1)
 
 printVarSummary(EMtsr)

 
 sm_smon = EMtsr

 delete(EMtsr)


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


  do i = 0, 4

      st1 = (styr(i) - 1920)*12 
      end1 = st1 + 360 -1

      temp11 = sm_smean_anom(:, st1:end1, :, :)

      temp121 = esacr_n(temp11, 12, 1)

      sm_memory(:, i, :, :) = (/temp121(:, 12, :, :)/)

      delete(temp11)
      delete(temp121)


  end do

 ;exit 

   dirout = "/glade/work/kumar34/hydroclimate_ncc/"


   outfile = dirout+"H2OSOI_CESM-LE_five30yr_period_sm_memory_lag12.nc"
   system("/bin/rm -f "+outfile)
   ncdf1 = addfile(outfile,"c")
   ncdf1->sm_memory = sm_memory 
   
 
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
