;*****pec_1.ncl
;
; Concepts illustrated:
;   - Calculating and plotting spectra
;************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
;************************************************

;************************************************
; variable and file handling
;************************************************

 dirIn = "/glade/work/kumar34/hydroclimate_ncc/"

 nyrp = 30
 yrstr = (/1940,1970,2000,2030,2060/)
 yrend = yrstr + nyrp - 1

 region = (/"Southwest","Southeast","SouthernUSPlains","NorthernUSPlains","CanadianPlains_rev1"/)
 nrgs = dimsizes(region)

 model = (/"H2OSOI_CESM-LE", "SOILM_GFDLCM3-LE"/) 

 regN = 0
 modN = 1
 
 f1 = addfile(dirIn+model(modN)+"_1940-2089_30yr_TimeSeries_"+region(regN)+".nc", "r")
 print(f1)

 var_mod = f1->gfdlcm3_soilm(:, :, :)    ; cesm_soilm(periods, ens, time) 
 
 var_syn = f1->gfdlcm3_synt_soilm(:, :, :, :)   ; cesm_synt_soilm(periods, iter, ens, time)  

 f2 = addfile(dirIn+"ERA5_obs_40ens_30yr_TimeSeries_"+region(regN)+".nc", "r")
 obs_mod = f2->soilm(0:19, :)     ;(nEns, time)

 obs_syn = f2->synt_soilm(:, 0:19, :)  ; (iter, ens, time)
 
 
 ;******************************
 
  var_mod(0, :, :) = (/obs_mod/)
  var_syn(0, :, :, :) = (/obs_syn/)
 
;************************

   d   = 0
   sm  = 1         ; periodogram
   pct = 0.10

  ;************************************************
  ; calculate mean spectrum spectrum and lag1 auto cor
  ;************************************************
  
  ; loop over each segment of length ntim
  
  
  
 ;  spcavg = new(5, floattoint(ntim/2), typeof(var_mod)) 
  
   ntim = dimsizes(var_mod&time) 

   nseg = dimsizes(var_mod&ens)
  
   dlen = ntim   

   spcavg_itr = new((/1000, 5, floattoint(ntim/2)/), typeof(var_mod), var_mod@_FillValue)
  
   spctemp_itr = new((/5, floattoint(ntim/2)/), typeof(var_mod), var_mod@_FillValue)


   spcavg_itr_le = new((/1000, 5, floattoint(ntim/2)/), typeof(var_mod), var_mod@_FillValue)

   ;spctemp_itr_le = new((/5, floattoint(ntim/2)/), typeof(var_mod), var_mod@_FillValue)



   spcavg_itr = 0.0
   spctemp_itr = 0.0


   printVarSummary(var_syn)
  

  var_mod_hist_var = var_mod(:, 0, 0)

  do i = 0, 4
 
     var_mod_hist = var_mod(i, :, :)
     var_mod_hist_1D = ndtooned(var_mod_hist)
     var_mod_hist_var(i) = dim_variance_n_Wrap(var_mod_hist_1D, 0)

  end do
 
  print(var_mod_hist_var)

  
  do i = 0, 0
  do itr = 0, 999
     x = var_syn(i, itr, :, :)
     spctemp_itr = 0.0
     do n = 0, nseg-1
        dof := specx_anal(x(n,:),d,sm,pct)
        temp1 := dof@spcx
        ;temp2 := temp1*100/(dim_variance_n(x(n, :), 0)*dlen)
        spctemp_itr(i, :) = spctemp_itr(i, :) + temp1
    end do
       temp11 := ndtooned(x)
       temp11_var := dim_variance_n(temp11, 0)
       ;spctemp_itr(i, :) = spctemp_itr(i, :) * 100 / (temp11_var*dlen*nseg)
       spctemp_itr(i, :) = spctemp_itr(i, :) * 100 / (var_mod_hist_var(0)*dlen*nseg)       

       spcavg_itr(itr, i, :) = spctemp_itr(i, :)

  end do
  end do

  
 
  random_setallseed(36484749, 9494848)

  do i = 0, 0
  do itr = 0, 999
     ;random_setallseed(36484749, 9494848)  
     unf = toint(random_uniform(0, nseg, nseg))
     x = var_mod(i, unf, :)
     spctemp_itr = 0.0
     do n = 0, nseg-1
        dof := specx_anal(x(n,:),d,sm,pct)
        temp1 := dof@spcx
        ;temp2 := temp1*100/(dim_variance_n(x(n, :), 0)*dlen)
        spctemp_itr(i, :) = spctemp_itr(i, :) + temp1
    end do
       temp11 := ndtooned(x)
       temp11_var := dim_variance_n(temp11, 0)
       spctemp_itr(i, :) = spctemp_itr(i, :) * 100 / (var_mod_hist_var(0)*dlen*nseg)
       ;spctemp_itr(i, :) = spctemp_itr(i, :) * 100 / (temp11_var*dlen*nseg)     
 
       ;spctemp_itr(i, :) = spctemp_itr(i, :) / nseg
       ;spcavg_itr(itr, i, :) = dim_avg_n(spctemp_itr(i, :), 0)

       spcavg_itr_le(itr, i, :) = spctemp_itr(i, :)

  end do
  end do


  do i = 1, 4
  do itr = 0, 999
     x = var_syn(i, itr, :, :)
     spctemp_itr = 0.0
     do n = 0, nseg-1
        dof := specx_anal(x(n,:),d,sm,pct)
        temp1 := dof@spcx
        ;temp2 := temp1*100/(dim_variance_n(x(n, :), 0)*dlen)
        spctemp_itr(i, :) = spctemp_itr(i, :) + temp1
    end do
       temp11 := ndtooned(x)
       temp11_var := dim_variance_n(temp11, 0)
       ;spctemp_itr(i, :) = spctemp_itr(i, :) * 100 / (temp11_var*dlen*nseg)
        spctemp_itr(i, :) = spctemp_itr(i, :) * 100 / (var_mod_hist_var(1)*dlen*nseg)       


       spcavg_itr(itr, i, :) = spctemp_itr(i, :)

  end do
  end do

  
 
  random_setallseed(36484749, 9494848)

  do i = 1, 4
  do itr = 0, 999
     ;random_setallseed(36484749, 9494848)  
     unf = toint(random_uniform(0, nseg, nseg))
     x = var_mod(i, unf, :)
     spctemp_itr = 0.0
     do n = 0, nseg-1
        dof := specx_anal(x(n,:),d,sm,pct)
        temp1 := dof@spcx
        ;temp2 := temp1*100/(dim_variance_n(x(n, :), 0)*dlen)
        spctemp_itr(i, :) = spctemp_itr(i, :) + temp1
    end do
       temp11 := ndtooned(x)
       temp11_var := dim_variance_n(temp11, 0)
       spctemp_itr(i, :) = spctemp_itr(i, :) * 100 / (var_mod_hist_var(1)*dlen*nseg)
       ;spctemp_itr(i, :) = spctemp_itr(i, :) * 100 / (temp11_var*dlen*nseg)     
 
       ;spctemp_itr(i, :) = spctemp_itr(i, :) / nseg
       ;spcavg_itr(itr, i, :) = dim_avg_n(spctemp_itr(i, :), 0)

       spcavg_itr_le(itr, i, :) = spctemp_itr(i, :)

  end do
  end do










  spcavg_le = dim_median_n(spcavg_itr_le, 0)
  
  spcavg = dim_median_n(spcavg_itr, 0)

  spcavg_low = spcavg

  spcavg_high = spcavg


   print(sum(spcavg)/5)


   ip1 = dim_pqsort_n(spcavg_itr, 2, 0)

   spcavg_low = spcavg_itr(24, :, :)
   spcavg_high = spcavg_itr(974, :, :) 


  ;************************************************
  ; plotting
  ;************************************************
   
      typPlt = "x11"
       ;typPlt@wkWidth  = 3000
       ;typPlt@wkHeight = 3000
 

    wks  = gsn_open_wks(typPlt, region(regN)+"_"+model(modN)+"syn_unc2_rev5")              ; Opens a ps file
  
    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@gsnMaximize = True

    res@vpHeightF                = 0.4
    res@vpWidthF                 = 0.6



    res@tiXAxisString = "Frequency (cycles/"+dlen+"year)"  ; xaxis
    res@tiYAxisString = "% variance explained"                  ; yaxis
  


    ;res@vpHeightF= 0.2                    ; change aspect ratio of plot
    ;res@vpWidthF = 0.35                  
  
    res@trYMinF  = 2.0                   ; min value on y-axis
    res@trYMaxF  = 16.0 +19.0            ; max value on y-axis

    res@trXMinF  = 1.0                   ; min value on y-axiss
    res@trXMaxF  = 9.0                   ; max value on y-axis

     cmap = read_colormap_file("MPL_rainbow")
     idc = (/2, 20, 40, 93, 127/)
     clevs = cmap(idc,:)

     print(clevs)

     ;exit
     
      lcs = (/"black", "deepskyblue3","gold3","green", "gold3"/)


     idc2 = idc-2
     clevs2 = cmap(idc2, :)


     res@xyDashPattern  = 0                         ; Make curves all solid

     res@xyLineThicknesses = (/  5.0,   5.0, 5.0, 5.0, 5.0/) + 10.0

     res@xyLineColors      = lcs

     printVarSummary(clevs)

     ;spcavg(0, :) = spcavg@_FillValue
     spcavg(2:3, :) = spcavg@_FillValue


     plot = gsn_csm_xy(wks, dof@frq*dlen, spcavg, res)      ; create plot


     res@xyDashPattern  = 2
     
     ;spcavg_le(0, :) = spcavg@_FillValue
     spcavg_le(2:3, :) = spcavg@_FillValue

     clevsM = clevs
    ;clevsM(0) = "deepskyblue"
    ;clevsM(4) = "firebrick3"


     resM = res
     delete(resM@xyDashPattern)
     resM@xyMarkLineMode          = "MarkLines"
     resM@xyLineThicknesses = (/  5.0,   5.0, 5.0, 5.0, 5.0/) - 2.0
     resM@xyMarker                = 12
     resM@xyMarkerSizeF            = 0.020
     resM@xyMarkerThicknessF       = 5
     ;resM@xyMarkerColors = (/"blue", "deepskyblue4","darkgoldenrod2","green", "darkgoldenrod2"/)   
       resM@xyMarkerColors = (/"black", "deepskyblue4","darkgoldenrod2","green", "darkgoldenrod4"/)
     resM@xyLineColors = (/"black", "deepskyblue4","darkgoldenrod2","green", "darkgoldenrod4"/)
     


    
     plot2 = gsn_csm_xy(wks, dof@frq*dlen, spcavg_le, resM)

     overlay(plot, plot2)

     bot_plot = new(5, graphic)

   res@gsnXYFillOpacities = 0.5 


    do i = 0, 0

      minmax_fill = spcavg_low(0:1, :)
      minmax_fill(0, :) = spcavg_low(i, :)
      minmax_fill(1, :) = spcavg_high(i, :)


      delete(res@xyLineColors)
      ;res@gsnXYFillColors = cmap(idc(i), :)
      res@gsnXYFillColors = "gray60"
      res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
      bot_plot(i)  = gsn_csm_xy (wks, dof@frq*dlen, minmax_fill, res)  ; Create filled XY plot.
      overlay(plot, bot_plot(i))

   end do





    do i = 1, 1

      minmax_fill = spcavg_low(0:1, :)
      minmax_fill(0, :) = spcavg_low(i, :) 
      minmax_fill(1, :) = spcavg_high(i, :) 

 
      delete(res@xyLineColors)
      ;res@gsnXYFillColors = cmap(idc(i), :)
      res@gsnXYFillColors = "deepskyblue"
      res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
      bot_plot(i)  = gsn_csm_xy (wks, dof@frq*dlen, minmax_fill, res)  ; Create filled XY plot.
      overlay(plot, bot_plot(i))

   end do

     do i = 4, 4

      minmax_fill = spcavg_low(0:1, :)
       minmax_fill(0, :) = spcavg_low(i, :)
      minmax_fill(1, :) = spcavg_high(i, :)

      ;res@gsnXYFillColors = cmap(idc(i), :)
      res@gsnXYFillColors = "gold"
      res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
      bot_plot(i)  = gsn_csm_xy (wks, dof@frq*dlen, minmax_fill, res)  ; Create filled XY plot.
      overlay(plot, bot_plot(i))

   end do




    draw(plot)

    res_text               = True
    res_text@txFontHeightF = 0.025
    res_text@txJust        = "CenterLeft"

    res_lines                   = True        ; polyline mods desired

    yy_ar = (/4.0, 5.0, 6.0, 7.0, 7.0/) + 15.0

    ;yy_ar = (/4.0, 5.0, 6.0, 7.0, 8.0/) * 3.0

    txt_arr = (/"1940-69", "1970-99", "2000-29", "2030-59", "2060-89"/)

   ;txt_arr = (/"1970-99", "2060-89"/)


    do i = 1, 1

    res_lines@gsLineColor       = lcs(i) 
    res_lines@gsLineThicknessF  = 5 + 3.0           ; line thicker
    res_lines@gsLineDashPattern = 0   ; dash pattern

    xx = (/1.4,2.0/) 
    xx = xx
    yy = (/yy_ar(i), yy_ar(i)/)
    ;gsn_polyline(wks,plot,xx,yy,res_lines)
    ;gsn_text(wks,plot,txt_arr(i),2.2,yy_ar(i),res_text)

   end do


   do i = 4, 4

    res_lines@gsLineColor       = lcs(i) 
    res_lines@gsLineThicknessF  = 5 + 3.0           ; line thicker
    res_lines@gsLineDashPattern = 0   ; dash pattern

    xx = (/1.4,2.0/) 
    xx = xx
    yy = (/yy_ar(i), yy_ar(i)/)
    ;gsn_polyline(wks,plot,xx,yy,res_lines)
    ;gsn_text(wks,plot,txt_arr(i),2.2 ,yy_ar(i),res_text)

   end do




    ;yy = (/52.5,52.5/)
    ;res_lines@gsLineColor       = line_col2   ; line color
    ;res_lines@gsLineDashPattern = line_pat2   ; dash pattern
    ;gsn_polyline(wks,plot,xx,yy,res_lines)
    ;gsn_text(wks,plot,m2,1900,52.5,res_text)  

   frame(wks) 









