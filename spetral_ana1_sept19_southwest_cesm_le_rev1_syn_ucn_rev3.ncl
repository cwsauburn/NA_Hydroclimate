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

 region = (/"Southwest","Southeast","SouthernUSPlains","NorthernUSPlains","CanadianPlains"/)
 nrgs = dimsizes(region)

 model = (/"H2OSOI_CESM-LE", "SOILM_GFDLCM3-LE"/) 

 regN = 0
 modN = 0
 
 f1 = addfile(dirIn+model(modN)+"_1940-2089_30yr_TimeSeries_"+region(regN)+".nc", "r")
 print(f1)

 var_mod = f1->cesm_soilm(:, :, :)    ; cesm_soilm(periods, ens, time) 
 
 var_syn = f1->cesm_synt_soilm(:, :, :, :)   ; cesm_synt_soilm(periods, iter, ens, time)  
 
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

 


  do i = 0, 4
  do itr = 0, 999
     x = var_syn(i, itr, :, :)
     spctemp_itr = 0.0
     do n = 0, nseg-1
        dof := specx_anal(x(n,:),d,sm,pct)
        temp1 := dof@spcx
        temp2 := temp1*100/(dim_variance_n(x(n, :), 0)*dlen)
        spctemp_itr(i, :) = spctemp_itr(i, :) + temp2
    end do
       spctemp_itr(i, :) = spctemp_itr(i, :) / nseg
       ;spcavg_itr(itr, i, :) = dim_avg_n(spctemp_itr(i, :), 0)

       spcavg_itr(itr, i, :) = spctemp_itr(i, :)

  end do
  end do


  do i = 0, 4
  do itr = 0, 999
     unf = toint(random_uniform(0, nseg, nseg))
     x = var_mod(i, unf, :)
     spctemp_itr = 0.0
     do n = 0, nseg-1
        dof := specx_anal(x(n,:),d,sm,pct)
        temp1 := dof@spcx
        temp2 := temp1*100/(dim_variance_n(x(n, :), 0)*dlen)
        spctemp_itr(i, :) = spctemp_itr(i, :) + temp2
    end do
       spctemp_itr(i, :) = spctemp_itr(i, :) / nseg
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

   spcavg_low = spcavg_itr(25, :, :)
   spcavg_high = spcavg_itr(975, :, :) 

  ;************************************************
  ; plotting
  ;************************************************
    wks  = gsn_open_wks("x11", region(regN)+"_"+model(modN)+"syn_unc2")              ; Opens a ps file
  
    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@gsnMaximize = True


     res@tiXAxisString = "Frequency (cycles/"+dlen+"year)"  ; xaxis
    res@tiYAxisString = "% variance explained"                  ; yaxis
  


    res@vpHeightF= 0.2                    ; change aspect ratio of plot
    res@vpWidthF = 0.35                  
  
    res@trYMinF  = 2.0                   ; min value on y-axis
    res@trYMaxF  = 16.0 + 4.0           ; max value on y-axis

    res@trXMinF  = 1.0                   ; min value on y-axiss
    res@trXMaxF  = 10.0                   ; max value on y-axis

     cmap = read_colormap_file("MPL_rainbow")
     idc = (/2, 20, 40, 93, 127/)
     clevs = cmap(idc,:)

     idc2 = idc-2
     clevs2 = cmap(idc2, :)


     res@xyDashPattern  = 0                         ; Make curves all solid

     res@xyLineThicknesses = (/  5.0,   5.0, 5.0, 5.0, 5.0/) + 3.0

       res@xyLineColors      = clevs

    printVarSummary(clevs)


     plot = gsn_csm_xy(wks, dof@frq*dlen, spcavg, res)      ; create plot


     res@xyDashPattern  = 1
     
     plot2 = gsn_csm_xy(wks, dof@frq*dlen, spcavg_le, res)

     overlay(plot, plot2)

     bot_plot = new(5, graphic)

   res@gsnXYFillOpacities = 0.2 


    do i = 0, 0

      minmax_fill = spcavg_low(0:1, :)
      minmax_fill(0, :) = spcavg_low(i, :) 
      minmax_fill(1, :) = spcavg_high(i, :) 

 
      delete(res@xyLineColors)
        res@gsnXYFillColors = cmap(idc(i), :)
      res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
      bot_plot(i)  = gsn_csm_xy (wks, dof@frq*dlen, minmax_fill, res)  ; Create filled XY plot.
      overlay(plot, bot_plot(i))

   end do

     do i = 4, 4

      minmax_fill = spcavg_low(0:1, :)
       minmax_fill(0, :) = spcavg_low(i, :)
      minmax_fill(1, :) = spcavg_high(i, :)

      res@gsnXYFillColors = cmap(idc(i), :)
      res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
      bot_plot(i)  = gsn_csm_xy (wks, dof@frq*dlen, minmax_fill, res)  ; Create filled XY plot.
      overlay(plot, bot_plot(i))

   end do




    draw(plot)

    res_text               = True
    res_text@txFontHeightF = 0.015
    res_text@txJust        = "CenterLeft"

    res_lines                   = True        ; polyline mods desired

    yy_ar = (/4.0, 5.0, 6.0, 7.0, 8.0/) + 10.0

    ;yy_ar = (/4.0, 5.0, 6.0, 7.0, 8.0/) * 3.0

    txt_arr = (/"1940-69", "1970-99", "2000-29", "2030-59", "2060-89"/)

   ;txt_arr = (/"1970-99", "2060-89"/)


    do i = 0, 4

    res_lines@gsLineColor       = clevs(i, :)   ; line color
    res_lines@gsLineThicknessF  = 5 + 3.0           ; line thicker
    res_lines@gsLineDashPattern = 0   ; dash pattern

    xx = (/1.4,2.0/) + 6.0
    xx = xx
    yy = (/yy_ar(i), yy_ar(i)/)
    gsn_polyline(wks,plot,xx,yy,res_lines)
    gsn_text(wks,plot,txt_arr(i),2.2 + 6.0,yy_ar(i),res_text)

   end do

    ;yy = (/52.5,52.5/)
    ;res_lines@gsLineColor       = line_col2   ; line color
    ;res_lines@gsLineDashPattern = line_pat2   ; dash pattern
    ;gsn_polyline(wks,plot,xx,yy,res_lines)
    ;gsn_text(wks,plot,m2,1900,52.5,res_text)  

   frame(wks) 









