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

 regN = 4 
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

   spcavg = new((/5, floattoint(ntim/2)/), typeof(var_mod), var_mod@_FillValue)

   spcavg_ens = new((/5, floattoint(ntim/2), nseg/), typeof(var_mod), var_mod@_FillValue)  




   spcavg = 0.0


   do i = 0, 4
      x = var_mod(i, :, :)
      
      do n = 0, nseg-1
         dof = specx_anal(x(n,:),d,sm,pct)         
         temp1 = dof@spcx
         temp2 = temp1*100/(dim_variance_n(x(n, :), 0)*dlen)
         spcavg(i, :) = spcavg(i, :) + temp2
         spcavg_ens(i, :, n) = (/temp2/)         
      end do
   
   end do

   spcavg = spcavg/nseg

   print(sum(spcavg)/5)

   spcavg_se = 2.0* dim_stddev_n_Wrap(spcavg_ens, 2)/(nseg^0.5)

 
   spcavg_hf = spcavg(0:1, :)
   spcavg_hf(0, :) = spcavg(1, :)
   spcavg_hf(1, :) = spcavg(4, :)

   
   spcavg_se_hf = spcavg_se(0:1, :)
   spcavg_se_hf(0, :) = spcavg_se(1, :)
   spcavg_se_hf(1, :) = spcavg_se(4, :)

  ;print(spcavg_hf)
 ;exit 
 ;print(spcavg_se_hf)






  ;************************************************
  ; plotting
  ;************************************************
  ;    wks  = gsn_open_wks("x11","spec")              ; Opens a ps file
  
    wks  = gsn_open_wks("x11", region(regN)+"_"+model(modN)+"uncn")              ; Opens a ps file


    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@gsnMaximize = True


    ;res@tiMainString = "Mean Spectra: "+nseg+" segments, dof="+df                   ; title
    ;res@tiXAxisString = "Frequency (cycles/month)"  ; xaxis
    res@tiXAxisString = "Frequency (cycles/"+dlen+"year)"  ; xaxis
    res@tiYAxisString = "Variance"                  ; yaxis
  
    ;splt = specx_ci(df, 0.05, 0.95)                 ; confidence interval
  
    ;plot = gsn_csm_xy(wks, df@frq, splt ,res)      ; create plot

    res@vpHeightF= 0.4                    ; change aspect ratio of plot
    res@vpWidthF = 0.6

    res@trYMinF  = 2.0                   ; min value on y-axis
    res@trYMaxF  = 18.0  + 22.0                 ; max value on y-axis

    res@trXMinF  = 1.0                   ; min value on y-axiss
    res@trXMaxF  = 10.0                   ; max value on y-axis



     cmap = read_colormap_file("MPL_rainbow")
     ;idc = (/2, 20, 40, 93, 127/)
      idc = (/20, 127/)
     clevs = cmap(idc,:)

     res@xyDashPattern  = 0                         ; Make curves all solid

     ;res@xyLineThicknesses = (/  4.0,   4.0, 4.0, 4.0, 4.0/) 
          ; make second line thicker

     res@xyLineThicknesses = (/7.0, 7.0/)
     res@xyLineColors      = clevs

     plot = gsn_csm_xy(wks, dof@frq*dlen, spcavg_hf, res)      ; create plot

      
    minmax_fill = spcavg_se_hf
    minmax_fill(0, :) = spcavg_hf(0, :) - spcavg_se_hf(0, :)
    minmax_fill(1, :) = spcavg_hf(0, :) + spcavg_se_hf(0, :)

     delete(res@xyLineColors)
     res@gsnXYFillColors = "LightBlue"
     res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
     bot_plot  = gsn_csm_xy (wks, dof@frq*dlen, minmax_fill, res)  ; Create filled XY plot.

    minmax_fill = spcavg_se_hf
    minmax_fill(0, :) = spcavg_hf(1, :) - spcavg_se_hf(1, :)
    minmax_fill(1, :) = spcavg_hf(1, :) + spcavg_se_hf(1, :)


    res@gsnXYFillColors = "LightPink"
    res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
    mid_plot  = gsn_csm_xy (wks,dof@frq*dlen,minmax_fill,res)  ; Create another filled XY plot.


   overlay(plot,bot_plot)
   overlay(plot,mid_plot)

   draw(plot)



   res_text               = True
    res_text@txFontHeightF = 0.015
    res_text@txJust        = "CenterLeft"

    res_lines                   = True        ; polyline mods desired

    yy_ar = (/4.0, 5.0, 6.0, 7.0, 8.0/) * 3.0

    ;yy_ar = yy_ar + 9.5

    ;txt_arr = (/"1940-69", "1970-99", "2000-29", "2030-59", "2060-89"/)


   txt_arr = (/"1970-99", "2060-89"/)

    do i = 0, 1

    res_lines@gsLineColor       = clevs(i, :)   ; line color
    res_lines@gsLineThicknessF  = 5           ; line thicker
    res_lines@gsLineDashPattern = 0   ; dash pattern

       xx = (/1.4,2.0/) + 6.0
    yy = (/yy_ar(i), yy_ar(i)/)
    gsn_polyline(wks,plot,xx,yy,res_lines)
    gsn_text(wks,plot,txt_arr(i),2.2+6.0,yy_ar(i),res_text)

   end do
  
  frame(wks)


