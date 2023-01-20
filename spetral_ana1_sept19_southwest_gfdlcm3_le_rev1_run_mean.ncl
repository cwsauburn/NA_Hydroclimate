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

 modN = 1
 
 f1 = addfile(dirIn+model(modN)+"_1940-2089_30yr_TimeSeries_"+region(regN)+".nc", "r")
 print(f1)

 var_mod = f1->gfdlcm3_soilm(:, :, :)    ; cesm_soilm(periods, ens, time) 
 
 var_syn = f1->gfdlcm3_synt_soilm(:, :, :, :)   ; cesm_synt_soilm(periods, iter, ens, time)  


 ave_yr = 5

 var_mod_run = runave_n_Wrap(var_mod, ave_yr, 0, 2)
 var_syn_run = runave_n_Wrap(var_syn, ave_yr, 0, 3)

;**************************************

 var_mod_VAR = dim_variance_n_Wrap(var_mod, 2)
 var_syn_VAR = dim_variance_n_Wrap(var_syn, 3)

 var_mod_run_VAR = dim_variance_n_Wrap(var_mod_run, 2)
 var_syn_run_VAR = dim_variance_n_Wrap(var_syn_run, 3)

;*******************************************

 var_mod_VAR_sum = dim_sum_n_Wrap(var_mod_VAR, 1)
 var_syn_VAR_sum = dim_sum_n_Wrap(dim_sum_n_Wrap(var_syn_VAR, 2), 1)


 var_mod_run_VAR_sum = dim_sum_n_Wrap(var_mod_run_VAR, 1)
 var_syn_run_VAR_sum = dim_sum_n_Wrap(dim_sum_n_Wrap(var_syn_run_VAR, 2), 1)


;******************************************
print("I am printing model")
 print(var_mod_run_VAR_sum/var_mod_VAR_sum)
print("I am printing synthetic")
 print(var_syn_run_VAR_sum/var_syn_VAR_sum)

exit






;************************************************************* 
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
  

   spcavg = 0.0


   do i = 0, 4
      x = var_mod(i, :, :)
      
      do n = 0, nseg-1
         dof = specx_anal(x(n,:),d,sm,pct)         
         temp1 = dof@spcx
         temp2 = temp1*100/(dim_variance_n(x(n, :), 0)*dlen)
         spcavg(i, :) = spcavg(i, :) + temp2         
      end do
   
   end do

   spcavg = spcavg/nseg

  print(sum(spcavg)/5)

  ;************************************************
  ; plotting
  ;************************************************
  ;    wks  = gsn_open_wks("x11","spec")              ; Opens a ps file
  
    wks  = gsn_open_wks("x11", region(regN)+"_"+model(modN))              ; Opens a ps file



    res = True
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
    res@trYMaxF  = 18.0                   ; max value on y-axis

    res@trXMinF  = 1.0                   ; min value on y-axiss
    res@trXMaxF  = 10.0                   ; max value on y-axis



     cmap = read_colormap_file("MPL_rainbow")
     idc = (/2, 20, 40, 93, 127/)
     clevs = cmap(idc,:)

     res@xyDashPattern  = 0                         ; Make curves all solid

     res@xyLineThicknesses = (/  4.0,   4.0, 4.0, 4.0, 4.0/)          ; make second line thicker
     res@xyLineColors      = clevs

     plot = gsn_csm_xy(wks, dof@frq*dlen, spcavg, res)      ; create plot


       res_text               = True
    res_text@txFontHeightF = 0.015
    res_text@txJust        = "CenterLeft"

    res_lines                   = True        ; polyline mods desired

    yy_ar = (/4.0, 5.0, 6.0, 7.0, 8.0/)

    yy_ar = yy_ar + 8.0

    txt_arr = (/"1940-69", "1970-99", "2000-29", "2030-59", "2060-89"/)


    do i = 0, 4

    res_lines@gsLineColor       = clevs(i, :)   ; line color
    res_lines@gsLineThicknessF  = 5           ; line thicker
    res_lines@gsLineDashPattern = 0   ; dash pattern

       xx = (/1.4,2.0/)
    yy = (/yy_ar(i), yy_ar(i)/)
    gsn_polyline(wks,plot,xx,yy,res_lines)
    gsn_text(wks,plot,txt_arr(i),2.2,yy_ar(i),res_text)

   end do
  
  frame(wks)


