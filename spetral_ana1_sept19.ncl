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
begin
;************************************************
; variable and file handling
;************************************************
   data1  = asciiread("/home/szk0139/CLIMATE_ANA/AMO_DATA.csv", (/163, 12/), "float")
   
   printVarSummary(data1)    ; this is 1856 to 2018 data 

   data1@_FillValue = -99.0

                               ; open netcdf file
   soi  = dim_avg_n(data1, 1)

   dlen = dimsizes(soi)

   print(dlen)


   printVarSummary(soi)

   
 
   print(dim_variance_n(soi, 0))

                                        ; get data
;************************************************
; set function arguments
;************************************************
; detrending opt: 0=>remove mean 1=>remove mean and detrend
  d = 0
; smoothing periodogram: (0 <= sm <= ??.) should be at least 3 and odd
  sm = 1
; percent tapered: (0.0 <= pct <= 1.0) 0.10 common. 
  pct = 0.10
;************************************************
; calculate spectrum
;************************************************
  spec = specx_anal(soi,d,sm,pct)

  splt = specx_ci(spec, 0.05, 0.95)             ; calc confidence interval


;************************************************
; plotting
;************************************************
   wks  = gsn_open_wks("x11","spec")               ; send graphics to PNG file 

   res = True
   ;res@gsnDraw       = False                       ; do not draw
   ;res@gsnFrame      = False                       ; do not advance frame					   ; plot mods desired
   res@tiMainString = "AMO"		           ; title
   res@tiXAxisString = "Frequency (cycles/"+dlen+"year)"  ; xaxis
   res@tiYAxisString = "Variance"                  ; yaxis

   


   res@trXMinF  = 1.0                   ; min value on y-axis
   res@trXMaxF  = 21.0                   ; max value on y-axis


   res@xyLineThicknesses   = (/3.,1.,1.,1./)       ; Define line thicknesses 
   res@xyDashPatterns      = (/0,0,1,1/)           ; Dash patterns   ;plot=gsn_csm_xy(wks,(spec@frq)*dlen,spec@spcx,res)     ; create plot   
   res@xyLineColors        = (/"foreground","blue","blue","blue"/)

   plot = gsn_csm_xy(wks,(spec@frq)*dlen, splt,res)

;***********************************************
end


