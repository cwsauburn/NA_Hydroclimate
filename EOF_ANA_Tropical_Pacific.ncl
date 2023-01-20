inDir = "/glade/work/kumar34/GlobalSST/"

f1 = addfile(inDir+"/Land_and_Ocean_LatLong1.nc", "r")

st1 = (1950 - 1850)* 12

end1 = (2022 - 1850) * 12 -1

lm = f1->land_mask

sst = f1->temperature(st1:end1, :, :)

sst = mask(sst, lm.eq.0, True)

sst@missing_value = sst@_FillValue

sst = lonFlip(sst)
lm = lonFlip(lm)
printVarSummary(sst)

var1 = sst

var1 = runave_n_Wrap(var1, 12, 0, 0)


var1_dtrend = var1


nyr = floattointeger(dimsizes(var1&time)/12)



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

;*****************************

 lonRegL = (/45.0, 0.5/)
 lonRegR = (/275.0, 359.5/)

 latRegL = (/-30.0, -89.0/)
 latRegH = (/30.0, 89.0/)

  SI = 0

  obs_mask = lm
  obs_mask = 1.0

  obs_mask1 = obs_mask
  obs_mask1 = 0.0
  obs_mask1({latRegL(SI):latRegH(SI)}, {lonRegL(SI):lonRegR(SI)}) = 1.0

  obs_mask = mask(obs_mask, obs_mask1.eq.1.0, True)

  var1 = mask(var1, obs_mask.eq.1.0, True)



;**********************
 

;**Use  deterned data*********************


 neof = 3

 optEOF = True

 optEOF@jopt = 0

 optEOF@pcrit = 50

 optETS = False

;*******************************


 wgt1 = latRegWgt(var1&latitude, "float", 0)


 var1_wgt = var1

 var1_wgt = (/var1 *conform_dims(dimsizes(var1), wgt1, 1)/)



 eof1 = eofunc_n_Wrap(var1_wgt, neof, optEOF, 0)

 printVarSummary(eof1)

 eof_ts1 = eofunc_ts_n_Wrap(var1_wgt, eof1, optETS, 0)



 printVarSummary(eof_ts1)


 print(esacr(eof_ts1(0, :), 24)) 







;***********************************************************


 ;****ploting


 wks = gsn_open_wks("x11","EC_Earth_SST_EOF_1950_2021")

 gsn_define_colormap(wks, "BlueWhiteOrangeRed")

 plot = new(4*neof,graphic)                ; create graphic array

 ; EOF patterns

  res                      = True
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@gsnMaximize = True

   res@gsnSpreadColors      = True         ; spread out color table

  res@gsnAddCyclic         = True        ; plotted dataa are not cyclic

  res@mpFillOn             = False        ; turn off map fill

  res@mpLimitMode          = "LatLon"
  res@mpMinLatF            = latRegL(1)        ; zoom in on map
  res@mpMaxLatF            = latRegH(1)
  res@mpMinLonF            = lonRegL(1)
  res@mpMaxLonF            = lonRegR(1)

  res@cnFillOn             = True         ; turn on color fill
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@lbLabelBarOn         = False        ; turn off individual lb's


  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels    = (/-0.25, -0.20, -0.15, -0.10, -0.05, -0.03, -0.01, 0.0, 0.01, 0.03, 0.05, 0.10, 0.15, 0.20, 0.25/)*0.1

  res@mpProjection      = "Robinson"       ; choose projection
 res@mpDataBaseVersion    = "LowRes"    ; better map outlines
 res@mpGridAndLimbOn      = True           ; turn on limb and grid lines
 res@mpGridLineColor      = -1             ; ...but don't draw grid lines...
 res@mpPerimOn            = False          ; turn off map perimeter
 res@mpFillDrawOrder      = "PreDraw"     ; draw map fill last
 res@mpOutlineOn            = True
  ;res@mpPerimOn = True
  ;res@mpLandFillColor       = 2

  res@mpCenterLonF = 180.0

  res@cnFillOn = True
  res@cnLinesOn = False
  ;res@cnFillMode = "RasterFill"
  res@cnLineLabelsOn = False
  res@cnInfoLabelOn = False
  res@cnMissingValFillColor = "gray"
  res@lbBoxLinesOn         = False           ; Turn off labelbar

 ;res@gsnAddCyclic = True
  res@gsnAddCyclic = True
 


  resP                     = True         ; modify the panel plot
  resP@gsnMaximize         = True         ; large format
  resP@gsnPanelLabelBar    = True         ; add common colorbar
  resP@lbLabelAutoStride   = True         ; auto stride on labels
  resP@pmLabelBarWidthF   = 0.8    ; label bar width
  resP@gsnPanelMainString = "Linear Trend Removed from EC Earth 1950-2021"
  
  k = 0
  
  do n = 0, neof-1
     res@gsnAddCyclic = True
     res@gsnLeftString  = "EOF "+(n+1)+" EC EARTH" 
     res@gsnRightString = sprintf("%5.1f", eof1@pcvar(n)) +"%"
     plot(n)=gsn_csm_contour_map_ce(wks,eof1(n,:,:),res)
  end do

  gsn_panel(wks,plot(0:2),(/neof,1/),resP)


  ;*****second plot
;*******************************************
; EOF time series  [bar form]

  rts           = True
  rts@gsnDraw   = False       ; don't draw yet
  rts@gsnFrame  = False       ; don't advance frame yet
  rts@gsnScale  = True        ; force text scaling

; these four rtsources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.

  rts@vpHeightF = 0.40        ; Changes the aspect ratio
  rts@vpWidthF  = 0.85
  rts@vpXF      = 0.10        ; change start locations
  rts@vpYF      = 0.75        ; the plot


  rts@tiYAxisString = "normalized"                    ; y-axis label

  rts@gsnYRefLine           = 0.              ; reference line
  rts@gsnXYBarChart         = True            ; create bar chart
  rts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
  rts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
  rts@trYMinF  = -2.0                   ; min value on y-axis
  rts@trYMaxF  =  2.0                   ; max value on y-axis

; panel plot only resources
  rtsP                      = True            ; modify the panel plot
  rtsP@gsnMaximize          = True            ; large format
  rtsP@txString             = "pi_Control_CMIP5(21X200)_CESM-LE(1801)_LIM(11000)_OBS(105)"


  do n=0,neof-1
     rts@gsnLeftString  = "EOF "+(n+1)+" EC EARTH"
     rts@gsnRightString = sprintf("%5.1f", eof1@pcvar(n)) +"%"
     year1 = var1&time
      plot(n) = gsn_csm_xy (wks,year1,eof_ts1(n,:),rts)
  end do

  gsn_panel(wks,plot(0:2),(/neof,1/),rtsP)     ; now draw as one plot


























 











 









