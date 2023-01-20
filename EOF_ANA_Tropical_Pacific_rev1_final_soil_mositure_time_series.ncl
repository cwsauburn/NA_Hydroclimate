inDir = "/glade/work/kumar34/GlobalSST/"

f1 = addfile(inDir+"/Land_and_Ocean_LatLong1.nc", "r")

st1 = (1950 - 1850)* 12

end1 = (2022 - 1850) * 12 -1

lm = f1->land_mask

sst = f1->temperature(st1:end1, :, :)

;sst := tofloat(sst)

;sst = mask(sst, lm.eq.0, True)
 sst = mask(sst, lm.eq.1, False)

sst@missing_value = sst@_FillValue

 nNaN = num(isnan_ieee(sst))
 print(nNaN)
 if (nNaN.gt.0) then
     sst = where(isnan_ieee(sst), sst@_FillValue, sst)
 end if

 nBig = num(sst.gt.1e10)

 print(nBig)
 if (nBig.gt.0) then
     sst = where(sst.gt.1e30, sst@_FillValue, sst)
 end if





sst = lonFlip(sst)
lm = lonFlip(lm)
printVarSummary(sst)

var1 = sst

var1 = runave_n_Wrap(var1, 12, 0, 0)


var1_dtrend = var1


nyr = floattointeger(dimsizes(var1&time)/12)

wgt1 = latRegWgt(var1&latitude, "float", 0)
  
temp1 = var1(0:nyr-1, :, :)

;exit

do i = 0, 11
   do j = 0, nyr-1
      temp1(j, :, :) = var1(j*12+i, :, :)
   end do

   temp1_dtrend = temp1
 
   temp1_dtrend = dtrend_msg_n(ispan(1, nyr, 1), temp1, True, False, 0)  

   ;print(avg(temp1_dtrend))

   do j = 0, nyr-1
      temp4 = temp1_dtrend(j, :, :)
      temp4 = mask(temp4, lm.eq.0, True)
      ;ocn_avg = 0.0
      ocn_avg = wgt_areaave(temp4, wgt1, 1.0, 0)
      if (.not.ismissing(ocn_avg)) then
          var1_dtrend(j*12+i, :, :) = (/temp1_dtrend(j, :, :)-ocn_avg/)
      else
         var1_dtrend(j*12+i, :, :) = (/temp1_dtrend(j, :, :)/)

      end if


 
   end do

end do

;exit


;**********************

;**Use  deterned data*********************
 var1 = var1_dtrend

;*****************************

 lonRegL = (/45.0, 0.5/)
 lonRegR = (/275.0, 359.5/)

 ;lonRegL = (/0.5, 0.5/)
 ;lonRegR = (/359.5, 359.5/)



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

 var1_sst = var1


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

 eof_ts1@_FillValue = -9999.0

 printVarSummary(eof_ts1)

 eof_ts2 = tofloat(eof_ts1)


 P1 = new((/dimsizes(var1&time), 2/), "double", -9999.0)
 P1 = P1@_FillValue

 P1(:, 0) = (/var1&time/)
 P1(:, 1) = eof_ts2(0, :)

 opt1 = True
 opt1@fout = "Tropical_IndoPacific_30S_30N_PC1_1950_2021_EE_Earth_final.txt"
 fmtx = "2f15.8"
 write_matrix(P1, fmtx, opt1) 


  data2 = asciiread("./NINO34_NOAA_MONTHLY.csv", (/65, 13/), "float")

  printVarSummary(data2)

  enso_noaa = P1(0:779, 1) 

  do i = 0, 64

    enso_noaa(i*12:i*12+11) = data2(i, 1:12)

  end do
   enso_noaa = runave_n_Wrap(enso_noaa, 12, 0, 0)


 print(escorc(enso_noaa, P1(0:779, 1)))





 ;print(esacr(eof_ts1(0, :), 24)) 







;***********************************************************


 ;****ploting



  plottype = "png"
 ;plottype@wkWidth  = 2000
 ;plottype@wkHeight = 2000





 wks = gsn_open_wks(plottype,"Berkley_Earth_SST_EOF_1950_2021_Soil_Mositure_TS")

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
  res@mpMinLatF            = latRegL(0)        ; zoom in on map
  res@mpMaxLatF            = latRegH(0)
  res@mpMinLonF            = lonRegL(0)
  res@mpMaxLonF            = lonRegR(0)

  res@cnFillOn             = True         ; turn on color fill
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@lbLabelBarOn         = False        ; turn off individual lb's


  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels    = (/-0.25, -0.20, -0.15, -0.10, -0.05, -0.03, -0.01, 0.0, 0.01, 0.03, 0.05, 0.10, 0.15, 0.20, 0.25/)*0.1

;res@mpProjection      = "Robinson"       ; choose projection
 res@mpDataBaseVersion    = "LowRes"    ; better map outlines
 res@mpGridAndLimbOn      = True           ; turn on limb and grid lines
 res@mpGridLineColor      = -1             ; ...but don't draw grid lines...
 res@mpPerimOn            = False          ; turn off map perimeter
 res@mpFillDrawOrder      = "PreDraw"     ; draw map fill last
 res@mpOutlineOn            = True
  ;res@mpPerimOn = True
  ;res@mpLandFillColor       = 2
 res@gsnMaximize = True

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
  res@mpLabelFontHeightF  = 0.2 


  resP                     = True         ; modify the panel plot
  resP@gsnMaximize         = True         ; large format
  resP@gsnPanelLabelBar    = True         ; add common colorbar
  resP@lbLabelAutoStride   = True         ; auto stride on labels
  resP@pmLabelBarWidthF   = 0.8    ; label bar width
  resP@gsnPanelMainString = "Linear Trend Removed from EC Earth 1950-2021"
  resP@gsnPanelFigureStringsFontHeightF = 0.2
 
  res@gsnScale  = True        ; force text scaling
 
  k = 0
  
  do n = 0, neof-1
     res@gsnAddCyclic = True
     res@gsnLeftString  = "EOF "+(n+1)+" Berkeley Earth" 
     res@gsnRightString = sprintf("%5.1f", eof1@pcvar(n)) +"%"
     plot(n)=gsn_csm_contour_map_ce(wks,eof1(n,:,:),res)
  end do

  ;gsn_panel(wks,plot(0:2),(/neof,1/),resP)

  gsn_panel(wks,plot(0),(/1,1/),resP)


  ;*****second plot


  flis = systemfunc("ls /glade/scratch/yanand/ERA5_SM_NorthAmerica/ERA52SubX_SM_monthly_NorthAmerica*.nc")
  print(flis)

 ; exit

  f2 = addfiles(flis, "r")

  ListSetType(f2, "cat")        ; concatenate (=default)


  sm1 = f2[:]->swvl1(:, :, :)

  sm2 = f2[:]->swvl2(:, :, :)

  sm3 = f2[:]->swvl3(:, :, :)

  sm = sm1

  sm = (/(sm1 * 7.0 + sm2 * 21 + sm3 * 72)/100.0/)

  sm@missing_value = sm@_FillValue

 nNaN = num(isnan_ieee(sm))
 print(nNaN)
 if (nNaN.gt.0) then
     sm = where(isnan_ieee(sm), sm@_FillValue, sm)
 end if

 nBig = num(sm.gt.1e10)

 print(nBig)
 if (nBig.gt.0) then
     sm = where(sm.gt.1e30, sm@_FillValue, sm)
 end if


printVarSummary(sm)

 ;exit


 sm_climo = clmMonTLL(sm(12:371, :, :))

 sm_anom = calcMonAnomTLL(sm, sm_climo)



;**************************************************

  
  var1 := sm_anom


  var1 = runave_n_Wrap(var1, 12, 0, 0)


;  elete(H2OSOI_ANOM)
 delete(sm_anom)

  nyr = floattointeger(dimsizes(var1&time)/12)

  var1_dtrend := var1

  temp1 := var1(0:nyr-1, :, :)

  do i = 0, 11
     do j = 0, nyr-1
        temp1(j, :, :) = var1(j*12+i, :, :)
     end do

     temp1_dtrend := temp1
     temp1_dtrend = dtrend_msg_n(ispan(1, nyr, 1), temp1, False, False, 0)
     do j = 0, nyr-1
        var1_dtrend(j*12+i, :, :) = (/temp1_dtrend(j, :, :)/)
     end do

  end do

 sm_anom = var1_dtrend

 ;*********** 

 latS := (/30.,25.,53./)
 latN := (/40.,35.,61./)
 lonW := (/239.,265.,242./)
 lonE := (/252.,282.,257./)

 nregs = 3

 lacoords = new((/nregs,5/),"float")
 locoords = lacoords
 do r = 0,nregs-1
   lacoords(r,:) = (/latS(r),latN(r),latN(r),latS(r),latS(r)/)
   locoords(r,:) = (/lonW(r),lonW(r),lonE(r),lonE(r),lonW(r)/)
 end do

  latwt = latRegWgt(sm&lat, "float", 0)

  latwt!0 = "lat"
  latwt&lat = sm&lat

;  print("printing the skill for year 1")

 P2 = new((/864, 3/), "float", -9999.0)

 
  do r = 0, 2
     P2(:, r) =      wgt_areaave(sm_anom(:, {latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)
  end do


  P2_std = dim_stddev_n_Wrap(P2, 0)

  do r = 0, 2
     P2(:, r) = P2(:, r) / P2_std(r)
  end do 



  ;do r = 0, 2
  ;   eof_ts1(r, :) = (/P2(:, r)/)
  ;end do








;*******************************************
; EOF time series  [bar form]


 resb = True
 resb@gsnScale = True  
 resb@gsnDraw = False
 resb@gsnFrame = False


 resb@vpHeightF = 0.2        ; Changes the aspect ratio
 resb@vpWidthF  = 0.80
 ;resb@gsnMaximize = True

 resb@trYMinF   = -3.0        ; min value on y-axis
 resb@trYMaxF   =  3.0        ; max value on y-axis

 ;resb@tiYAxisString = "Anomalies" ; y-axis label      
 ;resb@tiXAxisString = ""
 ;resb@tiMainString  = "Darwin Southern Oscillation Index" ; title

 resb@gsnYRefLine           = 0.              ; reference line   
 resb@gsnXYBarChart         = True            ; create bar chart 
 resb@gsnAboveYRefLineColor = "deepskyblue"           ; above ref line fill red
 resb@gsnBelowYRefLineColor = "brown"          ; below ref line fill blue
 resb@gsnXYBarChartOutlineOnly = False
 resb@gsnXYBarChartOutlineThicknessF = -1.
 resb@xyLineColors = (/"deepskyblue","brown"/)         ; colors for the above/blow

 
  resb@tiYAxisString = "normalized"  

  plots1 = new(neof, graphic)


  do n=0,neof-1
      ;resb@gsnLeftString  = "EOF "+(n+1)+" Berkeley Earth"
      ;resb@gsnRightString = sprintf("%5.1f", eof1@pcvar(n)) +"%"
      year1 = var1_sst&time
      ;plots1(n) = gsn_csm_xy (wks,year1,eof_ts1(n,:),resb)
       plots1(n) = gsn_csm_xy (wks,year1,P2(:, n),resb)
  end do

  resV = True
  resV@gsnDraw = False
  resV@gsnFrame = False  
  resV@xyLineThicknesses = 5.0
  resV@xyLineColor = "Black"
  resV@xyDashPatterns = 2
  
  x1 = fspan(1, 10, 10)
  x1 = 1977.5
  y1 = fspan(-3.0, 3.0, 10)

  plotv1 = gsn_csm_xy(wks, x1, y1, resV)
  overlay(plots1(0), plotv1)

  x1 = 1999.0
  plotv2 = gsn_csm_xy(wks, x1, y1, resV)
  overlay(plots1(0), plotv2)




   rtsP                      = True            ; modify the panel plot
   rtsP@gsnMaximize          = True            ; large format
   gsn_panel(wks,plots1(0),(/1,1/),rtsP)
  
   exit








  rts           = True
  rts@gsnDraw   = False       ; don't draw yet
  rts@gsnFrame  = False       ; don't advance frame yet
  ;rts@gsnScale  = True        ; force text scaling

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
  rts@trYMinF  = -1.5                   ; min value on y-axis
  rts@trYMaxF  =  1.5                   ; max value on y-axis

; panel plot only resources
  rtsP                      = True            ; modify the panel plot
  rtsP@gsnMaximize          = True            ; large format
  ;rtsP@txString             = "pi_Control_CMIP5(21X200)_CESM-LE(1801)_LIM(11000)_OBS(105)"


  do n=0,neof-1
     rts@gsnLeftString  = "EOF "+(n+1)+" Berkeley Earth"
     rts@gsnRightString = sprintf("%5.1f", eof1@pcvar(n)) +"%"
     year1 = var1&time
      plot(n) = gsn_csm_xy (wks,year1,eof_ts1(n,:),rts)
  end do

  ;gsn_panel(wks,plot(0:2),(/neof,1/),rtsP)     ; now draw as one plot

  gsn_panel(wks,plot(0),(/1,1/),rtsP)

  draw(plot(0))
  frame(wks)























 











 









