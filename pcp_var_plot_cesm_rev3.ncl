; Here I test how many years and/or ensemble members are needed to achieve a stable climatology,
; as measured by a high signal-to-total ratio (see Dirmeyer et al. 2013, eqs 1-4)

; For faster computations, analysis is restricted to North America
 latmin = 10
 latmax = 75
 lonmin = 360-170
 lonmax = 360-50

;****************************************************************************
 fm = addfile("/glade/work/kumar34/prj_data/surf_data/CLM4.0CN.h0.fixed.nc", "r")
;fm = addfile("/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/datafiles/CLM4.0CN.h0.fixed.nc", "r")
 lm = fm->landmask({latmin:latmax},{lonmin:lonmax})
 area = fm->area({latmin:latmax},{lonmin:lonmax})
 landfrac = fm->landfrac({latmin:latmax},{lonmin:lonmax})
 area@_FillValue = 1.0*10^36
 landfrac@_FillValue = 1*0*10^36
 arwt = area
 arwt = area*landfrac/100
 lat = fm->lat({latmin:latmax})
 lon = fm->lon({lonmin:lonmax})


 ;*****************************************

 inDir = "/glade/campaign/cesm/collections/cesmLE/CESM-CAM5-BGC-LE/lnd/proc/tseries/monthly/H2OSOI/"
 
 run = (/"B20TRC5CNBDRD","BRCP85C5CNBDRD"/)

 file1 = systemfunc("ls "+inDir+"b.e11."+run(0)+".f09_g16.002.clm2.h0.H2OSOI.192001-200512.nc")

 ft1 = addfile(file1,"r")

 ;var1 = ft1->H2OSOI(0, 0, {latmin:latmax}, {lonmin:lonmax})  

 var1 = ft1->H2OSOI(0, 0, :, :)


 printVarSummary(var1)


  dirout = "/glade/work/kumar34/hydroclimate_ncc/"

  f2 = addfile(dirout+"/PRECT_CESM-LE_five30yr_period_prect_var_rev2.nc", "r") 

  
  var2 = f2->prect_var

  var2 = (/var2*1000.0/) 

  var2 = mask(var2, ismissing(var1), False)  
  

  print(avg(var2))

  delete(var2@long_name)
  delete(var2@units)

  




  var2_var = var2
  var2_var = (/var2^2.0/)

  var2_avg_var = dim_avg_n_Wrap(var2_var, 0)


  var2_avg = var2_avg_var  
  var2_avg = (/var2_avg_var^0.5/)

  print(avg(var2_avg))

  


  ;print(avg(var2_std))

  var2_diff = var2_avg


  printVarSummary(var2)


  var2_diff_ens = var2
  var2 = mask(var2, var2.ge.0.001, True)
  do i = 0, 4
     var2_diff_ens(:, i, :, :) = ((var2(:, i, :, :) - var2(:, 1, :, :))*100.0/var2(:, 1, :, :))
  end do
 
  var2_diff_ens_std = dim_stddev_n_Wrap(var2_diff_ens, 0)


  






 ;exit

  var2_avg(1, :, :) = mask(var2_avg(1, :, :), var2_avg(1, :, :).ge.0.01, True)


  do i = 0, 4

     ;var2_avg(1, :, :) = mask(var2_avg(1, :, :), var2_avg(1, :, :).ge.0.5, True)

     var2_diff(i, :, :) = (var2_avg(i, :, :) - var2_avg(1, :, :))*100/var2_avg(1, :, :)

  end do

  print(avg(var2_diff))



  var2_sig = var2_diff

  do i = 0, 4

     alpha1 = 100.*(1. - ftest(var2_avg_var(i, :, :), 39*30, var2_avg_var(1, :, :), 39*30, 0))

     var2_sig(i, :, :) = (/alpha1/)

     delete(alpha1)
 
 end do

 ;exit

 ;********************************************************************

  plottype = "png"
  plottype@wkWidth  = 1000
  plottype@wkHeight = 1000



  wks = gsn_open_wks(plottype, "PCP_VAR_Scale_Change_CESM-LE_rev5")

  gsn_define_colormap(wks, "precip_11lev")
 

 res = True
 res1 = True

 res@gsnDraw            = False             ; don't draw yet
 res@gsnFrame           = False             ; don't advance frame yet
 res@gsnMaximize = True


 res@cnInfoLabelOn = False

 ;;res@mpLimitMode = "Corners"
 ;;res@mpLeftCornerLatF = lat(55)
 ;res@mpLeftCornerLonF = lon(143)
 ;;res@mpRightCornerLatF = lat(165)
 ;res@mpRightCornerLonF = lon(145)

 res@mpLimitMode          = "LatLon"
 res@mpMinLatF            = 15
 res@mpMaxLatF            = 72
 res@mpMinLonF            = -360.0+192.0
 res@mpMaxLonF            =  -360.0+308.0


 NA_slat = 15.0

 NA_nlat = 72.0 

 NA_elon = 192.0 
 NA_wlon = 308.0

 res@mpOutlineBoundarySets = "Geophysical"
 res@mpGeophysicalLineColor = "black"

 res@pmTickMarkDisplayMode = "Always"

 res@cnFillOn = True
 res@cnLinesOn = False
 ;res@cnFillMode = "RasterFill"
 res@cnLineLabelsOn = False
  
 res@gsnSpreadColors = True
 
 res@lbOrientation = "Horizontal"

 res@lbBoxLinesOn         = False           ; Turn off labelbar box lines
 res@lbLabelAutoStride    = True            ; Control labelbar labels.
 
 res@lbLabelBarOn = False

 res@lbLabelFontHeightF = 0.010

 ;*******PLEASE CHECK THIS SECTION-A *************************
 ;*******SECTION A******************************************* 
 ; ***first plot by commenting all the res in section A****
 ; *** the decide which maximum and minimum value you want to assign *****
 ; **** accordingly decide the min, max and stide for manual cn level selction *****
 res@cnLevelSelectionMode = "ManualLevels"
 res@cnMinLevelValF =  5.0
 res@cnMaxLevelValF = 30.0
 res@cnLevelSpacingF = 5.0
 ;******************************************

   plots1 = new(1, graphic)

   print(avg(var2_avg(1, :, :)))

   ;exit 
  res@gsnAddCyclic = False


     res@tmYLLabelsOn               = False            ; turn off lat labels
    res@tmXBLabelsOn               = False            ; turn off lon labels
  
    res@tmYROn                     = False
    res@tmYLOn                     = False
    res@tmXTOn                     = False
    res@tmXBOn                     = False


 ;**************************************************

;  region = (/"Southwest","Southeast","South US Plains","North US Plains","Canadian Plains"/)
;  nregs = dimsizes(region)

; latS = (/30.,25.,28.,43.,53./)
; latN = (/40.,35.,38.,49.,61./)
; lonW = (/239.,265.,255.,254.,242./)
; lonE = (/252.,282.,263.,264.,257./)

 region = (/"Southwest", "Southeaset", "Canadian Plains"/)
 nregs = dimsizes(region)

  ;latS = (/30.,25.,53./)
  ;latN = (/40.,35.,61./)
  ;lonW = (/239.,265.,242./)
  ;lonE = (/252.,282.,257./)


  latS = (/30.,25.,49./)
  latN = (/40.,35.,60./)
  lonW = (/239.,265.,248./)
  lonE = (/252.,282.,263./)





 lacoords = new((/nregs,5/),"float")
 locoords = lacoords
 do r = 0,nregs-1
   lacoords(r,:) = (/latS(r),latN(r),latN(r),latS(r),latS(r)/)
   locoords(r,:) = (/lonW(r),lonW(r),lonE(r),lonE(r),lonW(r)/)
 end do



;*****************************************************



  plots1(0) = gsn_csm_contour_map(wks, var2_avg(1, :, :), res)
 

  bxres = True
   bxres@gsLineThicknessF = 2.
    bxres@gsLineColor = "black"

   nps = 8

   boxes = new((/nps*2,nregs/),graphic)

  do r = 0,4
   do pm = 0,0
     boxes(pm,r) = gsn_add_polyline(wks,plots1(pm),locoords(r,:),lacoords(r,:),bxres)
   end do
 end do



  pres = True
  pres@gsnFrame = False
  pres@gsnPanelLabelBar = True
  pres@lbOrientation     = "horizontal"
  ;pres@gsnPanelFigureStrings= (/"16-30 days", "31-60 days", "61-90 days"/) ; add strings to panel
  ;pres@txString = "1850"
 
  pres@amJust = "BottomLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.015
  pres@gsnPanelLeft  = 0.05
  pres@gsnPanelRight = 0.35 
  pres@pmLabelBarWidthF = 0.25

  MAXIMIZE_OUTPUT           = True 

  if(MAXIMIZE_OUTPUT) then
      pres@gsnDraw = False                           ; draw later, in maximize mode
      panelid1 = gsn_panel_return(wks,plots1,(/1,1/),pres) 
  else
      pres@gsnDraw = True
      gsn_panel(wks,plots1,(/1,1/),panres1)               ; now draw as one 
  end if

 plots = new(2, graphic)

 
 dres = res
 dres@gsnLeftString = "" 
 dres@gsnRightString = ""  


 dres@cnFillPalette = "GreenMagenta16"
  
 ;dres@cnFillPalette = "BlueDarkRed18"
 
 dres@cnMinLevelValF =  -25.0
 dres@cnMaxLevelValF =  25.0
 dres@cnLevelSpacingF = 5.0
 dres@cnSpanFillPalette = True

 
 plots(0) = gsn_csm_contour_map(wks, var2_diff(3, :, :), dres)
 plots(1) = gsn_csm_contour_map(wks, var2_diff(4, :, :), dres)

 printVarSummary(var2_diff&lat)
 printVarSummary(var2_diff&lon) 

 NA_slat = 10.0

 NA_nlat = 71.0

 NA_elon = 360.0-170.0
 NA_wlon = 360.0-50.0
 



 ;latwt = latGauWgt(dimsizes(var2_diff&lat), "lat", "wg wgt", "")


 latwt = latRegWgt(var2_diff&lat, "float", 0)
 latwt!0 = "lat"
 latwt&lat = var2_diff&lat

 ;***

 print("numbers print")

 print(wgt_areaave(var2_diff(3, {NA_slat:NA_nlat}, {NA_elon:NA_wlon}), latwt({NA_slat:NA_nlat}), 1.0, 0))

 print(wgt_areaave(var2_diff(4, {NA_slat:NA_nlat}, {NA_elon:NA_wlon}), latwt({NA_slat:NA_nlat}), 1.0, 0)) 
 print(wgt_areaave(var2_diff_ens_std(4, {NA_slat:NA_nlat}, {NA_elon:NA_wlon}), latwt({NA_slat:NA_nlat}), 1.0, 0))


 ;******



   bxres = True
   bxres@gsLineThicknessF = 2.
    bxres@gsLineColor = "black"

   nps = 8

   boxes2 = new((/nps*2,nregs/),graphic)

  do r = 0,4
   do pm = 0,1
     boxes2(pm,r) = gsn_add_polyline(wks,plots(pm),locoords(r,:),lacoords(r,:),bxres)
   end do
 end do






 
 
 ;***************************************************

  ;res2 = True
  ;res2@gsnDraw            = False             ; don't draw yet
  ;res2@gsnFrame           = False             ; don't advance fra

  ;res2@gsnAddCyclic = False
  ;res2@gsnMaximize = True

   ;res2@cnInfoLabelOn       = False       ; turn off info label
   ;res2@cnLinesOn           = False       ; do not draw contour lines
   ;res2@cnLineLabelsOn      = False       ; do not draw contour labels
   ;res2@cnFillScaleF        = 0.7         ; add extra density
   ;res2@cnFillMode = "RasterFill"

   ;res2@cnLevelSelectionMode = "ManualLevels"
   ;res2@cnMinLevelValF =  0.0
   ;res2@cnMaxLevelValF =  3.0
   ;res2@cnLevelSpacingF =  1.0

   opt = True
   opt@gsnShadeFillType = "pattern"      ; pattern fill
   opt@gsnShadeHigh = 17       ; use pattern #2
   ;opt@gsnShadeLow = 3

  res2 = True
  res2@gsnDraw            = False             ; don't draw yet
  res2@gsnFrame           = False             ; don't advance fra
  res2@gsnMaximize = True

   res2@cnInfoLabelOn       = False       ; turn off info label
   res2@cnLinesOn           = False       ; do not draw contour lines
   res2@cnFillOn             = True
   res2@cnLineLabelsOn      = False       ; do not draw contour labels
   res2@gsnAddCyclic = False

   res2@lbLabelBarOn = False
   res2@cnMonoFillColor = False                 ; default color is fground [black]
   res2@cnMonoFillPattern = False              ; want multiple patterns
   res2@cnMonoFillScale   = False              ; want patterns w/ diff densities

   res2@cnLevelSelectionMode = "ManualLevels"
   res2@cnMinLevelValF =  70.0
   res2@cnMaxLevelValF =  100.0
   res2@cnLevelSpacingF =  5.0



   res2@cnFillPatterns    = (/ -1,  -1,  -1, -1, -1, -1, 17, 17/)    ; -1=trasparent
                          ;(/less than 70, less than 75, less than 80, less than 85, less than 90, less than 95, less than 100, gretaer than 100/)
   res2@cnFillScales      = (/1., 1., 1., 1., 1., 1., 1., 1./)
   res2@cnFillColors      = (/"white", "white", "white", "white", "white", "white", "navy", "navy"/)
   ;res2@cnFillColors      = (/"white", "white", "white", "white", "white", "black", "black", "black"/)

   res2@cnFillDotSizeF       = 0.0035



   plots2 = new(2, graphic)


  
   plots2(0) = gsn_csm_contour(wks, var2_sig(3, :, :), res2)
   ;plots2(0) = gsn_contour_shade(plots2(0), -9999.0, 95.0, opt)  ; shade all areas greater than 95 CI

   overlay(plots(0), plots2(0))

   plots2(1) = gsn_csm_contour(wks, var2_sig(4, :, :), res2)
   ;plots2(1) = gsn_contour_shade(plots2(1), -9999.0, 95.0, opt)  ; shade all areas greater than 95 CI

   overlay(plots(1), plots2(1))
 


    pres@gsnPanelLeft  = 0.35
    pres@gsnPanelRight = .95
    delete(pres@gsnPanelFigureStrings)
    ;pres@gsnPanelFigureStrings= (/"16-30 days", "16-30 days", "31-60 days", "31-60 days", "61-90 days", "61-90 days"/)
    pres@pmLabelBarWidthF = 0.5
    txres = True 
    txres@txFontHeightF = 0.025 
    ;drawNDCGrid(wks)


    if(MAXIMIZE_OUTPUT) then
       pres@gsnDraw = False                           ; draw later, in maximize mode
       panelid2 = gsn_panel_return(wks,plots,(/1,2/),pres)
       gsn_text_ndc(wks,"1970-99",0.2,0.64,txres)
       gsn_text_ndc(wks,"2030-59",0.5,0.64,txres)
       gsn_text_ndc(wks,"2060-89",0.8,0.64,txres)
       mres = True
       maximize_output(wks,mres)
    else
       pres@gsnDraw = True
       gsn_panel(wks, plots, (/1,2/), pres)
       frame(wks)
    end if

 exit

