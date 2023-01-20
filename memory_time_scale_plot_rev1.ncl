
load "/gpfs/u/home/kumar34/GLEW_ANA/paper1/area_rectlinear_grid.ncl"
;******README for revision*************************
;***In this reviesion for mod 9 only 1st ensemble is used*************
;*************done reading with GLDAS data****************
;******************* adding GPCP precip data ***************************

   url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"

   varname = "H2OSOI"
   vartitle = "soil moisture"
   ocnvar = "ENSO"

   ;dirlnd = "/Projects/SERDP_regimeshifts/data/CESM-LE/lnd/postprocessed/"
   dirlnd = url+"CESM-LE/lnd/postprocessed/"
   filelnd = dirlnd+"NorthAmerica_"+varname+"_CESM-LE_1940-2089_30yr_MemoryENSOregressions.nc"
   fl = addfile(filelnd,"r")
   s2n = fl->sig2noise ;signal-to-noise of 3 regression models




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

 var1 = ft1->H2OSOI(0:4, 0, {latmin:latmax}, {lonmin:lonmax})  

 printVarSummary(var1)


  dirout = "/glade/work/kumar34/hydroclimate_ncc/"

  f2 = addfile(dirout+"/H2OSOI_CESM-LE_five30yr_period_sm_memory_timescale.nc", "r") 

  
  var2 = f2->sm_memory

  print(avg(var2))
  ;exit

  delete(var2@long_name)
  delete(var2@units)


  var2_avg = dim_avg_n_Wrap(var2, 0)

  print(avg(var2_avg))

  var2_std = dim_stddev_n_Wrap(var2, 0)


  print(avg(var2_std))

  var2_diff = var2_avg


  printVarSummary(var2)





 ;exit

  var2_avg(1, :, :) = mask(var2_avg(1, :, :), var2_avg(1, :, :).ge.0.5, True)


  do i = 0, 4

     ;var2_avg(1, :, :) = mask(var2_avg(1, :, :), var2_avg(1, :, :).ge.0.5, True)

     var2_diff(i, :, :) = (var2_avg(i, :, :) - var2_avg(1, :, :))*100/var2_avg(1, :, :)

  end do

  print(avg(var2_diff))




  var2_diff_ens = var2
  var2 = mask(var2, var2.ge.0.1, True)

  

  ;do i = 0, 4
  ;   var2_diff_ens(:, i, :, :) = ((var2(:, i, :, :) - var2(:, 1, :, :))*100.0/var2(:, 1, :, :))
  ;end do

  
  do k = 0, 39

   do i = 0, 4
     var2_diff_ens(k, i, :, :) = ((var2(k, i, :, :) - var2(k, 1, :, :))*100.0/var2_avg(1, :, :))
   end do

   end do


  var2_diff_ens_std = dim_stddev_n_Wrap(var2_diff_ens, 0)



 ;exit 



  var2_sig = var2_diff

  do i = 0, 4

     alpha1 = 100.*(1. - ttest(var2_avg(i, :, :), var2_std(i, :, :)^2.0, 39, var2_avg(1, :, :), var2_std(1, :, :)^2.0, 39, False, False))

     var2_sig(i, :, :) = (/alpha1/)

     delete(alpha1)
 
  end do

 ;exit

 ;********************************************************************


  plottype = "x11"
 plottype@wkWidth  = 1000
 plottype@wkHeight = 1000




  wks = gsn_open_wks(plottype, "SM_Memory_Scale_Change_CESM_rev1")

  gsn_define_colormap(wks, "precip_11lev")
 

   res = True
 res1 = True

 res@gsnDraw            = False             ; don't draw yet
 res@gsnFrame           = False             ; don't advance frame yet

 res@cnInfoLabelOn = False

 ;;res@mpLimitMode = "Corners"
 ;;res@mpLeftCornerLatF = lat(55)
 ;res@mpLeftCornerLonF = lon(143)
 ;;res@mpRightCornerLatF = lat(165)
 ;res@mpRightCornerLonF = lon(145)

 res@mpLimitMode          = "LatLon"
 res@mpMinLatF            = 10
 res@mpMaxLatF            = 75
 res@mpMinLonF            = -170
 res@mpMaxLonF            =  -50
 res@gsnMaximize = True

   ;res@mpLimitMode       = "LatLon"
   ;res@mpMinLatF         = min(s2n&lat)
   ;res@mpMaxLatF         = max(s2n&lat)
   ;res@mpMinLonF         = min(s2n&lon)
   ;res@mpMaxLonF         = max(s2n&lon)
   ;res@mpCenterLonF      = dim_median(s2n&lon) ;(min(lon)+max(lon))/2
   







res@mpOutlineBoundarySets = "Geophysical"
   res@mpGeophysicalLineColor = "black"

   res@mpLimitMode          = "LatLon"
 res@mpMinLatF            = 15
 res@mpMaxLatF            = 72
 res@mpMinLonF            = -360 + 192
 res@mpMaxLonF            =  -360 + 308




   a_lat = min(s2n&lat)
   b_lat = max(s2n&lat)

   c_lon = min(s2n&lon)
   d_lon = max(s2n&lon)


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

  latwt = latRegWgt(s2n&lat, "float", 0)

  latwt!0 = "lat"
  latwt&lat = s2n&lat

 ; print("printing the skill for year 1")

  P102 = new((/3, 4/), "float", -9999.0)






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
 res@cnMinLevelValF =  3.0
 res@cnMaxLevelValF = 21.0
 res@cnLevelSpacingF = 3.0
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


  

  plots1(0) = gsn_csm_contour_map(wks, var2_avg(1, :, :), res)
 

  do r = 0, 2
    P102(r, 0) = wgt_areaave(var2_avg(1, {latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)
  end do



    nps = 5

   nregs = 3

      bxres = True
    bxres@gsLineThicknessF = 4.
    bxres@gsLineColor = "black"



   ;delete(boxes)
   boxes = new((/nps*2,nregs/),graphic)

  do r = 0,2
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
  pres@gsnMaximize = True
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
  
 
 dres@cnMinLevelValF =  -15.0
 dres@cnMaxLevelValF =  15.0
 dres@cnLevelSpacingF = 3.0
 dres@cnSpanFillPalette = True


 
 plots(0) = gsn_csm_contour_map(wks, var2_diff(3, :, :), dres)
 plots(1) = gsn_csm_contour_map(wks, var2_diff(4, :, :), dres)

  do r = 0, 2
     P102(r, 1) = wgt_areaave(var2_diff(3, {latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)
     P102(r, 2) = wgt_areaave(var2_diff(4, {latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)
     P102(r, 3) = wgt_areaave(var2_diff_ens_std(4, {latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)


  end do

    nps = 5

   nregs = 3

      bxres = True
    bxres@gsLineThicknessF = 4.
    bxres@gsLineColor = "black"



   ;delete(boxes)
   boxes2 = new((/nps*2,nregs/),graphic)

  do r = 0,2
   do pm = 0,1
     boxes2(pm,r) = gsn_add_polyline(wks,plots(pm),locoords(r,:),lacoords(r,:),bxres)
   end do
 end do


 opt102 = True
 opt102@fout = "CESM-LE_sm_memory_scale_std.txt"
 fmtx = "4f15.8"
 write_matrix(P102, fmtx, opt102)

 
 
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



   res2@cnFillPatterns    = (/ -1,  -1,  -1, -1, -1, 17, 17, 17/)    ; -1=trasparent
                          ;(/less than 70, less than 75, less than 80, less than 85, less than 90, less than 95, less than 100, gretaer than 100/)
   res2@cnFillScales      = (/1., 1., 1., 1., 1., 1., 1., 1./)
   res2@cnFillColors      = (/"white", "white", "white", "white", "white", "navy", "navy", "navy"/)
   ;res2@cnFillColors      = (/"white", "white", "white", "white", "white", "black", "black", "black"/)

   res2@cnFillDotSizeF       = 0.0035



   plots2 = new(2, graphic)


  
   plots2(0) = gsn_csm_contour(wks, var2_sig(3, :, :), res2)
   ;plots2(0) = gsn_contour_shade(plots2(0), -9999.0, 95.0, opt)  ; shade all areas greater than 95 CI

   overlay(plots(0), plots2(0))

   plots2(1) = gsn_csm_contour(wks, var2_sig(4, :, :), res2)
   ;plots2(1) = gsn_contour_shade(plots2(1), -9999.0, 95.0, opt)  ; shade all areas greater than 95 CI

   overlay(plots(1), plots2(1))
 

    pres@gsnMaximize = True
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

