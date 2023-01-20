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


;*********************************************************************

  fm = addfile("/glade/work/kumar34/prj_data/surf_data/CLM4.0CN.h0.fixed.nc", "r")
  lm = fm->landmask
  area = fm->area
  printVarSummary(area)
  area = mask(area, lm.eq.1, True)

  landfrac = fm->landfrac
  printVarSummary(landfrac)
  area@_FillValue = 1.0*10^36
  landfrac@_FillValue = 1.0*10^36
  arwt = area
  arwt = area*landfrac/100
  lat = fm->lat
  lon = fm->lon
  print(sum(area))
  print(sum(arwt))

  delete(area)
  area = area_global_rectilinear_grid(lat, lon, False)
  area = mask(area, lm.eq.1, True)

  printVarSummary(area)

; Script to plot from netcdf file we just created 


f = addfile("/glade/work/kumar34/DPLE-DATA/DPLE_SM_ana_10_years_v2.nc","r")
;f2 = addfile("/gpfs/fs1/work/azpandey/CESM-DPLE/DPLE_SM_ana_v3_25_49.nc","r")
;f3 = addfile("/gpfs/fs1/work/azpandey/CESM-DPLE/DPLE_SM_ana_v3_50_74.nc","r")
;f4 = addfile("/gpfs/fs1/work/azpandey/CESM-DPLE/DPLE_SM_ana_v3_75_99.nc","r")
;f5 = addfile("/gpfs/fs1/work/azpandey/CESM-DPLE/DPLE_SM_ana_v3_99_121.nc","r") 

f1 = addfile("/glade/work/azpandey/CESM-DPLE/H2OSOI/b.e11.BDP.f09_g16.2015-11.028.clm2.h0.H2OSOI.201511-202512.nc","r")



lat = f1->lat
lon = f1->lon

tsa_DL_S2Tvar1 = f->tsa_DL_S2Tvar
tsa_DL_S2Nvar1 = f->tsa_DL_S2Nvar
tsa_DL_S2Tvar_sig1 = f->tsa_DL_S2Tvar_sig


;tsa_DL_S2Tvar2 = f2->tsa_DL_S2Tvar
;tsa_DL_S2Nvar2 = f2->tsa_DL_S2Nvar
;tsa_DL_S2Tvar_sig2 = f2->tsa_DL_S2Tvar_sig

;tsa_DL_S2Tvar3 = f3->tsa_DL_S2Tvar
;tsa_DL_S2Nvar3 = f3->tsa_DL_S2Nvar
;tsa_DL_S2Tvar_sig3 = f3->tsa_DL_S2Tvar_sig

;tsa_DL_S2Tvar4 = f4->tsa_DL_S2Tvar
;tsa_DL_S2Nvar4 = f4->tsa_DL_S2Nvar
;tsa_DL_S2Tvar_sig4 = f4->tsa_DL_S2Tvar_sig

;tsa_DL_S2Tvar5 = f5->tsa_DL_S2Tvar(1::,:,:)
;tsa_DL_S2Nvar5 = f5->tsa_DL_S2Nvar(1::,:,:)
;tsa_DL_S2Tvar_sig5 = f5->tsa_DL_S2Tvar_sig(1::,:,:)


;tsa_DL_S2Tvar_a1 = array_append_record(tsa_DL_S2Tvar1,tsa_DL_S2Tvar2, 0 )
;tsa_DL_S2Tvar_a2 = array_append_record(tsa_DL_S2Tvar3,tsa_DL_S2Tvar4, 0)
;tsa_DL_S2Tvar_a3 = array_append_record(tsa_DL_S2Tvar_a1,tsa_DL_S2Tvar_a2, 0)
;tsa_DL_S2Tvar    =  array_append_record(tsa_DL_S2Tvar_a3,tsa_DL_S2Tvar5, 0)

;tsa_DL_S2Nvar_a1 = array_append_record(tsa_DL_S2Nvar1,tsa_DL_S2Nvar2, 0 )
;tsa_DL_S2Nvar_a2 = array_append_record(tsa_DL_S2Nvar3,tsa_DL_S2Nvar4, 0)
;tsa_DL_S2Nvar_a3 = array_append_record(tsa_DL_S2Nvar_a1,tsa_DL_S2Nvar_a2, 0)
;tsa_DL_S2Nvar    =   array_append_record(tsa_DL_S2Nvar_a3,tsa_DL_S2Nvar5, 0)

;tsa_DL_S2Tvar_sig_a1 = array_append_record(tsa_DL_S2Tvar_sig1,tsa_DL_S2Tvar_sig2, 0 )
;tsa_DL_S2Tvar_sig_a2 = array_append_record(tsa_DL_S2Tvar_sig3,tsa_DL_S2Tvar_sig4, 0)
;tsa_DL_S2Tvar_sig_a3 = array_append_record(tsa_DL_S2Tvar_sig_a1,tsa_DL_S2Tvar_sig_a2, 0)
;tsa_DL_S2Tvar_sig    =   array_append_record(tsa_DL_S2Tvar_sig_a3,tsa_DL_S2Tvar_sig5, 0)

tsa_DL_S2Tvar = tsa_DL_S2Tvar1
tsa_DL_S2Nvar = tsa_DL_S2Nvar1
tsa_DL_S2Tvar_sig = tsa_DL_S2Tvar_sig1


;printVarSummary(tsa_DL_S2Tvar)
;printVarSummary(tsa_DL_S2Nvar)
;printVarSummary(tsa_DL_S2Tvar_sig)

;exit






;tsa_DL_S2Tvar = f->tsa_DL_S2Tvar

tsa_DL_S2Tvar!1 = "lat"
tsa_DL_S2Tvar!2 = "lon"
tsa_DL_S2Tvar&lat = lat
tsa_DL_S2Tvar&lon = lon

;tsa_DL_S2Nvar = f->tsa_DL_S2Nvar
tsa_DL_S2Nvar!1 = "lat"
tsa_DL_S2Nvar!2 = "lon"
tsa_DL_S2Nvar&lat = lat
tsa_DL_S2Nvar&lon = lon

;tsa_DL_S2Tvar_sig = f->tsa_DL_S2Tvar_sig
tsa_DL_S2Tvar_sig!1 = "lat"
tsa_DL_S2Tvar_sig!2 = "lon"
tsa_DL_S2Tvar_sig&lat = lat
tsa_DL_S2Tvar_sig&lon = lon

printVarSummary(tsa_DL_S2Tvar)
printMinMax(tsa_DL_S2Tvar, True)
printVarSummary(tsa_DL_S2Nvar)
printMinMax(tsa_DL_S2Nvar, True)
printVarSummary(tsa_DL_S2Tvar_sig)
printMinMax(tsa_DL_S2Tvar_sig, True)

;exit

; Take the avergae of 3 months

tsa_DL_S2Tvar_3mon = new((/10,192,288/),"float",9999.0)
tsa_DL_S2Nvar_3mon = new((/10,192,288/),"float",9999.0)
tsa_DL_S2Tvar_sig_3mon = new((/10,192,288/),"float",9999.0)

;do i = 0,9

;tsa_DL_S2Tvar_3mon(i,:,:) = dim_avg_n_Wrap(tsa_DL_S2Tvar(i*3+1:i*3+3,:,:),0)
;tsa_DL_S2Nvar_3mon(i,:,:) = dim_avg_n_Wrap(tsa_DL_S2Nvar(i*3+1:i*3+3,:,:),0)
;tsa_DL_S2Tvar_sig_3mon(i,:,:) = dim_avg_n_Wrap(tsa_DL_S2Tvar_sig(i*3+1:i*3+3,:,:),0)

;end do

tsa_DL_S2Tvar_3mon = tsa_DL_S2Tvar
tsa_DL_S2Nvar_3mon = tsa_DL_S2Nvar
tsa_DL_S2Tvar_sig_3mon = tsa_DL_S2Tvar_sig

printVarSummary(tsa_DL_S2Tvar_3mon)
printVarSummary(tsa_DL_S2Nvar_3mon)
printVarSummary(tsa_DL_S2Tvar_sig_3mon)

;exit



 tsa_DL_S2Tvar_3mon_R2 = tsa_DL_S2Tvar_3mon

 temp11 = tsa_DL_S2Tvar_3mon
 temp11 = 1.0

 T1 = temp11

 temp13 = temp11
 temp13 = temp11 - tsa_DL_S2Tvar_3mon

 temp13 = mask(temp13, temp13.gt.0, True)

 T1 = tsa_DL_S2Tvar_3mon/temp13

 T2 = 0.02564
 temp12 = temp11
 temp12 = (/T2/)

 T3 = temp11/temp13

 tsa_DL_S2Tvar_3mon_R2 = (/(T1 - temp12)/(T3)/)   ; (/---/)

 ;exit

; Now let's start plotting the data


 tsa_DL_S2Tvar_3mon_sqrt = tsa_DL_S2Tvar_3mon
 tsa_DL_S2Tvar_3mon_sqrt = (/tsa_DL_S2Tvar_3mon^0.5/)


;******************************************************************

;*****Reddedned ENSO model skill******************

 fre = addfile("/glade/work/kumar34/hydroclimate_ncc/NorthAmerica_obs_enso_h2osoi_regression_memoryENSO_5yrout_smnodt_noocdt_final_era5_1950_2021.nc", "r")

  ;s2n_re = fre->sig2noise(:, :, :)
  s2n_re = fre->actskill(:, :, :)
  s2n_sq = s2n_re
  s2n_sq = (/s2n_re^2.0/)
  s2t_sq = s2n_sq
  onemat = s2t_sq
  onemat = 1.0
  s2t_sq = s2n_sq / (s2n_sq + onemat)

  s2t_re = s2t_sq
  s2t_re = (/s2t_sq^0.5/) 


  latS = 15.
  latN = 72.
  lonW = 192.
  lonE = 308.

 tsa_DL_S2Tvar_3mon_sqrt(0, {latS:latN}, {lonW:lonE}) = (/s2n_re(2, :, :)/)   

 tsa_DL_S2Tvar_3mon_sqrt(1, {latS:latN}, {lonW:lonE}) = (/s2n_re(1, :, :)/)
 tsa_DL_S2Tvar_3mon_sqrt(2, {latS:latN}, {lonW:lonE}) = (/s2n_re(0, :, :)/)


;*****************************************************************************

  wks = gsn_open_wks("x11", "obs_enso_memory_model_act_skill_1950_2021_era5_lake_mead_rev3")

 ; gsn_define_colormap(wks, "sunshine_9lev")

 gsn_define_colormap(wks,"CBR_wet")

  

 res                     = True
 res@gsnDraw             = False
 res@gsnFrame            = False
 res@gsnSpreadColors     = True
 res@lbLabelBarOn        = False            ; turn off individual label bars
 res@mpFillOn            = False            ; no need
 res@gsnAddCyclic        = True
 res@mpLimitMode          = "LatLon"
 res@mpMinLatF            = 15.0
 res@mpMaxLatF            = 75.0
 ;res@mpMinLonF            = elon1(0)
 res@mpMinLonF            = -165.0
 res@mpMaxLonF            =  -55.0
 ;res@cnLabelBarEndStyle   = "ExcludeOuterBoxes"
; res@mpProjection = "Robinson"
  res@gsnTickMarksOn = True
res@gsnMaximize = False

  latS = res@mpMinLatF
 latN = res@mpMaxLatF
 lonE = res@mpMinLonF + 360.0
 lonW = res@mpMaxLonF  + 360.0

 ;res@mpLimitMode          = "LatLon"
 ;res@mpMinLatF            = 10
 ;res@mpMaxLatF            = 75
 ;res@mpMinLonF            = -170
 ;res@mpMaxLonF            =  -50


   res@mpLimitMode       = "LatLon"
   res@mpMinLatF         = min(s2n&lat)
   res@mpMaxLatF         = max(s2n&lat)
   res@mpMinLonF         = min(s2n&lon)
   res@mpMaxLonF         = max(s2n&lon)
   res@mpCenterLonF      = dim_median(s2n&lon) ;(min(lon)+max(lon))/2
   res@mpOutlineBoundarySets = "Geophysical"
   res@mpGeophysicalLineColor = "black"

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
  
  print("printing the skill for year 1") 

  do r = 0, 2
     print(wgt_areaave(tsa_DL_S2Tvar_3mon_sqrt(0, {latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0))

  end do


 ;tsa_DL_S2Tvar_3mon( i, :, :)



 res@mpOutlineBoundarySets = "Geophysical"
 ;res@mpOutlineBoundarySets = "National"
 res@mpGeophysicalLineColor = "black"

 ;res@pmTickMarkDisplayMode = "Always"

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
 ;res@cnMissingValFillColor = "gray"

 ;*******PLEASE CHECK THIS SECTION-A *************************
 ;*******SECTION A*******************************************
 ; ***first plot by commenting all the res in section A****
 ; *** the decide which maximum and minimum value you want to assign *****
 ; **** accordingly decide the min, max and stide for manual cn level selction *****
 ;res@cnLevelSelectionMode = "ManualLevels"
 ;res@cnMinLevelValF =  0.03
 ;res@cnMaxLevelValF =  0.27
 ;res@cnLevelSpacingF = 0.03
 
  ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  ;res@cnLevels    = (/ 0.03, .05, 0.1, 0.15 , 0.2 , 0.25, 0.3 , 0.35, 0.4, 0.45, 0.5 /)   ; set levels

  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = (/.1,.2,.3,.4,.5,.6, 0.7/)




 
 ;******************************************

;res@cnLevelSelectionMode = "ExplicitLevels"
 ;res@cnLevels = (/-40.0, -35.0, -30.0, -25.0, -21.0, -18.0, -15.0, -12.0, -10.0, -9.0, -8.0, -7.0, -6.0, -4.0, -3.0, -2.0, -1.0, 0.0/)*1
 ;res@cnLevels = (/1.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0/)*1
 ;res@cnExplicitLabelBarLabelsOn = True

  plots = new(10, graphic)


  ;delete(tsa_DL_S2Tvar_4S@long_name)
  ;delete(tsa_LC_S2Tvar_4S@long_name)
  ;delete(tsa_S2S_var_4S@long_name)

  ;delete(tsa_DL_S2Tvar_4S@units)
  ;delete(tsa_LC_S2Tvar_4S@units)
  ;delete(tsa_S2S_var_4S@units)

  ;delete(tsa_DL_S2Tvar_sig_4S@long_name)
  ;delete(tsa_LC_S2Tvar_sig_4S@long_name)
  ;delete(tsa_S2S_var_sig_4S@long_name)

  ;delete(tsa_DL_S2Tvar_sig_4S@units)
  ;delete(tsa_LC_S2Tvar_sig_4S@units)
  ;delete(tsa_S2S_var_sig_4S@units)


   res@tmYLLabelsOn               = False            ; turn off lat labels
   res@tmXBLabelsOn               = False            ; turn off lon labels

   res@tmYROn                     = False
   res@tmYLOn                     = False
   res@tmXTOn                     = False
   res@tmXBOn                     = False



;   tsa_DL_S2Tvar_3mon_sqrt = tsa_DL_S2Tvar_3mon
;   tsa_DL_S2Tvar_3mon_sqrt = (/tsa_DL_S2Tvar_3mon^0.5/)


  FirstID  = addfile("/glade/work/yanand/USDA_FACT/SMYLE/Corre_ERA5_SMYLE_10CM_SM_MAY_NOV_6_18_lead_1970_2019_NorthAmerica.nc","r")
  corre_sm = FirstID->corre_sm(:,:,0) 
  printVarSummary(corre_sm) 

   
 var_new1 = tsa_DL_S2Tvar_3mon_sqrt(0, :, :)

 lat_new = corre_sm&lat
 lon_new = corre_sm&lon

 var_new1({lat_new(0):lat_new(59)}, {lon_new(0):lon_new(92)}) = (/corre_sm/)

 tsa_DL_S2Tvar_3mon_sqrt(3, :, :) = var_new1

 do i = 0, 2

    plots(i) = gsn_csm_contour_map(wks, tsa_DL_S2Tvar_3mon_sqrt( i, :, :), res)

 end do
  
 plots(3) = gsn_csm_contour_map(wks, var_new1, res)

  print("SMYLE skill")

  do r = 0, 2
     print(wgt_areaave(var_new1({latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0))
  end do




 
 do i = 4, 9

    plots(i) = gsn_csm_contour_map(wks, tsa_DL_S2Tvar_3mon_sqrt( i, :, :), res)

 end do
 

  var_1d = ndtooned(tsa_DL_S2Tvar_3mon_sqrt(0, {a_lat:b_lat}, {c_lon:d_lon}))

  P102 = new((/dimsizes(var_1d), 2/), "double", -9999.0)

  P102(:, 0) = (/var_1d/)

  opt102 = True
  opt102@fout = "DPLE_pot_sqrt_NA_skill.txt"
  fmtx = "2f15.8"
  write_matrix(P102, fmtx, opt102)
 

  gsres                 = True
  gsres@gsMarkerIndex  = 15 
  gsres@gsMarkerColor = "red"
  gsres@gsMarkerThicknessF = 5.0
  gsres@gsMarkerSizeF = 0.018  
  ;gres@gsFillColor     = "yellow"
  ;gres@tfPolyDrawOrder = "Draw"
  dum20 = gsn_add_polymarker(wks,plots(1),-114.73, 36.016,gsres) 
 



   nps = 5

   nregs = 3

      bxres = True
    bxres@gsLineThicknessF = 4.
    bxres@gsLineColor = "black"



   ;delete(boxes)
   boxes = new((/nps*2,nregs/),graphic)

  do r = 0,2 ;---------ATTENTION
   do pm = 0,5
     boxes(pm,r) = gsn_add_polyline(wks,plots(pm),locoords(r,:),lacoords(r,:),bxres)
   end do
 end do





; Comment everything above this
   res2 = True
   res2@gsnDraw            = False             ; don't draw yet
   res2@gsnFrame           = False             ; don't advance fra

  ;res2@gsnAddCyclic = False
   res2@gsnMaximize = False

   res2@cnInfoLabelOn       = False       ; turn off info label
   res2@cnLinesOn           = False       ; do not draw contour lines
   res2@cnLineLabelsOn      = False       ; do not draw contour labels
   res2@cnFillScaleF        = 1.0         ; add extra density
   ;res2@cnFillMode = "RasterFill"
   res2@cnFillDotSizeF       = 0.0035
   ;res2@cnLineThicknessF       = 2.0

   ;res2@cnFillColor = "navyblue"


   res2@cnLevelSelectionMode = "ManualLevels"
   res2@cnMinLevelValF =  0.0
   res2@cnMaxLevelValF =  3.0
   res2@cnLevelSpacingF =  0.3
   res2@cnFillScaleF        = 1.5        ; add extra density


   opt = True
   opt@gsnShadeFillType = "pattern"      ; pattern fill
   opt@gsnShadeHigh = 17       ; use pattern #2
   ;opt@gsnShadeHigh = 9
   opt@gsnShadeLineThicknessF = 3.0
   ;opt@gsnShadeFillScales = 0.6
   ;opt@gsnFillDotSizeF       = 0.055  
 
   plots3 = plots

  tsa_DL_S2Tvar_sig_3mon = where(tsa_DL_S2Tvar_3mon_sqrt.ge.0.25, 2.5, 0.05)  

  var_new2 = var_new1
  var_new2 = (/where(var_new1.ge.0.25, 2.5, 0.05)/)  
  
  tsa_DL_S2Tvar_sig_3mon(3,:,:) = var_new2
  
  do i = 0, 2

    plots3(i) = gsn_csm_contour(wks, tsa_DL_S2Tvar_sig_3mon(i, :, :), res2)
    plots3(i) = gsn_contour_shade(plots3(i), -9999.0, 2.4, opt)  ; shade all areas greater than 95 CI
    overlay(plots(i), plots3(i))

  end do

    plots3(3) = gsn_csm_contour(wks, var_new2, res2)
    plots3(3) = gsn_contour_shade(plots3(3), -9999.0, 2.4, opt)  ; shade all areas greater than 95 CI
    overlay(plots(3), plots3(3))
	
  do i = 4, 9

    plots3(i) = gsn_csm_contour(wks, tsa_DL_S2Tvar_sig_3mon(i, :, :), res2)
    plots3(i) = gsn_contour_shade(plots3(i), -9999.0, 2.4, opt)  ; shade all areas greater than 95 CI
    overlay(plots(i), plots3(i))

  end do
  
  P101 = new((/4, 3/), "double", -9999.0)

  


delete(temp13)

do iii = 0, 3
  do jjj = 0, 2
   Marea = area
   Marea = (/mask(area, lm.eq.1, True)/)
   total_area = sum(Marea({latS(jjj):latN(jjj)}, {lonE(jjj):lonW(jjj)}))
   
   temp133 = tsa_DL_S2Tvar_3mon_sqrt(iii, :, :)
   temp13  = temp133
   temp13 = (/mask(temp133, lm.eq.1, True)/)
   temp_area = temp133
   print(typeof(temp_area))
   print(typeof(Marea*temp13))
   temp_area = doubletofloat((/Marea*temp13/))
   ;temp_area = mask(temp_area, temp13.ge.1.8, True)
   P101(iii, jjj) = (sum(temp_area({latS(jjj):latN(jjj)}, {lonE(jjj):lonW(jjj)})))*100.0/total_area



   ;plots3(i) = gsn_csm_contour(wks, temp13, res2)
   ;plots3(i) = gsn_contour_shade(plots3(i), -9999.0,1.8, opt)  ; shade all areas greater than 95 CI
   ;overlay(plots(i), plots3(i))
   delete(temp_area)
 end do
end do

 opt101 = True
 opt101@fout = "sm_pred_dple_10_year.txt"
 fmtx = "3f15.8"
 write_matrix(P101, fmtx, opt101)
 

;   plots3(1) = gsn_csm_contour(wks, tsa_DL_S2Tvar_sig_4S(s_id, 0, :, :), res2)
 ;  plots3(1) = gsn_contour_shade(plots3(1), -9999.0, 1.8, opt)  ; shade all areas greater than 95 CI
 ;  overlay(plots(1), plots3(1))

 ;  plots3(2) = gsn_csm_contour(wks, tsa_S2S_var_sig_4S(s_id, 0, :, :), res2)
 ;  plots3(2) = gsn_contour_shade(plots3(2), -9999.0, 1.8, opt)  ; shade all areas greater than 95 CI
 ;  overlay(plots(2), plots3(2))

  ;mon = ispan(1,120,1)


  resP                  = True                   ; modify the panel plot
  resP@gsnFrame         = False                  ; don't advance the frame, so we can use gsn_text_ndc
  resP@gsnPanelLabelBar = True                   ; add common colorbar
  resP@lbLabelFontHeightF = 0.015                ; set font height of Label Bar labels
  resP@gsnPanelBottom   = 0.1                   ; shrink panel plot by setting bottom edge of plot
  resP@gsnPanelTop      = 0.9                    ; shrink panel plot by setting top edge of plot
  resP@gsnPanelYWhiteSpacePercent = 1.           ; increase spacing along Y-axis between panel plots
  resP@gsnPanelMainString = "Soil moisture (100cm)"     ; set panel title
    
  ;resP@gsnPanelFigureStrings = (/"1-3", "4-6", "7-9", "10-12", "13-15", "16-18", "19-21", "22-24", "25-27", "28-30", "31-33", "34-36", "37-39", "40-42", "43-45", "46-48", "49-51", "52-54", "55-57", "58-60", "61-63", "64-66", "67-69", "70-72", "73-75", "76-78", "79-81", "82-84","85-87", "88-90", "91-93", "94-96", "97-99","100-102", "103-105", "106-108", "109-111", "112-114", "115-117", "118-120" /)
 
  ;resP@gsnPanelFigureStrings = (/"Yr1","Yr2","Yr3","Yr4","Yr5","Yr6","Yr7","Yr8","Yr9","Yr10"/) 


  resP@gsnPanelFigureStrings = (/"[a] SMYLE", "[b] Reddened ENSO", "[c] ENSO", "[d] Memory"/)

 resP@lbBoxMinorExtentF = 0.8
  resP@vpWidthF = 0.5
  resP@pmLabelBarOrthogonalPosF = -0.03
 ; resP@lbLabelStride         = 2
  resP@gsnPanelFigureStringsFontHeightF = 0.015
;  resP@gsnPanelFigureStringsFontColor   = "gray8"
;  resP@gsnPanelFigureStringsBackgroundFillColor = -1
  resP@gsnPanelFigureStringsPerimOn = "False"
  resP@amJust   = "BottomLeft"
  resP@gsnPanelRowSpec = True
  resP@gsnPanelCenter = False
  res@gsnPanelLeft = 0.05 
  res@gsnPanelRight = 0.95
  plots2 = new(4, graphic)
  plots2(1:3) = plots(0:2) 
  plots2(0) = plots(3) 
  gsn_panel(wks,plots2,(/2,2/),resP)               ; now draw as one plot



  gres = True
  gres@gsMarkerIndex = 4
  gres@gsMarkerColor          = "blue"       ;-- marker color

  ;dum20 = gsn_add_polymarker(wks,plots(0),360.0-114.41,36.14,gres)

  ;txres               = True
  ;txres@txFontHeightF = 0.015
  ;gsn_text_ndc(wks,"Figure 16: A smaller panel plot",0.5,0.16,txres)
  frame(wks)








 
