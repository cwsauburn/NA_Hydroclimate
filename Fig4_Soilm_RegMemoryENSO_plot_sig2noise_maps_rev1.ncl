
url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"


;***************************************************************************************
; plot signal-to-noise maps for Memory+ENSO linear model
; one row per period, one column for each GCM

 plottype = "x11"
 ;plottype@wkWidth  = 4200
 ;plottype@wkHeight = 4200

 ;dirout = "/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/plots/"


 dirout = "/glade/work/kumar34/hydroclimate_ncc/"

 plotname = dirout+"Soilm_RegMemoryENSO_sig2noise_maps"
 system("rm -f "+plotname+"."+plottype)
 wks  = gsn_open_wks (plottype,plotname)

;***************************************************************************************

; specify 30-yr periods
 nyrp = 30
 yrstr = (/1940,1970,2000,2030,2060/)
 yrend = yrstr + nyrp - 1
; periods = (/"1940-~C~1969","1970-~C~1999","2000-~C~2029","2030-~C~2059","2060-~C~2089"/)
 nps = dimsizes(yrstr)
 
;***************************************************************************************

; selected regions
 
 region = (/"Southwest","Southeast","South US Plains","North US Plains","Canadian Plains"/)
 nregs = dimsizes(region)

 latS = (/30.,25.,28.,43.,53./)
 latN = (/40.,35.,38.,49.,61./)
 lonW = (/239.,265.,255.,254.,242./)
 lonE = (/252.,282.,263.,264.,257./)

 lacoords = new((/nregs,5/),"float")
 locoords = lacoords
 do r = 0,nregs-1
   lacoords(r,:) = (/latS(r),latN(r),latN(r),latS(r),latS(r)/)
   locoords(r,:) = (/lonW(r),lonW(r),lonE(r),lonE(r),lonW(r)/)
 end do
 
;***************************************************************************************

; read in regression metrics from CESM-LE
; dircesm = "/Projects/SERDP_regimeshifts/data/CESM-LE/lnd/postprocessed/"

 dircesm = url+"CESM-LE/lnd/postprocessed/"
 filecesm = dircesm+"NorthAmerica_H2OSOI_CESM-LE_1940-2089_30yr_MemoryENSOregressions.nc"
 flc = addfile(filecesm,"r")

 lat_c = flc->lat
 nlatc = dimsizes(lat_c)
 lon_c = flc->lon
 nlonc = dimsizes(lon_c)

 s2n_cesm = flc->sig2noise(2,:,:,:)
 s2n_cesm = (/(s2n_cesm / (s2n_cesm + 1))/)
 
 printVarSummary(s2n_cesm)
 printMinMax(s2n_cesm,True)

;***************************************************************************************

; read in regression metrics from GFDL-CM3-LE
 ;dirgfdl = "/Projects/SERDP_regimeshifts/data/GFDL-CM3-LE/SOILM/"

 dirgfdl = url+"GFDL-CM3-LE/SOILM/"
 filegfdl = dirgfdl+"NorthAmerica_SOILM_GFDL-CM3-LE_1940-2089_30yr_MemoryENSOregressions.nc"
 flg = addfile(filegfdl,"r")

 print(flg)

 exit


 lat_g = flg->lat
 nlatg = dimsizes(lat_g)
 lon_g = flg->lon
 nlong = dimsizes(lon_g)

 s2n_gfdl = flg->sig2noise(2,:,:,:)
 s2n_gfdl = (/(s2n_gfdl / (s2n_gfdl + 1))/)
 printVarSummary(s2n_gfdl)
 printMinMax(s2n_gfdl,True)

;***************************************************************************************

; plot maps for all periods for both models

 gsn_define_colormap(wks,"cmp_haxby")
 gsn_reverse_colormap(wks)

 ; resources for CESM maps
 resc                   = True
 resc@gsnDraw           = False
 resc@gsnFrame          = False
 resc@gsnAddCyclic      = False
 resc@gsnLeftString     = ""
 resc@gsnRightString    = ""

 resc@mpLimitMode       = "LatLon"
 resc@mpMinLatF         = min(lat_c)
 resc@mpMaxLatF         = max(lat_c)
 resc@mpMinLonF         = min(lon_c)
 resc@mpMaxLonF         = max(lon_c)
 resc@mpCenterLonF      = dim_median(lon_c)
 resc@mpOutlineBoundarySets = "Geophysical"
 resc@mpGeophysicalLineColor = "black"
   
 resc@pmTickMarkDisplayMode = "Always"

 resc@cnFillOn          = True
 resc@cnFillMode        = "RasterFill"
 resc@cnRasterSmoothingOn = True
 resc@cnLinesOn         = False
 resc@cnLineLabelsOn    = False
 resc@cnInfoLabelOn     = False
 resc@cnLevelSelectionMode = "ExplicitLevels"
 resc@cnLevels = (/.1,.2,.3,.4,.5,.6,.7,.8,.9/)

 resc@lbLabelBarOn      = False

 ; resources for GFDL maps
 resg = resc
 resg@mpMinLatF         := min(lat_g)
 resg@mpMaxLatF         := max(lat_g)
 resg@mpMinLonF         := min(lon_g)
 resg@mpMaxLonF         := max(lon_g)
 resg@mpCenterLonF      := dim_median(lon_g)

 ; resources for regional box outlines
 bxres = True
 bxres@gsLineThicknessF = 4.
 bxres@gsLineColor = "black"

 ; plot!
 plot_maps = new(nps*2,graphic)
 boxes = new((/nps*2,nregs/),graphic)

 pc = (/0,2,4,6,8/)
 pg = (/1,3,5,7,9/)
 do p = 0,nps-1
   plot_maps(pc(p)) = gsn_csm_contour_map(wks,s2n_cesm(p,:,:),resc)
   plot_maps(pg(p)) = gsn_csm_contour_map(wks,s2n_gfdl(p,:,:),resg)
 end do

 do r = 0,nregs-1
   do pm = 0,nps*2-1
     boxes(pm,r) = gsn_add_polyline(wks,plot_maps(pm),locoords(r,:),lacoords(r,:),bxres)
   end do
 end do

 pres = True
 pres@gsnFrame            = False
 pres@gsnPanelLabelBar    = True
 pres@gsnPanelTop         = .98
 pres@pmLabelBarWidthF    = 0.5
 pres@pmLabelBarHeightF   = 0.05
; pres@lbLabelFontHeightF  = 0.01
; pres@gsnPanelDebug       = True
 gsn_panel(wks,plot_maps,(/nps,2/),pres)

 tres = True
 tres@txFontHeightF = 0.018
 tres@txAngleF      = 90.
 xpos = .1
 ypos = (/.9,.715,.525,.335,.15/)
 do p = 0,nps-1
   gsn_text_ndc(wks,yrstr(p)+"-"+yrend(p),xpos,ypos(p),tres)
 end do
 
 tres@txAngleF      = 0.
 gsn_text_ndc(wks,"CESM-LE",.33,.99,tres)
 gsn_text_ndc(wks,"GFDL-CM3-LE",.71,.99,tres)

; drawNDCGrid(wks)

 frame(wks)
 system("convert -trim "+plotname+"."+plottype+" "+plotname+"."+plottype)
