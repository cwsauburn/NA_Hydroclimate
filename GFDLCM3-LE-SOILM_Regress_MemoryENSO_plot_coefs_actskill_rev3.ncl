
   url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"


   varname = "SOILM"
   vartitle = "soil moisture"
   ocnvar = "ENSO"

   dirlnd = url+"/GFDL-CM3-LE/SOILM/"
   filelnd = dirlnd+"NorthAmerica_"+varname+"_GFDL-CM3-LE_1940-2089_30yr_MemoryENSOregressions.nc"
   fl = addfile(filelnd,"r")
   s2n = fl->sig2noise
   coefs = fl->coefs_ens
   askill = fl->actskill
   pskill = fl->predskill

   nyrp = 30
   yrstr = (/1940,1970,2000,2030,2060/)
   yrend = yrstr + nyrp - 1
   nps = dimsizes(yrstr)

   printMinMax(s2n,True)
   printMinMax(coefs,True)
   printMinMax(askill,True)
   printMinMax(pskill,True)

   ; periods to plot
   ; second and last 30yr periods
   p2p = (/1,4/)
   
;***************************************************************************************

   plottype = "x11"
   ;plottype@wkWidth  = 4200
   ;plottype@wkHeight = 4200

   ;dirout = "/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/plots/"

  dirout = "/glade/work/kumar34/hydroclimate_ncc/"


 region = (/"Southwest","Southeast", "Canadian Plains"/)
 nregs = dimsizes(region)

 ;latS = (/30.,25.,28.,43.,53./)
 ;latN = (/40.,35.,38.,49.,61./)
 ;lonW = (/239.,265.,255.,254.,242./)
 ;lonE = (/252.,282.,263.,264.,257./)


 latS = (/30.,25.,53./)
 latN = (/40.,35.,61./)
 lonW = (/239.,265.,242./)
 lonE = (/252.,282.,257./)



  lacoords = new((/nregs,5/),"float")
 locoords = lacoords
 do r = 0,nregs-1
   lacoords(r,:) = (/latS(r),latN(r),latN(r),latS(r),latS(r)/)
   locoords(r,:) = (/lonW(r),lonW(r),lonE(r),lonE(r),lonW(r)/)
 end do





   ; Plot regression coefficients **********************
   
   plotname = dirout+"GFDLCM3-LE_"+varname+"_MemoryENSO_Coefs_30yr_maps"
   system("rm -f "+plotname+"."+plottype)
   wks  = gsn_open_wks (plottype,plotname)

   gsn_define_colormap(wks,"BlueWhiteOrangeRed")

   res                   = True
   res@gsnDraw           = False
   res@gsnFrame          = False
   res@gsnAddCyclic      = False
   res@gsnSpreadColors   = True

   res@mpLimitMode       = "LatLon"
   res@mpMinLatF         = min(s2n&lat)
   res@mpMaxLatF         = max(s2n&lat)
   res@mpMinLonF         = min(s2n&lon)
   res@mpMaxLonF         = max(s2n&lon)
   res@mpCenterLonF      = dim_median(s2n&lon) ;(min(lon)+max(lon))/2
   res@mpOutlineBoundarySets = "Geophysical"
   res@mpGeophysicalLineColor = "black"
   
   res@pmTickMarkDisplayMode = "Always"

   res@cnInfoLabelOn     = False
   res@cnFillOn          = True
   res@cnLinesOn         = False
   res@cnFillMode        = "RasterFill"
   res@cnRasterSmoothingOn = True
   res@cnLineLabelsOn    = False

   res@cnLevelSelectionMode = "ExplicitLevels"
   res@cnLevels = (/-.7,-.6,-.5,-.4,-.3,-.2,.2,.3,.4,.5,.6,.7/)

   res@lbLabelBarOn      = False
;   res@vpXF              = .45

   txres = True
   txres@txFontHeightF = 0.035

   adres = True
   adres@amJust           = "BottomLeft"
   adres@amParallelPosF   = -0.36    ; 0.5 is the right edge of the plot.
   adres@amOrthogonalPosF = 0.26   ; -0.5 is the top edge of the plot.

   plot = new(8,graphic)
   tx = new(4,graphic)
   ad = new(4,graphic)
 

     res@tmYLLabelsOn               = False            ; turn off lat labels
    res@tmXBLabelsOn               = False            ; turn off lon labels

    res@tmYROn                     = False
    res@tmYLOn                     = False
    res@tmXTOn                     = False
    res@tmXBOn                     = False




   plot(0) = gsn_csm_contour_map(wks,coefs(0,p2p(0),0,:,:),res)
   plot(1) = gsn_csm_contour_map(wks,coefs(0,p2p(1),0,:,:),res)
   tx(0) = gsn_create_text(wks,"a",txres)
   ad(0) = gsn_add_annotation(plot(0),tx(0),adres)

   plot(2) = gsn_csm_contour_map(wks,coefs(1,p2p(0),1,:,:),res)
   plot(3) = gsn_csm_contour_map(wks,coefs(1,p2p(1),1,:,:),res)
   tx(1) = gsn_create_text(wks,"b",txres)
   ad(1) = gsn_add_annotation(plot(2),tx(1),adres)

   plot(4) = gsn_csm_contour_map(wks,coefs(2,p2p(0),0,:,:),res)
   plot(5) = gsn_csm_contour_map(wks,coefs(2,p2p(1),0,:,:),res)
   tx(2) = gsn_create_text(wks,"~F33~a",txres)
   ad(2) = gsn_add_annotation(plot(4),tx(2),adres)

   plot(6) = gsn_csm_contour_map(wks,coefs(2,p2p(0),1,:,:),res)
   plot(7) = gsn_csm_contour_map(wks,coefs(2,p2p(1),1,:,:),res)
   tx(3) = gsn_create_text(wks,"~F33~b",txres)
   ad(3) = gsn_add_annotation(plot(6),tx(3),adres)



   P101 = new((/3, 8/), "double", -9999.0) 
   nReg = 3  

   latwt = latRegWgt(coefs&lat, "float", 0) 

  latwt! = "lat"
  latwt&lat = coefs&lat

   print(latwt)
  

   ; exit


   do r = 0, nReg-1
      P101(r, 0) = wgt_areaave(coefs(2, p2p(0),0, {latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)
      P101(r, 1) = wgt_areaave(coefs(2, p2p(0),1, {latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)
      
      P101(r, 4) = wgt_areaave(coefs(2, p2p(1),0, {latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)
      P101(r, 5) = wgt_areaave(coefs(2, p2p(1),1, {latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)

   end do








     bxres = True
   bxres@gsLineThicknessF = 4.
    bxres@gsLineColor = "black"

   nps = 8

   boxes = new((/nps*2,nregs/),graphic)

  do r = 0,4
   do pm = 0,7
     boxes(pm,r) = gsn_add_polyline(wks,plot(pm),locoords(r,:),lacoords(r,:),bxres)
   end do
 end do







   pres = True
   pres@gsnFrame = False
;   pres@gsnPanelRight = 0.995
;   pres@gsnPanelLeft  = 0.064
   pres@gsnPanelTop = .97
;   pres@gsnPanelDebug = True
;   pres@gsnPanelYF = (/0.775,0.775,0.526,0.526/)
;   pres@gsnPanelXF = (/.07,.51,.07,.51/)
   pres@gsnPanelLabelBar    = True
   pres@pmLabelBarWidthF    = 0.6
   pres@pmLabelBarHeightF   = 0.05
;   pres@lbLabelFontHeightF  = 0.01
   gsn_panel(wks,plot,(/4,2/),pres)
  
   tres = True
   tres@txFontHeightF = 0.015
   gsn_text_ndc(wks,yrstr(p2p(0))+"-"+yrend(p2p(0)),.3,.985,tres)
   gsn_text_ndc(wks,yrstr(p2p(1))+"-"+yrend(p2p(1)),.75,.985,tres)

;   tres@txAngleF      = 90.
;   gsn_text_ndc(wks,"Memory",.025,.85,tres)
;   gsn_text_ndc(wks,"ENSO",.02,.65,tres)
;   gsn_text_ndc(wks,"Memory+ENSO",.02,.3,tres)

;   drawNDCGrid(wks)

   frame(wks)
   system("convert -trim "+plotname+"."+plottype+" "+plotname+"."+plottype)
   
   delete([/wks,res,tres,pres,plot/])

   ; Plot actual skill **********************
   
   plotname = dirout+"GFDLCM3-LE_"+varname+"_MemoryENSO_ActSkill_30yr_maps"
   system("rm -f "+plotname+"."+plottype)
   wks  = gsn_open_wks (plottype,plotname)

   gsn_define_colormap(wks,"CBR_wet")
;   gsn_reverse_colormap(wks)

   res                   = True
   res@gsnDraw           = False
   res@gsnFrame          = False
   res@gsnAddCyclic      = False
   res@gsnSpreadColors   = True

   res@mpLimitMode       = "LatLon"
   res@mpMinLatF         = min(s2n&lat)
   res@mpMaxLatF         = max(s2n&lat)
   res@mpMinLonF         = min(s2n&lon)
   res@mpMaxLonF         = max(s2n&lon)
   res@mpCenterLonF      = dim_median(s2n&lon) ;(min(lon)+max(lon))/2
   res@mpOutlineBoundarySets = "Geophysical"
   res@mpGeophysicalLineColor = "black"
   
   res@pmTickMarkDisplayMode = "Always"

   res@cnInfoLabelOn     = False
   res@cnFillOn          = True
   res@cnLinesOn         = False
   res@cnFillMode        = "RasterFill"
   res@cnRasterSmoothingOn = True
   res@cnLineLabelsOn    = False

   res@cnLevelSelectionMode = "ExplicitLevels"
   res@cnLevels = (/.1,.2,.3,.4,.5,.6,.7/)

   res@lbLabelBarOn      = False

      res@tmYLLabelsOn               = False            ; turn off lat labels
    res@tmXBLabelsOn               = False            ; turn off lon labels

    res@tmYROn                     = False
    res@tmYLOn                     = False
    res@tmXTOn                     = False
    res@tmXBOn                     = False





   plot = new(6,graphic)
   plot(0) = gsn_csm_contour_map(wks,askill(0,p2p(0),:,:),res)
   plot(1) = gsn_csm_contour_map(wks,askill(0,p2p(1),:,:),res)
   plot(2) = gsn_csm_contour_map(wks,askill(1,p2p(0),:,:),res)
   plot(3) = gsn_csm_contour_map(wks,askill(1,p2p(1),:,:),res)
   plot(4) = gsn_csm_contour_map(wks,askill(2,p2p(0),:,:),res)
   plot(5) = gsn_csm_contour_map(wks,askill(2,p2p(1),:,:),res)



   do r = 0, nReg-1
      P101(r, 2) = wgt_areaave(askill(2, p2p(0),{latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)
      P101(r, 6) = wgt_areaave(askill(2, p2p(1),{latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)
   end do




   tres = True
   tres@txFontHeightF = 0.015
   gsn_text_ndc(wks,yrstr(p2p(0))+"-"+yrend(p2p(0)),.33,.885,tres)
   gsn_text_ndc(wks,yrstr(p2p(1))+"-"+yrend(p2p(1)),.79,.885,tres)

   ;tres@txAngleF      = 90.
   ;gsn_text_ndc(wks,"Memory",.05,.765,tres)
   ;gsn_text_ndc(wks,"ENSO",.05,.54,tres)
   ;gsn_text_ndc(wks,"Memory+ENSO",.05,.31,tres)
   
;   drawNDCGrid(wks)

    bxres = True
   bxres@gsLineThicknessF = 4.
    bxres@gsLineColor = "black"

   nps = 8

   boxes = new((/nps*2,nregs/),graphic)

  do r = 0,4
   do pm = 0,5
     boxes(pm,r) = gsn_add_polyline(wks,plot(pm),locoords(r,:),lacoords(r,:),bxres)
   end do
 end do






   pres = True
   pres@gsnPanelRight = 0.995
   pres@gsnPanelLeft  = 0.064
;   pres@gsnPanelDebug = True
;   pres@gsnPanelYF = (/0.775,0.775,0.526,0.526/)
;   pres@gsnPanelXF = (/.07,.51,.07,.51/)
   pres@gsnPanelLabelBar    = True
   pres@pmLabelBarWidthF    = 0.6
   pres@pmLabelBarHeightF   = 0.05
;   pres@lbLabelFontHeightF  = 0.01
   gsn_panel(wks,plot,(/3,2/),pres)
  
   system("convert -trim "+plotname+"."+plottype+" "+plotname+"."+plottype)
   
   delete([/wks,res,tres,pres,plot/])

   ; Plot predicted skill **********************
   
   plotname = dirout+"GFDLCM3-LE_"+varname+"_MemoryENSO_PredSkill_30yr_maps"
   system("rm -f "+plotname+"."+plottype)
   wks  = gsn_open_wks (plottype,plotname)

   gsn_define_colormap(wks,"CBR_wet")
;   gsn_reverse_colormap(wks)

   res                   = True
   res@gsnDraw           = False
   res@gsnFrame          = False
   res@gsnAddCyclic      = False
   res@gsnSpreadColors   = True

   res@mpLimitMode       = "LatLon"
   res@mpMinLatF         = min(s2n&lat)
   res@mpMaxLatF         = max(s2n&lat)
   res@mpMinLonF         = min(s2n&lon)
   res@mpMaxLonF         = max(s2n&lon)
   res@mpCenterLonF      = dim_median(s2n&lon) ;(min(lon)+max(lon))/2
   res@mpOutlineBoundarySets = "Geophysical"
   res@mpGeophysicalLineColor = "black"
   
   res@pmTickMarkDisplayMode = "Always"

   res@cnInfoLabelOn     = False
   res@cnFillOn          = True
   res@cnLinesOn         = False
   res@cnFillMode        = "RasterFill"
   res@cnRasterSmoothingOn = True
   res@cnLineLabelsOn    = False

   res@cnLevelSelectionMode = "ExplicitLevels"
   res@cnLevels = (/.1,.2,.3,.4,.5,.6,.7/)

   res@lbLabelBarOn      = False

     res@tmYLLabelsOn               = False            ; turn off lat labels
    res@tmXBLabelsOn               = False            ; turn off lon labels

    res@tmYROn                     = False
    res@tmYLOn                     = False
    res@tmXTOn                     = False
    res@tmXBOn                     = False




   plot = new(6,graphic)
   plot(0) = gsn_csm_contour_map(wks,pskill(0,p2p(0),:,:),res)
   plot(1) = gsn_csm_contour_map(wks,pskill(0,p2p(1),:,:),res)
   plot(2) = gsn_csm_contour_map(wks,pskill(1,p2p(0),:,:),res)
   plot(3) = gsn_csm_contour_map(wks,pskill(1,p2p(1),:,:),res)
   plot(4) = gsn_csm_contour_map(wks,pskill(2,p2p(0),:,:),res)
   plot(5) = gsn_csm_contour_map(wks,pskill(2,p2p(1),:,:),res)

   do r = 0, nReg-1
      P101(r, 3) = wgt_areaave(pskill(2, p2p(0),{latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)
      P101(r, 7) = wgt_areaave(pskill(2, p2p(1),{latS(r):latN(r)}, {lonW(r):lonE(r)}), latwt({latS(r):latN(r)}), 1.0, 0)
   end do


   fmtx = "8f15.8"
   opt101 = True
   opt101@fout = "GFDL_M_E_model_alpha_beta_act_pot_skill_1985_2075.txt"
   write_matrix(P101, fmtx, opt101)
    


   tres = True
   tres@txFontHeightF = 0.015
   gsn_text_ndc(wks,yrstr(p2p(0))+"-"+yrend(p2p(0)),.33,.885,tres)
   gsn_text_ndc(wks,yrstr(p2p(1))+"-"+yrend(p2p(1)),.79,.885,tres)

   ;tres@txAngleF      = 90.
   ;gsn_text_ndc(wks,"Memory",.05,.765,tres)
   ;gsn_text_ndc(wks,"ENSO",.05,.54,tres)
   ;gsn_text_ndc(wks,"Memory+ENSO",.05,.31,tres)
   
;   drawNDCGrid(wks)

   bxres = True
   bxres@gsLineThicknessF = 4.
    bxres@gsLineColor = "black"

   nps = 8

   boxes = new((/nps*2,nregs/),graphic)

  do r = 0,4
   do pm = 0,5
     boxes(pm,r) = gsn_add_polyline(wks,plot(pm),locoords(r,:),lacoords(r,:),bxres)
   end do
 end do








   pres = True
   pres@gsnPanelRight = 0.995
   pres@gsnPanelLeft  = 0.064
;   pres@gsnPanelDebug = True
;   pres@gsnPanelYF = (/0.775,0.775,0.526,0.526/)
;   pres@gsnPanelXF = (/.07,.51,.07,.51/)
   pres@gsnPanelLabelBar    = True
   pres@pmLabelBarWidthF    = 0.6
   pres@pmLabelBarHeightF   = 0.05
;   pres@lbLabelFontHeightF  = 0.01
   gsn_panel(wks,plot,(/3,2/),pres)
  
   system("convert -trim "+plotname+"."+plottype+" "+plotname+"."+plottype)
   
   delete([/wks,res,tres,pres,plot/])

 
