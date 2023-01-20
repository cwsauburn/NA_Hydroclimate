   url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"



   varname = "H2OSOI"
   vartitle = "soil moisture"
   ocnvar = "ENSO"

   ;dirlnd = "/Projects/SERDP_regimeshifts/data/CESM-LE/lnd/postprocessed/"
   dirlnd = url+"CESM-LE/lnd/postprocessed/"
   filelnd = dirlnd+"NorthAmerica_"+varname+"_CESM-LE_1940-2089_30yr_MemoryENSOregressions.nc"
   fl = addfile(filelnd,"r")
   s2n = fl->sig2noise ;signal-to-noise of 3 regression models
   coefs = fl->coefs_ens ; all coefficients of the models
   askill = fl->actskill ; actual (cross-validated) skill
   pskill = fl->predskill ; rho_infinity

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

   plottype = "png"
   plottype@wkWidth  = 4200
   plottype@wkHeight = 4200


   dirout = "/glade/work/kumar34/hydroclimate_ncc/"

;   dirout = "/Users/mnewman/DewesetalNCC/"
;   dirout = "/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/plots/"

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


   ; Plot regression coefficients **********************
   
   plotname = dirout+"CESM-LE_"+varname+"_MemoryENSO_Coefs_30yr_maps"
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

   tres@txAngleF      = 90.
   gsn_text_ndc(wks,"Memory",.025,.85,tres)
   gsn_text_ndc(wks,"ENSO",.02,.65,tres)
   gsn_text_ndc(wks,"Memory+ENSO",.02,.3,tres)

;   drawNDCGrid(wks)

   frame(wks)
   system("convert -trim "+plotname+"."+plottype+" "+plotname+"."+plottype)
   
   delete([/wks,res,txres,adres,tres,pres,plot,tx,ad/])

   ; Plot actual skill **********************
   
   plotname = dirout+"CESM-LE_"+varname+"_MemoryENSO_ActSkill_30yr_maps"
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

   plot = new(6,graphic)
   plot(0) = gsn_csm_contour_map(wks,askill(0,p2p(0),:,:),res)
   plot(1) = gsn_csm_contour_map(wks,askill(0,p2p(1),:,:),res)
   plot(2) = gsn_csm_contour_map(wks,askill(1,p2p(0),:,:),res)
   plot(3) = gsn_csm_contour_map(wks,askill(1,p2p(1),:,:),res)
   plot(4) = gsn_csm_contour_map(wks,askill(2,p2p(0),:,:),res)
   plot(5) = gsn_csm_contour_map(wks,askill(2,p2p(1),:,:),res)

   tres = True
   tres@txFontHeightF = 0.015
   gsn_text_ndc(wks,yrstr(p2p(0))+"-"+yrend(p2p(0)),.33,.885,tres)
   gsn_text_ndc(wks,yrstr(p2p(1))+"-"+yrend(p2p(1)),.79,.885,tres)

   tres@txAngleF      = 90.
   gsn_text_ndc(wks,"Memory",.05,.765,tres)
   gsn_text_ndc(wks,"ENSO",.05,.54,tres)
   gsn_text_ndc(wks,"Memory+ENSO",.05,.31,tres)
   
;   drawNDCGrid(wks)

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
   
   plotname = dirout+"CESM-LE_"+varname+"_MemoryENSO_PredSkill_30yr_maps"
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

   plot = new(6,graphic)
   plot(0) = gsn_csm_contour_map(wks,pskill(0,p2p(0),:,:),res)
   plot(1) = gsn_csm_contour_map(wks,pskill(0,p2p(1),:,:),res)
   plot(2) = gsn_csm_contour_map(wks,pskill(1,p2p(0),:,:),res)
   plot(3) = gsn_csm_contour_map(wks,pskill(1,p2p(1),:,:),res)
   plot(4) = gsn_csm_contour_map(wks,pskill(2,p2p(0),:,:),res)
   plot(5) = gsn_csm_contour_map(wks,pskill(2,p2p(1),:,:),res)


   P1 = new((/3, 2/), "float", -9999.0)

   r = 0


  do i = 0, 2
     do j = 0, 1
        temp11 = pskill(i, p2p(j), {latS(r):latN(r)},{lonW(r):lonE(r)})
        P1(i, j) = dim_avg_n_Wrap(dim_avg_n_Wrap(temp11, 0), 0) 
      end do
  end do
 

  opt1 = True
  opt1@fout = "pred_skill_southwest_fig3.txt"
  fmtx = "2f15.8"
  write_matrix(P1, fmtx, opt1)
  

    bxres = True
    bxres@gsLineThicknessF = 4.
    bxres@gsLineColor = "black"

   nps = 5

   boxes = new((/nps*2,nregs/),graphic)

  do r = 0,0
   do pm = 0,5
     boxes(pm,r) = gsn_add_polyline(wks,plot(pm),locoords(r,:),lacoords(r,:),bxres)
   end do
 end do





   tres = True
   tres@txFontHeightF = 0.015
   gsn_text_ndc(wks,yrstr(p2p(0))+"-"+yrend(p2p(0)),.33,.885,tres)
   gsn_text_ndc(wks,yrstr(p2p(1))+"-"+yrend(p2p(1)),.79,.885,tres)

   tres@txAngleF      = 90.
   gsn_text_ndc(wks,"Memory",.05,.765,tres)
   gsn_text_ndc(wks,"ENSO",.05,.54,tres)
   gsn_text_ndc(wks,"Memory+ENSO",.05,.31,tres)
   
;   drawNDCGrid(wks)

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

  
