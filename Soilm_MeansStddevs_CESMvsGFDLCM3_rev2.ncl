
;url = "https://psl.noaa.gov/thredds/catalog/Projects/SERDP_regimeshifts/data/"

url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"


;***************************************************************************************
; plot soil moisture means and standard deviations

 plottype = "x11"
 plottype@wkWidth  = 4200
 plottype@wkHeight = 4200

 ;dirout = "/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/plots/"
 dirout = "/glade/work/kumar34/hydroclimate_ncc/"
 plotname = dirout+"SoilmMeansStdv_1970-99_30yrShifts_rev4"
 system("rm -f "+plotname+"."+plottype)
 wks  = gsn_open_wks (plottype,plotname)

;***************************************************************************************

; specify start and end dates for each 30-yr period
 nyrp = 30
 s0 = 194001
 wstep = nyrp*100
 nps = 5
 strdates = new((/nps/),"integer")
 do p = 0,nps-1
   strdates(p) = s0
   s0 = s0 + wstep
 end do
 enddates = strdates + wstep - 89
 ntimp = nyrp*12
 yrstr = (/1940,1970,2000,2030,2060/)
 yrend = yrstr + nyrp - 1
 periods = (/"1940-~C~1969","1970-~C~1999","2000-~C~2029","2030-~C~2059","2060-~C~2089"/)

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

; read in soil moisture from CESM-LE, computing means and stddevs on map and per region
 ;dircesm = "/Projects/SERDP_regimeshifts/data/CESM-LE/lnd/postprocessed/"
 dircesm = url+"CESM-LE/lnd/postprocessed/"
 filecesm = dircesm+"H2OSOI_CESM-LE_40em_192001-210012_12moRunningMeans.nc"
 flc = addfile(filecesm,"r")
 varname_cesm = "H2OSOI"

 print(flc)

 ;exit

 lat_c = flc->lat
 nlatc = dimsizes(lat_c)
 lon_c = flc->lon
 nlonc = dimsizes(lon_c)
 datec = flc->date

 ensci = ispan(1,40,1)
 nEnsc = dimsizes(ensci)

 mean_cesm = new((/nps,nlatc,nlonc/),"float",1e+36)
 stdv_cesm = mean_cesm
 mean_cesm_reg = new((/nps,nregs,nEnsc/),"float",1e+36)
 stdv_cesm_reg = mean_cesm_reg
 do p = 0,nps-1
   print_table([/yrstr(p)+"-"+yrend(p)/],"%s")
   sti = ind(datec.eq.strdates(p))
   endi = ind(datec.eq.enddates(p))
   means = flc->$varname_cesm$(sti:endi,:,:,:)
   means = (/means*1000.0/)
   mn = dim_avg_n_Wrap(means,(/0,1/))
   mean_cesm(p,:,:) = mn
   sd = dim_stddev_n_Wrap(means,(/0,1/))
   stdv_cesm(p,:,:) = sd
   do r = 0,nregs-1
     means_reg = dim_avg_n(means(:,:,{latS(r):latN(r)},{lonW(r):lonE(r)}),(/2,3/)) 
     mn_rg = dim_avg_n(means_reg,0)
     mean_cesm_reg(p,r,:) = mn_rg
     sd_rg = dim_stddev_n(means_reg,0)
     stdv_cesm_reg(p,r,:) = sd_rg
   end do 
   delete([/means,mn,sd,means_reg,mn_rg,sd_rg/])
 end do

 printVarSummary(mean_cesm)
 printMinMax(mean_cesm,True)
 printVarSummary(stdv_cesm)
 printMinMax(stdv_cesm,True)

 ; compute ensemble mean for reference period (1970-99)
 em_cesm_mn_ref = dim_avg_n(mean_cesm_reg(1,:,:),1)
 em_cesm_sd_ref = dim_avg_n(stdv_cesm_reg(1,:,:),1)

; compute shifts in mean and stdv for each ensemble member, as differences from
; the reference period ensemble mean
 cesm_mn_shfts = mean_cesm_reg
 cesm_sd_shfts = stdv_cesm_reg
 do r = 0,nregs-1
   do i = 0,nEnsc-1
     do p = 0,nps-1
       cesm_mn_shfts(p,r,i) = mean_cesm_reg(p,r,i) - em_cesm_mn_ref(r)
       cesm_sd_shfts(p,r,i) = stdv_cesm_reg(p,r,i) - em_cesm_sd_ref(r)
     end do
   end do
 end do

 em_cesm_mn_shfts = dim_avg_n(cesm_mn_shfts,2)
 em_cesm_sd_shfts = dim_avg_n(cesm_sd_shfts,2)

 em_cesm_mn_shfts_std = dim_stddev_n(cesm_mn_shfts,2)
 em_cesm_sd_shfts_std = dim_stddev_n(cesm_sd_shfts,2)


 P2 = new((/nps, nregs*4/), "float", -9999.0)

 opt2 = True
 opt2@fout = "sm_cesm_mean_std_chnage.txt"
 fmtx = "20f15.8"

 do r = 0, nregs-1
    P2(:, 4*r) = em_cesm_mn_shfts(:, r)
    P2(:, 4*r+1) = 2*em_cesm_mn_shfts_std(:, r)/(40.0^0.5)
    P2(:, 4*r+2) = em_cesm_sd_shfts(:, r)
    P2(:, 4*r+3) = 2*em_cesm_sd_shfts_std(:, r)/(40.0^0.5)
 end do


 write_matrix(P2, fmtx, opt2)





 printMinMax(cesm_mn_shfts,True)
 printMinMax(cesm_sd_shfts,True)

;***************************************************************************************

; read in soil moisture from GFDL-CM3-LE, computing means and stddevs on map and per region
 ;dirgfdl = "/Projects/SERDP_regimeshifts/data/GFDL-CM3-LE/SOILM/"
 dirgfdl = url+"GFDL-CM3-LE/SOILM/"
 filegfdl = dirgfdl+"SOILM_GFDL-CM3-LE_20em_192001-210012_12moRunningMeans.nc"
 flg = addfile(filegfdl,"r")
 varname_gfdl = "SOILM"

 lat_g = flg->lat
 nlatg = dimsizes(lat_g)
 lon_g = flg->lon
 nlong = dimsizes(lon_g)
 dateg = flg->date

 ensgi = ispan(1,20,1)
 nEnsg = dimsizes(ensgi)

 mean_gfdl = new((/nps,nlatg,nlong/),"float",1e+36)
 stdv_gfdl = mean_gfdl
 mean_gfdl_reg = new((/nps,nregs,nEnsg/),"float",1e+36)
 stdv_gfdl_reg = mean_gfdl_reg
 do p = 0,nps-1
   print_table([/yrstr(p)+"-"+yrend(p)/],"%s")
   sti = ind(dateg.eq.strdates(p))
   endi = ind(dateg.eq.enddates(p))
   means = flg->$varname_gfdl$(sti:endi,:,:,:)
   means = means * 1000.0
   mn = dim_avg_n_Wrap(means,(/0,1/))
   mean_gfdl(p,:,:) = mn
   sd = dim_stddev_n_Wrap(means,(/0,1/))
   stdv_gfdl(p,:,:) = sd
   do r = 0,nregs-1
     means_reg = dim_avg_n(means(:,:,{latS(r):latN(r)},{lonW(r):lonE(r)}),(/2,3/)) 
     mn_rg = dim_avg_n(means_reg,0)
     mean_gfdl_reg(p,r,:) = mn_rg
     sd_rg = dim_stddev_n(means_reg,0)
     stdv_gfdl_reg(p,r,:) = sd_rg
   end do 
   delete([/means,mn,sd,means_reg,mn_rg,sd_rg/])
 end do

 printVarSummary(mean_gfdl)
 printMinMax(mean_gfdl,True)
 printVarSummary(stdv_gfdl)
 printMinMax(stdv_gfdl,True)
 
 ; compute ensemble mean for reference period (1970-99)
 em_gfdl_mn_ref = dim_avg_n(mean_gfdl_reg(1,:,:),1)
 em_gfdl_sd_ref = dim_avg_n(stdv_gfdl_reg(1,:,:),1)

 ; compute shifts in mean and stdv, departing from the reference period ensemble mean
 gfdl_mn_shfts = mean_gfdl_reg
 gfdl_sd_shfts = stdv_gfdl_reg
 do r = 0,nregs-1
   do i = 0,nEnsg-1
     do p = 0,nps-1
       gfdl_mn_shfts(p,r,i) = mean_gfdl_reg(p,r,i) - em_gfdl_mn_ref(r)
       gfdl_sd_shfts(p,r,i) = stdv_gfdl_reg(p,r,i) - em_gfdl_sd_ref(r)
     end do
   end do
 end do

 em_gfdl_mn_shfts = dim_avg_n(gfdl_mn_shfts,2)
 em_gfdl_sd_shfts = dim_avg_n(gfdl_sd_shfts,2)
 
 printMinMax(gfdl_mn_shfts,True)
 printMinMax(gfdl_sd_shfts,True)


 em_gfdl_mn_shfts_std = dim_stddev_n(gfdl_mn_shfts,2)
 em_gfdl_sd_shfts_std = dim_stddev_n(gfdl_sd_shfts,2)


 P2 = new((/nps, nregs*4/), "float", -9999.0)

 opt2 = True
 opt2@fout = "sm_gfdl_mean_std_chnage.txt"
 fmtx = "20f15.8"

 do r = 0, nregs-1
    P2(:, 4*r) = em_gfdl_mn_shfts(:, r)
    P2(:, 4*r+1) = 2*em_gfdl_mn_shfts_std(:, r)/(20.0^0.5)
    P2(:, 4*r+2) = em_gfdl_sd_shfts(:, r)
    P2(:, 4*r+3) = 2*em_gfdl_sd_shfts_std(:, r)/(20.0^0.5)
 end do


 write_matrix(P2, fmtx, opt2)





;***************************************************************************************

; plot maps of mean and stdv for 1970-1999 for both models
 cmap_mn = read_colormap_file("MPL_YlGn")
 cmap_sd = read_colormap_file("MPL_gist_earth")
 ;cmap_mn = read_colormap_file("CBR_drywet")
 ;cmap_sd = read_colormap_file("CBR_coldhot")

 cmap_sd = cmap_sd(::-1,:)

 ; resources for CESM maps
 res1c                   = True
 res1c@gsnDraw           = False
 res1c@gsnFrame          = False
 res1c@gsnAddCyclic      = False
 res1c@gsnLeftString     = ""
 res1c@gsnRightString    = ""

 res1c@mpLimitMode       = "LatLon"
 res1c@mpMinLatF         = min(lat_c)
 res1c@mpMaxLatF         = max(lat_c)
 res1c@mpMinLonF         = min(lon_c)
 res1c@mpMaxLonF         = max(lon_c)
 res1c@mpCenterLonF      = dim_median(lon_c)
 res1c@mpOutlineBoundarySets = "Geophysical"
 res1c@mpGeophysicalLineColor = "black"
   
 res1c@pmTickMarkDisplayMode = "Always"

 res1c@cnFillOn          = True
 res1c@cnFillMode        = "RasterFill"
 res1c@cnRasterSmoothingOn = True
 res1c@cnLinesOn         = False
 res1c@cnLineLabelsOn    = False
 res1c@cnInfoLabelOn     = False

 res1c@lbLabelBarOn      = False

 res2c = res1c

 res1c@cnLevelSelectionMode = "ExplicitLevels"
 res1c@cnLevels = tofloat(ispan(10,40,5))*10.0
 ;res1c@cnLevels = tofloat(ispan(-15,15,3)) 
 res1c@cnFillPalette     = cmap_mn

 res2c@cnLevelSelectionMode = "ExplicitLevels"
 res2c@cnLevels = tofloat(ispan(5,30,5))
 ;res2c@cnLevels = tofloat(ispan(-50,50,10))
 res2c@cnFillPalette    = cmap_sd

 ; resources for GFDL maps
 res1g = res1c
 res1g@mpMinLatF         := min(lat_g)
 res1g@mpMaxLatF         := max(lat_g)
 res1g@mpMinLonF         := min(lon_g)
 res1g@mpMaxLonF         := max(lon_g)
 res1g@mpCenterLonF      := dim_median(lon_g)

 res1g@lbLabelBarOn      = True
 res1g@lbLabelFontHeightF = 0.022
 res1g@lbLeftMarginF     = -0.5
 res1g@lbRightMarginF    = -0.5
 
 res2g = res1g
 res2g@cnLevels := tofloat(ispan(5,30,5))
 ;res2g@cnLevels = tofloat(ispan(-50,50,10))
 res2g@cnFillPalette    := cmap_sd

 ; resources for regional box outlines
 bxres = True
 bxres@gsLineThicknessF = 4.
 bxres@gsLineColor = "black"

 ; Panel #1 ===========================================
 plot_maps = new(4,graphic)
 boxes = new((/4,nregs/),graphic)

 mean_cesm2 = mean_cesm
 stdv_cesm2 = stdv_cesm
 mean_gfdl2 = mean_gfdl
 stdv_gfdl2 = stdv_gfdl

 mean_cesm(1,:,:) = mask(mean_cesm(1,:,:), mean_cesm(1,:,:).gt.0., True)
 stdv_cesm(1,:,:) = mask(stdv_cesm(1,:,:), stdv_cesm(1,:,:).gt.0., True) 
 mean_gfdl(1,:,:) = mask(mean_gfdl(1,:,:), mean_gfdl(1,:,:).gt.0., True)
 stdv_gfdl(1,:,:) = mask(stdv_gfdl(1,:,:), stdv_gfdl(1,:,:).gt.0., True) 


 do n = 0, nps-1
    mean_cesm2(n, :, :) =(/ ((mean_cesm(n, :, :) - mean_cesm(1, :, :))*100.0)/mean_cesm(1, :, :)/)
    mean_gfdl2(n, :, :) =(/ ((mean_gfdl(n, :, :) - mean_gfdl(1, :, :))*100.0)/mean_gfdl(1, :, :)/)
   
    stdv_cesm2(n, :, :) = (/((stdv_cesm(n, :, :) - stdv_cesm(1, :, :))*100.0)/stdv_cesm(1, :, :)/)
    stdv_gfdl2(n, :, :) = (/((stdv_gfdl(n, :, :) - stdv_gfdl(1, :, :))*100.0)/stdv_gfdl(1, :, :)/)
  
 end do 



 plot_maps(0) = gsn_csm_contour_map(wks,mean_cesm(1,:,:),res1c)
 plot_maps(1) = gsn_csm_contour_map(wks,stdv_cesm(1,:,:),res2c)
 plot_maps(2) = gsn_csm_contour_map(wks,mean_gfdl(1,:,:),res1g)
 plot_maps(3) = gsn_csm_contour_map(wks,stdv_gfdl(1,:,:),res2g)


  ;plot_maps(0) = gsn_csm_contour_map(wks,mean_cesm2(4,:,:),res1c)
  ;plot_maps(1) = gsn_csm_contour_map(wks,stdv_cesm2(4,:,:),res2c)
  ;plot_maps(2) = gsn_csm_contour_map(wks,mean_gfdl2(4,:,:),res1g)
  ;plot_maps(3) = gsn_csm_contour_map(wks,stdv_gfdl2(4,:,:),res2g)




 do r = 0,nregs-1
   do pm = 0,3
     boxes(pm,r) = gsn_add_polyline(wks,plot_maps(pm),locoords(r,:),lacoords(r,:),bxres)
   end do
 end do

 tres = True
 tres@txFontHeightF = 0.018
 tres@txAngleF      = 90.
 gsn_text_ndc(wks,"CESM-LE",.02,.875,tres)
 gsn_text_ndc(wks,"GFDL-CM3-LE",.02,.65,tres)
   
 tres@txAngleF      = 0.
 gsn_text_ndc(wks,"mean",.31,.99,tres)
 gsn_text_ndc(wks,"std dev",.77,.99,tres)

 tres@txFontHeightF = 0.02
 gsn_text_ndc(wks,"1970-1999",.53,.99,tres)

 tres@txFontHeightF = 0.014
 gsn_text_ndc(wks,"mm",.535,.5,tres)

 p1res = True
 p1res@gsnFrame            = False
 p1res@gsnPanelLabelBar    = False
 p1res@gsnPanelLeft        = .03
 p1res@gsnPanelTop         = .975
 p1res@gsnPanelBottom      = .52
 gsn_panel(wks,plot_maps,(/2,2/),p1res)
  
;***************************************************************************************

; plot 30-year shifts in mean and stdv by region for each model (one panel per model)

 res                   = True
 res@gsnDraw           = False
 res@gsnFrame          = False

 res@xyMarkLineMode    = "Markers"
 res@xyMarker          = 4

 Xmin = -0.0075*1000.0 ; -0.03
 Xmax = 0.0075*1000.0 ; 0.03
; XtickInt = 0.01
 Ymin = -0.025*1000.0 ; -0.03
 Ymax = 0.025*1000.0 ; 0.03
; YtickInt = 0.01
 ;XYvalues = (/-.015,-.005,.00, .005,.015/)*1000.0
 ;XYlabels = (/"-.015","-.005","0",".005",".015"/)*1000.0

 XYvalues = (/-15, -5, 0,  5, 15/)
 XYlabels = (/-15, -5, 0,  5, 15/)

; XYlabels = sprintf("%3.2f",XYvalues)
 
 res@trXMinF = Xmin
 res@trXMaxF = Xmax
 res@trYMinF = Ymin
 res@trYMaxF = Ymax

; res@tmXBMode = "Manual"
; res@tmXBTickStartF   = Xmin
; res@tmXBTickSpacingF = XtickInt
; res@tmXBTickEndF     = Xmax
; res@tmXBFormat       = "f4.2"
; res@tmYLMode = "Manual"
; res@tmYLTickStartF   = Ymin
; res@tmYLTickSpacingF = YtickInt
; res@tmYLTickEndF     = Ymax
; res@tmYLFormat       = "f4.2"

 res@tmXBMode = "Explicit"
 res@tmXBValues = XYvalues
 res@tmXBLabels = XYlabels
 res@tmYLMode = "Explicit"
 res@tmYLValues = XYvalues
 res@tmYLLabels = XYlabels
 
 res@tmXBMajorOutwardLengthF = 0.0               ; draw tickmarks inward
 res@tmXBMinorOutwardLengthF = 0.0               ; draw minor ticks inward
 res@tmYLMajorOutwardLengthF = 0.0               ; draw tickmarks inward
 res@tmYLMinorOutwardLengthF = 0.0               ; draw minor ticks inward

 res@tmXMajorGrid                = True          ; implement x grid 
 res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
 res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
 res@tmYMajorGrid                = True          ; implement y grid
 res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
 res@tmYMajorGridLineDashPattern = 2             ; select short dash lines

 res@tmXBLabelFontHeightF  = 0.035
 res@tmYLLabelFontHeightF  = 0.035

 res@tiXAxisFontHeightF    = 0.035
 res@tiYAxisFontHeightF    = 0.035
 res@tiMainFontHeightF     = 0.05
  
 res@tiXAxisString  = "~F33~D~F21~ std dev"
 res@tiYAxisString  = "~F33~D~F21~ mean"

 cmap = read_colormap_file("MPL_rainbow")
 idc = (/2,40,93,127/)
 clevs = cmap(idc,:)

 eblnres = True
 eblnres@tfPolyDrawOrder = "PreDraw"
 eblnres@gsLineOpacityF = .3
 eblnres@gsLineThicknessF = 4.

 emlnres = True
 emlnres@tfPolyDrawOrder = "PostDraw"
 emlnres@gsLineThicknessF = 25.

 txres = True
 txres@txFontHeightF = 0.035

 admnres = True
 admnres@amJust           = "TopRight"
 admnres@amParallelPosF   =  0.47    ; 0.5 is the right edge of the plot.
 admnres@amOrthogonalPosF = -0.47    ; -0.5 is the top edge of the plot.
 
 adsdres = admnres
 adsdres@amOrthogonalPosF = -0.39
 
 ; Panel #2 ===========================================
 plot_scat_cesm = new(nregs,graphic)
 oplot_cesm = plot_scat_cesm
 emblinesegs = new((/nregs,nps-1,nEnsc/),graphic)
 emnlinesegs = new((/nregs,nps-1/),graphic)

 txmn = new(nregs,graphic)
 txsd = new(nregs,graphic)
 admn = new(nregs,graphic)
 adsd = new(nregs,graphic)

 do r = 0,nregs-1
;   res@tiMainString = region(r)
   res@xyMarkerColor = clevs(0,:) ;"purple"
   plot_scat_cesm(r) = gsn_csm_xy(wks,cesm_sd_shfts(0,r,:),cesm_mn_shfts(0,r,:),res)
    ;plot_scat_cesm(r) = gsn_csm_xy(wks,cesm_sd_shfts(0,r,:)*10.0,cesm_mn_shfts(0,r,:),res)

   do p = 0,nps-2
     eblnres@gsLineColor = clevs(p,:)
     do i = 0,nEnsc-1
       emblinesegs(r,p,i) = gsn_add_polyline(wks,plot_scat_cesm(r),(/cesm_sd_shfts(p,r,i),cesm_sd_shfts(p+1,r,i)/),\
                                             (/cesm_mn_shfts(p,r,i),cesm_mn_shfts(p+1,r,i)/),eblnres)
     end do
   end do
   do p = 0,nps-2
     emlnres@gsLineColor = clevs(p,:)
     emnlinesegs(r,p) = gsn_add_polyline(wks,plot_scat_cesm(r),(/em_cesm_sd_shfts(p,r),em_cesm_sd_shfts(p+1,r)/),\
                                             (/em_cesm_mn_shfts(p,r),em_cesm_mn_shfts(p+1,r)/),emlnres)
   end do

   res@xyMarkerColor = clevs(nps-2,:) ;"red"
   oplot_cesm(r) = gsn_csm_xy(wks,cesm_sd_shfts(nps-1,r,:),cesm_mn_shfts(nps-1,r,:),res)
   overlay(plot_scat_cesm(r),oplot_cesm(r))

   txmn(r) = gsn_create_text(wks,"~F33~m ~F21~= "+sprintf("%5.1f",em_cesm_mn_ref(r)),txres)
   txsd(r) = gsn_create_text(wks,"~F33~s ~F21~= "+sprintf("%5.1f",em_cesm_sd_ref(r)),txres)

   admn(r) = gsn_add_annotation(plot_scat_cesm(r),txmn(r),admnres)
   adsd(r) = gsn_add_annotation(plot_scat_cesm(r),txsd(r),adsdres)
 end do
 
 p2res = True
 p2res@gsnFrame            = False
 p2res@gsnPanelLeft        = .04
 p2res@gsnPanelRight       = .975
 p2res@gsnPanelTop         = .45
 p2res@gsnPanelBottom      = .26
 gsn_panel(wks,plot_scat_cesm,(/1,nregs/),p2res)

 tres@txFontHeightF = 0.015
 ypos = .455
 xpos = (/.157,.34,.53,.716,.904/)
 do r = 0,nregs-1
   gsn_text_ndc(wks,region(r),xpos(r),ypos,tres)
 end do
  
 ; Panel #3 ===========================================
 plot_scat_gfdl = new(nregs,graphic)
 oplot_gfdl = plot_scat_gfdl
 emblinesegs := new((/nregs,nps-1,nEnsg/),graphic)
 emnlinesegs = new((/nregs,nps-1/),graphic)

 txmn = new(nregs,graphic)
 txsd = new(nregs,graphic)
 admn = new(nregs,graphic)
 adsd = new(nregs,graphic)

 do r = 0,nregs-1
   res@xyMarkerColor = clevs(0,:) ;"purple"
   plot_scat_gfdl(r) = gsn_csm_xy(wks,gfdl_sd_shfts(0,r,:),gfdl_mn_shfts(0,r,:),res)

   do p = 0,nps-2
     eblnres@gsLineColor = clevs(p,:)
     do i = 0,nEnsg-1
       emblinesegs(r,p,i) = gsn_add_polyline(wks,plot_scat_gfdl(r),(/gfdl_sd_shfts(p,r,i),gfdl_sd_shfts(p+1,r,i)/),\
                                             (/gfdl_mn_shfts(p,r,i),gfdl_mn_shfts(p+1,r,i)/),eblnres)
     end do
   end do
   do p = 0,nps-2
     emlnres@gsLineColor = clevs(p,:)
     emnlinesegs(r,p) = gsn_add_polyline(wks,plot_scat_gfdl(r),(/em_gfdl_sd_shfts(p,r),em_gfdl_sd_shfts(p+1,r)/),\
                                             (/em_gfdl_mn_shfts(p,r),em_gfdl_mn_shfts(p+1,r)/),emlnres)
   end do

   res@xyMarkerColor = clevs(nps-2,:) ;"red"
   oplot_gfdl(r) = gsn_csm_xy(wks,gfdl_sd_shfts(nps-1,r,:),gfdl_mn_shfts(nps-1,r,:),res)
   overlay(plot_scat_gfdl(r),oplot_gfdl(r))

   txmn(r) = gsn_create_text(wks,"~F33~m ~F21~= "+sprintf("%5.1f",em_gfdl_mn_ref(r)),txres)
   txsd(r) = gsn_create_text(wks,"~F33~s ~F21~= "+sprintf("%5.1f",em_gfdl_sd_ref(r)),txres)

   admn(r) = gsn_add_annotation(plot_scat_gfdl(r),txmn(r),admnres)
   adsd(r) = gsn_add_annotation(plot_scat_gfdl(r),txsd(r),adsdres)
 end do

 p3res = True
 p3res@gsnFrame            = False
 p3res@gsnPanelLeft        = .04
 p3res@gsnPanelRight       = .975
 p3res@gsnPanelTop         = .27
 p3res@gsnPanelBottom      = .08
 gsn_panel(wks,plot_scat_gfdl,(/1,nregs/),p3res)

 ; Labelbar for scatterplots ===========================================
 lbres                    = True
 lbres@vpWidthF           = 0.5
 lbres@vpHeightF          = 0.06
 lbres@lbOrientation      = "Horizontal"
 lbres@lbLabelAlignment   = "ExternalEdges"
 lbres@lbFillColors       = clevs
 lbres@lbMonoFillPattern  = True             ; Fill them all solid.
 lbres@lbLabelFontHeightF = 0.012
 lbres@lbTitleString      = "Shifts between 30-year periods"
 lbres@lbTitleFontHeightF = 0.012
 lbres@lbTitlePosition    = "Top"
 lbres@lbPerimOn          = False
 gsn_labelbar_ndc(wks,nps-1,periods,0.28,0.06,lbres)  

 tres@txAngleF      = 90.
 tres@txFontHeightF = 0.018
 gsn_text_ndc(wks,"CESM-LE",.02,.375,tres)
 gsn_text_ndc(wks,"GFDL-CM3-LE",.02,.185,tres)
   
; drawNDCGrid(wks)

 frame(wks)
 system("convert -trim "+plotname+"."+plottype+" "+plotname+"."+plottype)
