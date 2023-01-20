;***************************************************************************************
; plot 30-yr shifts in Nino3.4 SST means and standard deviations

 plottype = "png"
 ;plottype@wkPaperWidthF  =  8.5  ; in inches
 ;plottype@wkPaperHeightF = 14.0  ; in inches

 ;plottype@wkWidth  = 4200
 ;plottype@wkHeight = 4200

 ;dirout = "/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/plots/"
 dirout = "/glade/work/kumar34/hydroclimate_ncc/"
 plotname = dirout+"SM_SE_MeansStdv_30yrShifts_rev4"
 system("rm -f "+plotname+"."+plottype)
 wks  = gsn_open_wks (plottype,plotname)

;***************************************************************************************

; specify start and end dates for each 30-yr period
 nyrp = 30 ; length of period in yrs
 s0 = 194001 ; starting time
 wstep = nyrp*100
 nps = 5 ; number of nyrp-yr periods
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
; MAKE NINO INDICES FROM MODEL RUN
 
; identify Nino3 or Nino3.4 region

 ninoreg = "SW SM" ; "Nino3"
 ntilde  = "SW SM" ; 
 ninolab = "SW SM" ; "Ni"+ntilde+"o 3"

 latS = 25.
 latN = 35.

 ; GCMs have different longitude grids
 lonWc = 265 ; -170 for Nino3.4, -150 for Nino3
 lonEc = 282 ; -120 for Nino3.4, -90 for Nino3

 lonWg = -360.0+265. ; -170 for Nino3.4, -150 for Nino3
 lonEg = -360.0+282. ; -120 for Nino3.4, -90 for Nino3

;***************************************************************************************
 url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"


; read in region's SST from GFDL-CM3-LE, computing means and stddevs
 ;dirgfdl = "/Projects/SERDP_regimeshifts/data/GFDL-CM3-LE/SST/"
 ;dirgfdl = url+"GFDL-CM3-LE/SST/"
 ;filegfdl = dirgfdl+"SST_GFDL-CM3-LE_20em_192001-210012_12moRunningMeans.nc"

  dirgfdl = url+"GFDL-CM3-LE/SOILM/"
  filegfdl = dirgfdl+"SOILM_GFDL-CM3-LE_20em_192001-210012_12moRunningMeans.nc"



 print(filegfdl)

 ;exit

 flg = addfile(filegfdl,"r")

 lat_g = flg->lat
 nlatg = dimsizes(lat_g)
 lon_g = flg->lon
 nlong = dimsizes(lon_g)

 date = flg->date

 ensgi = ispan(1,20,1)
 nEnsg = dimsizes(ensgi)

 mean_gfdl = new((/nps,nEnsg/),"float",1e+36)
 stdv_gfdl = mean_gfdl
 do p = 0,nps-1
   print_table([/yrstr(p)+"-"+yrend(p)/],"%s")
   sti = ind(date.eq.strdates(p))
   endi = ind(date.eq.enddates(p))
   means = flg->SOILM(sti:endi,:,{latS:latN},{lonWc:lonEc})
   means = (/means*1000.0/)
   regtsr = dim_avg_n_Wrap(means,(/2,3/))
   mn = dim_avg_n(regtsr,0)
   mean_gfdl(p,:) = mn
   sd = dim_stddev_n(regtsr,0)
   stdv_gfdl(p,:) = sd
   delete([/means,regtsr,mn,sd/])
 end do

 printVarSummary(mean_gfdl)
 printMinMax(mean_gfdl,True)
 printVarSummary(stdv_gfdl)
 printMinMax(stdv_gfdl,True)
 
 ; compute ensemble mean for reference period (1970-99)
 em_gfdl_mn_ref = avg(mean_gfdl(1,:))
 em_gfdl_sd_ref = (avg(stdv_gfdl(1,:)^2))^0.5

 ; compute shifts in mean and stdv, departing from the reference period ensemble mean
 gfdl_mn_shfts = mean_gfdl
 gfdl_sd_shfts = stdv_gfdl
 do i = 0,nEnsg-1
   do p = 0,nps-1
     gfdl_mn_shfts(p,i) = mean_gfdl(p,i) - em_gfdl_mn_ref
     gfdl_sd_shfts(p,i) = (stdv_gfdl(p,i) - em_gfdl_sd_ref)*100.0/em_gfdl_sd_ref
   end do
 end do

 em_gfdl_mn_shfts = dim_avg_n(gfdl_mn_shfts,1)


 em_gfdl_sd_shfts = dim_avg_n(gfdl_sd_shfts,1)

 se_gfdl_mn_shfts = dim_stddev_n(gfdl_mn_shfts,1)*2.0/(20^0.5)
 se_gfdl_sd_shfts = dim_stddev_n(gfdl_sd_shfts,1)*2.0/(20^0.5)


 print(em_gfdl_mn_shfts)
 print(se_gfdl_mn_shfts) 
 printMinMax(gfdl_mn_shfts,True)
 

 print(em_gfdl_sd_shfts)
 print(se_gfdl_sd_shfts)
 printMinMax(gfdl_sd_shfts,True)

 ;exit

;***************************************************************************************

; read in SST from CESM-LE, computing means and stddevs on map and per region
 ;dircesm = "/Projects/SERDP_regimeshifts/data/CESM-LE/ocn/postprocessed/"

 ;dircesm = url+"CESM-LE/ocn/postprocessed/"
 ;filecesm = dircesm+"IndoPacific_SST_CESM-LE_40members_1920-2100_12moRunMeans.nc"

 dircesm = url+"CESM-LE/lnd/postprocessed/"
 filecesm = dircesm+"H2OSOI_CESM-LE_40em_192001-210012_12moRunningMeans.nc"

 ;flc = addfile(filecesm,"r")
 ;varname_cesm = "H2OSOI"


 flc = addfile(filecesm,"r")

 ;tlat = flc->tlat
 ;tlon = flc->tlon
 ;ji_reg = region_ind(tlat,tlon,latS,latN,lonWc,lonEc)
 ;jStrt = ji_reg(0)      ; lat start
 ;jLast = ji_reg(1)      ; lat last
 ;iStrt = ji_reg(2)      ; lon start
 ;iLast = ji_reg(3)      ; lon last

 ;nlatc = jLast - jStrt + 1
 ;nlonc = iLast - iStrt + 1
 
 lat_c = flc->lat
 nlatc = dimsizes(lat_c)
 lon_c = flc->lon
 nlonc = dimsizes(lon_c)
 datec = flc->date

 ensci = ispan(1,40,1)
 nEnsc = dimsizes(ensci)

 mean_cesm = new((/nps,nEnsc/),"float",1e+36)
 stdv_cesm = mean_cesm

 do p = 0,nps-1
   print_table([/yrstr(p)+"-"+yrend(p)/],"%s")
   sti = ind(datec.eq.strdates(p))
   endi = ind(datec.eq.enddates(p))
   ;means = flc->H2OSOI(sti:endi,:,jStrt:jLast,iStrt:iLast)
   ;regtsr = dim_avg_n_Wrap(means,(/2,3/))

   means = flc->H2OSOI(sti:endi,:,{latS:latN},{lonWc:lonEc})
   means = (/means*1000.0/)
   regtsr = dim_avg_n_Wrap(means,(/2,3/))


   mn = dim_avg_n(regtsr,0)
   mean_cesm(p,:) = mn
   sd = dim_stddev_n(regtsr,0)
   stdv_cesm(p,:) = sd
   delete([/means,regtsr,mn,sd/])
 end do

 printVarSummary(mean_cesm)
 printMinMax(mean_cesm,True)
 printVarSummary(stdv_cesm)
 printMinMax(stdv_cesm,True)

 ; compute ensemble mean for reference period (1970-99)
 em_cesm_mn_ref = avg(mean_cesm(1,:))
 em_cesm_sd_ref = (avg(stdv_cesm(1,:)^2))^0.5

; compute shifts in mean and stdv for each ensemble member, as differences from
; the reference period ensemble mean
 cesm_mn_shfts = mean_cesm
 cesm_sd_shfts = stdv_cesm
 do i = 0,nEnsc-1
   do p = 0,nps-1
     cesm_mn_shfts(p,i) = mean_cesm(p,i) - em_cesm_mn_ref
     cesm_sd_shfts(p,i) = (stdv_cesm(p,i) - em_cesm_sd_ref)*100/em_cesm_sd_ref
   end do
 end do

 em_cesm_mn_shfts = dim_avg_n(cesm_mn_shfts,1)
 em_cesm_sd_shfts = dim_avg_n(cesm_sd_shfts,1)


 se_cesm_mn_shfts = dim_stddev_n(cesm_mn_shfts,1)*2.0/(40^0.5)
 se_cesm_sd_shfts = dim_stddev_n(cesm_sd_shfts,1)*2.0/(40^0.5)



 print(em_cesm_mn_shfts)
 print(se_cesm_mn_shfts)
 printMinMax(cesm_mn_shfts,True)


 print(em_cesm_sd_shfts)
 print(se_cesm_sd_shfts)
 printMinMax(cesm_sd_shfts,True)

 ;exit

 
;****************************************************************************

; plot 30-year shifts in mean and stdv for each model

 res                   = True
 res@gsnDraw           = False
 res@gsnFrame          = False
 ;res@gsnMaximize = True


 res@xyMarkLineMode    = "Markers"
 res@xyMarker          = 4

 Xmin = -60.
 Xmax = 60.
; XtickInt = 0.5
 ;Xvalues = (/-20.,-.5,0.,.5,1./)
 ;Xlabels = (/"-1.","-.5","0.",".5","1."/)
 ;Xvalues = (/-50.0, -30.0, -40.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0/)
  Xvalues = (/-40.0, -20.0, 0.0, 20.0, 40.0/)
 ;Xlabels = (/-25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0/)
 ;Xlabels = (/-20.0, -10.0, 0.0, 10.0, 20.0/)
  ;Xlabels = (/' ', -40.0,' ', -20.0, ' ', 0.0, ' ', 20.0, ' ', 40.0, ' '/)
 
 Xlabels = (/-40.0, -20.0, 0.0, 20.0, 40.0/)



; Ymin = -.5
; Ymax = 4.
; YtickInt = .5
; Yvalues = (/.0,1.,2.,3.,4./)
; Ylabels = (/"0.","1.","2.","3.","4."/)
 ; XtickInt = 0.01
 Ymin = -0.010*1000.0 ; -0.03
 Ymax = 0.020*1000.0 ; 0.03
; YtickInt = 0.01
 ;XYvalues = (/-.015,-.005,.00, .005,.015/)*1000.0
 ;XYlabels = (/"-.015","-.005","0",".005",".015"/)*1000.0

 Yvalues = (/-10.0, 0.0, 10.0, 20.0/)
 Ylabels = (/-10.0, 0.0, 10.0, 20.0/)




 res@trXMinF = Xmin
 res@trXMaxF = Xmax
 res@trYMinF = Ymin
 res@trYMaxF = Ymax

; res@tmXBMode = "Manual"
; res@tmXBTickStartF   = Xmin
; res@tmXBTickSpacingF = XtickInt
; res@tmXBTickEndF     = Xmax
; res@tmXBFormat       = "f3.1"
; res@tmYLMode = "Manual"
; res@tmYLTickStartF   = Ymin
; res@tmYLTickSpacingF = YtickInt
; res@tmYLTickEndF     = Ymax
; res@tmYLFormat       = "f2.0"

 res@tmXBMode = "Explicit"
 res@tmXBValues = Xvalues
 res@tmXBLabels = Xlabels
 ;res@tmXBMinorValues = (/-.875,-.75,-.625,-.375,-.25,-.125,.125,.25,.375,.625,.75,.875/)
 res@tmYLMode = "Explicit"
 res@tmYLValues = Yvalues
 res@tmYLLabels = Ylabels
 ;res@tmYLMinorValues = (/.25,.5,.75,1.25,1.5,1.75,2.25,2.5,2.75,3.25,3.5,3.75/)
 
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

 res@tmXBLabelFontHeightF  = 0.025
 res@tmYLLabelFontHeightF  = 0.025

 res@tiXAxisFontHeightF    = 0.025
 res@tiYAxisFontHeightF    = 0.025
  
 res@tiXAxisString  = "~F33~D~F21~ std dev (%)"
 res@tiYAxisString  = "~F33~D~F21~ mean (mm)"

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
 txres@txFontHeightF = 0.025

 admnres = True
 admnres@amJust           = "TopLeft"
 admnres@amParallelPosF   =  -0.48    ; 0.5 is the right edge of the plot.
 admnres@amOrthogonalPosF = -0.47    ; -0.5 is the top edge of the plot.
 
 adsdres = admnres
 adsdres@amOrthogonalPosF = -0.39
 
 plot_scat = new(2,graphic)
 oplot = plot_scat
 Cemblinesegs = new((/nps-1,nEnsc/),graphic)
 Gemblinesegs = new((/nps-1,nEnsg/),graphic)
 emnlinesegs = new((/2,nps-1/),graphic)

 txmn = new(2,graphic)
 txsd = new(2,graphic)
 admn = new(2,graphic)
 adsd = new(2,graphic)

 ; CESM plot
; res@tiMainString = "CESM-LE"
 res@xyMarkerColor = clevs(0,:) ;"purple"
 plot_scat(0) = gsn_csm_xy(wks,cesm_sd_shfts(0,:),cesm_mn_shfts(0,:),res)

 do p = 0,nps-2
   eblnres@gsLineColor = clevs(p,:)
   do i = 0,nEnsc-1
     Cemblinesegs(p,i) = gsn_add_polyline(wks,plot_scat(0),(/cesm_sd_shfts(p,i),cesm_sd_shfts(p+1,i)/),\
                                         (/cesm_mn_shfts(p,i),cesm_mn_shfts(p+1,i)/),eblnres)
   end do
 end do

 do p = 0,nps-2
   emlnres@gsLineColor = clevs(p,:)
   emnlinesegs(0,p) = gsn_add_polyline(wks,plot_scat(0),(/em_cesm_sd_shfts(p),em_cesm_sd_shfts(p+1)/),\
                                      (/em_cesm_mn_shfts(p),em_cesm_mn_shfts(p+1)/),emlnres)
 end do

 res@xyMarkerColor = clevs(nps-2,:) ;"red"
 oplot(0) = gsn_csm_xy(wks,cesm_sd_shfts(nps-1,:),cesm_mn_shfts(nps-1,:),res)
 overlay(plot_scat(0),oplot(0))

 txmn(0) = gsn_create_text(wks,"~F33~m ~F21~= "+sprintf("%5.2f",em_cesm_mn_ref),txres)
 txsd(0) = gsn_create_text(wks,"~F33~s ~F21~= "+sprintf("%5.2f",em_cesm_sd_ref),txres)

 admn(0) = gsn_add_annotation(plot_scat(0),txmn(0),admnres)
 adsd(0) = gsn_add_annotation(plot_scat(0),txsd(0),adsdres)

 ; GFDLCM3 plot
; res@tiMainString = "GFDL-CM3-LE"
 res@xyMarkerColor = clevs(0,:) ;"purple"
 plot_scat(1) = gsn_csm_xy(wks,gfdl_sd_shfts(0,:),gfdl_mn_shfts(0,:),res)

 do p = 0,nps-2
   eblnres@gsLineColor = clevs(p,:)
   do i = 0,nEnsg-1
     Gemblinesegs(p,i) = gsn_add_polyline(wks,plot_scat(1),(/gfdl_sd_shfts(p,i),gfdl_sd_shfts(p+1,i)/),\
                                         (/gfdl_mn_shfts(p,i),gfdl_mn_shfts(p+1,i)/),eblnres)
   end do
 end do

 do p = 0,nps-2
   emlnres@gsLineColor = clevs(p,:)
   emnlinesegs(1,p) = gsn_add_polyline(wks,plot_scat(1),(/em_gfdl_sd_shfts(p),em_gfdl_sd_shfts(p+1)/),\
                                      (/em_gfdl_mn_shfts(p),em_gfdl_mn_shfts(p+1)/),emlnres)
 end do

 res@xyMarkerColor = clevs(nps-2,:) ;"red"
 oplot(1) = gsn_csm_xy(wks,gfdl_sd_shfts(nps-1,:),gfdl_mn_shfts(nps-1,:),res)
 overlay(plot_scat(1),oplot(1))

 txmn(1) = gsn_create_text(wks,"~F33~m ~F21~= "+sprintf("%5.2f",em_gfdl_mn_ref),txres)
 txsd(1) = gsn_create_text(wks,"~F33~s ~F21~= "+sprintf("%5.2f",em_gfdl_sd_ref),txres)

 admn(1) = gsn_add_annotation(plot_scat(1),txmn(1),admnres)
 adsd(1) = gsn_add_annotation(plot_scat(1),txsd(1),adsdres)

 pres = True
 pres@gsnFrame            = False
; pres@gsnPanelLeft        = .04
; pres@gsnPanelRight       = .975
; pres@gsnPanelTop         = .45
; pres@gsnPanelBottom      = .26
 gsn_panel(wks,plot_scat,(/1,2/),pres)

 ; Labelbar
 lbres                    = True
 lbres@vpWidthF           = 0.7
 lbres@vpHeightF          = 0.08
 lbres@lbOrientation      = "Horizontal"
 lbres@lbLabelAlignment   = "ExternalEdges"
 lbres@lbFillColors       = clevs
 lbres@lbMonoFillPattern  = True             ; Fill them all solid.
 lbres@lbLabelFontHeightF = 0.015
 lbres@lbTitleString      = "Shifts between 30-year periods"
 lbres@lbTitleFontHeightF = 0.018
 lbres@lbTitlePosition    = "Top"
 lbres@lbPerimOn          = False
 gsn_labelbar_ndc(wks,nps-1,periods,0.175,0.215,lbres)  

 tres = True
 tres@txFontHeightF = 0.025
 gsn_text_ndc(wks,"CESM-LE",.275,.76,tres)
 gsn_text_ndc(wks,"GFDL-CM3-LE",.775,.76,tres)

; drawNDCGrid(wks)

 frame(wks)
 system("convert -trim "+plotname+"."+plottype+" "+plotname+"."+plottype)

 
