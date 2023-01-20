;***************************************************************************************
; plot Tropical SST leading EOF mode in CESM-LE (first and last 20yr periods)
;***************************************************************************************


   yrstr = (/1950,1970,1990,2010,2030,2050/)
   yrend = (/1969,1989,2009,2029,2049,2069/)
   periods = (/19501969,19701989,19902009,20102029,20302049,20502069/)

    ;yrstr = (/1940,1970,2000,2030,2060/)
    ;yrend = (/1969,1999,2029,2059,2089/)
    ;periods = (/19401969,19701999,20002029,20302059,20602089/)

   url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"

;***************************************************************************************

; read in variable from CESM-LE
   
   dir = url+"/CESM-LE/ocn/postprocessed/" + \
         "EOFs_21yrSmoothEnsembClim/TropicalEOFs/20yrEOFs/spatial_mode/"
   
   filelens = dir+"sst.eof_anomaly_monthlymean_CESM-LE_125x125.ens001-040.jan01dec31.1950-1969.nc"
   fl = addfile(filelens,"r")

   lat = fl->lat
   nlat = dimsizes(lat)
   lon = fl->lon
   nlon = dimsizes(lon)
  
   exit
 
;***************************************************************************************

   EOF = new((/2,nlat,nlon/),"float",1e+30)
   varexp = new((/2/),"float",1e+30)
   pp = (/0,5/) ; only read first and last periods
   do p = 0,1
     print_table([/yrstr(pp(p))+"-"+yrend(pp(p))/],"%s")
     filelens = dir+"sst.eof_anomaly_monthlymean_CESM-LE_125x125.ens001-040.jan01dec31."+ \
                yrstr(pp(p))+"-"+yrend(pp(p))+".nc"
     fl = addfile(filelens,"r")
     eof = fl->eof(0,0,:,:)
     if(p.eq.1) then
       eof = eof *-1.
     end if
     EOF(p,:,:) = eof
     printMinMax(eof,True)
     varexp(p) = fl->var(0)
     delete([/fl,eof/])
   end do

   delete_VarAtts(EOF,(/"name","input_years","index","level"/))
   printVarSummary(EOF)
   print(varexp)
   
;***************************************************************************************

   plottype = "x11"
   plottype@wkWidth  = 3072
   plottype@wkHeight = 3072

   ;dirout = "/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/plots/"

   dirout = "/glade/work/kumar34/hydroclimate_ncc/"

   plotname = dirout+"CESM-LE_SST_TropicalLeadingEOF"
   system("rm -f "+plotname+"."+plottype)
   wks  = gsn_open_wks (plottype,plotname)

   cmap_mn = read_colormap_file("BlueDarkRed18")

   res                   = True
   res@gsnDraw           = False
   res@gsnFrame          = False
   res@gsnAddCyclic      = True

   res@mpShapeMode       = "FreeAspect"
   res@vpWidthF          = .5
   res@vpHeightF         = .2

   res@mpLimitMode       = "LatLon"
   res@mpMinLatF         = -30.
   res@mpMaxLatF         = 30.
;   res@mpMinLonF         = 30.
;   res@mpMaxLonF         = 20.
   res@mpCenterLonF      = 205.
   res@mpOutlineBoundarySets = "Geophysical"
   res@mpGeophysicalLineColor = "black"
   
;   res@pmTickMarkDisplayMode = "Always"

   res@cnInfoLabelOn     = False
   res@cnFillOn          = True
   res@cnLinesOn         = False
   res@cnFillMode        = "RasterFill"
   res@cnRasterSmoothingOn = True
   res@cnLineLabelsOn    = False

   res@lbLabelBarOn      = False
   res@lbLabelFontHeightF = 0.01
   res@lbLeftMarginF     = -0.5
   res@lbRightMarginF    = -0.5
   res@lbTopMarginF      = -0.2
   res@lbBottomMarginF   = 0.3

   res@cnLevelSelectionMode = "ManualLevels"
   res@cnMinLevelValF = -.03
   res@cnMaxLevelValF = .03
   res@cnLevelSpacingF = .005
   res@cnFillPalette     = cmap_mn

   ; set coords for Nino box
   la = (/-5.,5.,5.,-5.,-5./)
;   lo = (/210.,210.,270.,270.,210./) ; Nino3
   lo = (/190.,190.,240.,240.,190./) ; Nino3.4

   bxres = True
   bxres@gsLineThicknessF = 4.
   bxres@gsLineColor      = "black"

   lbres = True
   lbres@txPerimOn        = True
   lbres@txFontHeightF    = 0.012
   lbres@txBackgroundFillColor = "white"

   anres = True
   anres@amParallelPosF   = -0.5
   anres@amOrthogonalPosF = -0.5
   anres@amJust           = "TopLeft"
   
   plot = new(2,graphic)
   box = new(2,graphic)
   lbid = new(2,graphic)
   anid = new(2,graphic)
   do p = 0,1
     if (p.eq.1) then
       res@lbLabelBarOn = True
       res@tmXBLabelsOn = False
     end if
     plot(p) = gsn_csm_contour_map(wks,EOF(p,:,:),res)
     box(p) = gsn_add_polyline(wks,plot(p),lo,la,bxres)
     lbid(p) = gsn_create_text(wks,"var exp = "+sprintf("%4.1f",varexp(p))+"%",lbres)
     anid(p) = gsn_add_annotation(plot(p),lbid(p),anres)
   end do

   tres = True
   tres@txAngleF      = 90.          ; Rotate string clockwise 
   tres@txFontHeightF = 0.025
   gsn_text_ndc(wks,yrstr(0)+"-"+yrend(0),.04,.71,tres)
   gsn_text_ndc(wks,yrstr(5)+"-"+yrend(5),.04,.31,tres)
   
   tres@txAngleF      = 0.
   gsn_text_ndc(wks,"Leading EOF Mode of Tropical Oceans",.545,.915,tres)

;   drawNDCGrid(wks)

   pres = True
   pres@gsnPanelLabelBar    = False
   pres@gsnPanelLeft        = .05
   gsn_panel(wks,plot,(/2,1/),pres)

;   getvalues plot(0)
;    "vpWidthF"  : vpw
;    "vpHeightF" : vph
;   end getvalues
;   print(vpw)
;   print(vph)

   system("convert -trim "+plotname+"."+plottype+" "+plotname+"."+plottype)
   
   delete([/wks,res,tres,pres,plot/])

