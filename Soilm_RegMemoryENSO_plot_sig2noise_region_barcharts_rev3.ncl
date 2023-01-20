
url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"


;***************************************************************************************
; plot Fig. 5

; path to CESM-LE files
; dircesm = "/Projects/SERDP_regimeshifts/data/CESM-LE/lnd/postprocessed/regions_clean/"

; path to GFDL-CM3-LE files
; dirgfdl = "/Projects/SERDP_regimeshifts/data/GFDL-CM3-LE/SOILM/"

dircesm = url+"CESM-LE/lnd/postprocessed/regions_clean/"

dirgfdl =  url+"GFDL-CM3-LE/SOILM/"



;***************************************************************************************

; selected regions
 
; region = (/"Southwest","Southeast","SouthernUSPlains","NorthernUSPlains","CanadianPlains"/)
  region = (/"Southwest","Southeast", "CanadianPlains"/)
; regtitle = (/"Southwest","Southeast"," Southern~C~US Plains"," Northern~C~US Plains","Canadian~C~   Plains"/)
; regtitle = (/"Southwest","Southeast"," Southern US Plains"," Northern US Plains","Canadian Plains"/)
  regtitle = (/"Southwest","Southeast", "Canadian Plains"/)
 nregs = dimsizes(region)

;***************************************************************************************

; specify 30-yr periods
 nyrp = 30
 yrstr = (/1940,1970,2000,2030,2060/)
 yrend = yrstr + nyrp - 1
; periods = (/"1940-~C~1969","1970-~C~1999","2000-~C~2029","2030-~C~2059","2060-~C~2089"/)
 periods = (/"1940-69","1970-99","2000-29","2030-59","2060-89"/)
 nps = dimsizes(yrstr)
 
;***************************************************************************************
 
; plot signal-to-noise bar charts for all 3 linear models
; one row per period, one column for each GCM

 plottype = "x11"
 ;plottype@wkWidth  = 4200
 ;plottype@wkHeight = 4200

 ;dirout = "/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/plots/"
 dirout = "/glade/work/kumar34/hydroclimate_ncc/"
 plotname = dirout+"Soilm_RegMemoryENSO_sig2total_region_barcharts_rev4"
 system("rm -f "+plotname+"."+plottype)
 wks  = gsn_open_wks (plottype,plotname)

 res = True
 res@gsnDraw   = False
 res@gsnFrame  = False
 res@gsnXYBarChart = True
 res@gsnXYBarChartBarWidth = 0.4

 ; VP plotting positions (taken from the sig2noise maps panel plot and adjusted)
 vpx = (/.15,.48/)
; vpx = (/0.152784,0.53324/)
 vpy = (/.975,0.803,0.631,0.459,0.287/)
; vpy = (/0.978158,0.790263,0.602368,0.414474,0.226579/)
 vpwidth = .32 ; 0.343486
 vpheight = .165 ; 0.166074
 
 res@vpWidthF  = vpwidth
 res@vpHeightF = vpheight

 res@trXMinF   = 0.9
 res@trXMaxF   = 11.1
 res@trYMinF   = 0.
 res@trYMaxF   = 0.80

 res@tmXBMode = "Explicit"
 res@tmXBValues = (/2,4,6,8,10/)
 res@tmXBLabels = periods
 res@tmXBLabelFontHeightF = 0.009
 res@tmXBLabelDeltaF = -.6
; res@tmYLPrecision = 1
 res@tmYLLabelFontHeightF = 0.008
 res@tmYLLabelDeltaF = -.6
 res@tmXTOn = False
 
; colors = (/"gold","steelblue1","tomato"/)
 colors = (/"tan","cyan","blue"/)
 
 ; resources for error bars
 polyres               = True
 polyres@gsMarkerIndex = NhlNewMarker(wks,">",34,0.,0.,1.5,.7,0.)
 polyres@gsLineThicknessF = 5.

;***************************************************************************************

 random_setallseed(36494749,8494848) 

 do gcm = 0,1 
   res@vpXF = vpx(gcm)
   if(gcm.eq.1) then
     res@tmYLLabelsOn = False
   end if
   
   do r = 0,nregs-1 ; loop over regions
     res@vpYF = vpy(r)
     if(r.lt.nregs-1) then
       res@tmXBLabelsOn = False
     else
       res@tmXBLabelsOn = True
     end if
          
     if(gcm.eq.0) then
      ; read in regression metrics from CESM-LE
       filecesm = dircesm+"H2OSOI_CESM-LE_1940-2089_30yr_MemoryENSOregressions_"+region(r)+".nc"
       f = addfile(filecesm,"r")
       ;print(f)
       ;exit
     else
      ; read in regression metrics from GFDL-CM3-LE
       filegfdl = dirgfdl+"SOILM_GFDL-CM3-LE_1940-2089_30yr_MemoryENSOregressions_"+region(r)+".nc"
       f = addfile(filegfdl,"r")
     end if

     sig2noise := f->sig2noise ; signal-to-noise from each regression modeli
     sig2noise = sig2noise^2.0
     sig2noise = (/(sig2noise/(1+sig2noise))/)
      sig2noise = f->actskill
     lvmod := f->lvarmod  ; predicted soil moisture from each regression model
     error := f->error    ; error from each regression model
     ;error = (/(error/(1+error))/)

     nEns = dimsizes(error&ens)

     ; bootstrap s2n computations to get 95% confidence interval
     ; sample random selection of 40 ensemble members (with replacement), 1000 x
     niter = 1000

     ; get bootstrap sample for this iteration
     ens_sample = new((/niter,40/),"integer")
     do b = 0,niter-1
       sample = toint(random_uniform(0,nEns,40))
       ens_sample(b,:) = sample
       delete(sample)
     end do

     s2nci = new((/3,nps,2/),"float",sig2noise@_FillValue)
     do rg = 0,2
       do p = 0,nps-1
         s2nboot = new((/niter/),"float")
         do b = 0,niter-1
           sample = ens_sample(b,:)
           lm = lvmod(rg,p,:,sample)
           lm := ndtooned(lm)
           er = error(rg,p,:,sample)
           er := ndtooned(er)
           s2n = sqrt(variance(lm)/variance(er))
           s2n = (/(s2n/(1+s2n))/)
           s2nboot(b) = s2n
           delete([/sample,lm,er,s2n/])
         end do
         qsort(s2nboot)
         s2nci(rg,p,0) = s2nboot(toint(niter*.025)) ; 2.5%ile
         s2nci(rg,p,1) = s2nboot(toint(niter*.975)) ; 97.5%ile
         delete(s2nboot)
       end do
     end do
      
     res@gsnXYBarChartColors = colors(0)
     plot1 = gsn_csm_xy(wks,fspan(1.6,9.6,nps),sig2noise(0,:),res)
     res@gsnXYBarChartColors = colors(1)
     plot2 = gsn_csm_xy(wks,fspan(2.,10.,nps),sig2noise(1,:),res)
     res@gsnXYBarChartColors = colors(2)
     plot3 = gsn_csm_xy(wks,fspan(2.4,10.4,nps),sig2noise(2,:),res)

     cilo1 = new((/nps/),graphic)
     errbar1 = new((/nps/),graphic)
     cihi1 = new((/nps/),graphic)
     cilo2 = new((/nps/),graphic)
     errbar2 = new((/nps/),graphic)
     cihi2 = new((/nps/),graphic)
     cilo3 = new((/nps/),graphic)
     errbar3 = new((/nps/),graphic)
     cihi3 = new((/nps/),graphic)
     x01 = 1.6
     do p = 0,nps-1
       cilo1(p)   = gsn_add_polymarker(wks,plot1,x01,s2nci(0,p,0),polyres)
       errbar1(p) = gsn_add_polyline(wks,plot1,(/x01,x01/),(/s2nci(0,p,0),s2nci(0,p,1)/),polyres)
       cihi1(p)   = gsn_add_polymarker(wks,plot1,x01,s2nci(0,p,1),polyres)
     
       cilo2(p)   = gsn_add_polymarker(wks,plot2,x01+.4,s2nci(1,p,0),polyres)
       errbar2(p) = gsn_add_polyline(wks,plot2,(/x01+.4,x01+.4/),(/s2nci(1,p,0),s2nci(1,p,1)/),polyres)
       cihi2(p)   = gsn_add_polymarker(wks,plot2,x01+.4,s2nci(1,p,1),polyres)

       cilo3(p)   = gsn_add_polymarker(wks,plot3,x01+.8,s2nci(2,p,0),polyres)
       errbar3(p) = gsn_add_polyline(wks,plot3,(/x01+.8,x01+.8/),(/s2nci(2,p,0),s2nci(2,p,1)/),polyres)
       cihi3(p)   = gsn_add_polymarker(wks,plot3,x01+.8,s2nci(2,p,1),polyres)
     
       x01 = x01 + 2.
     end do

     draw(plot1)
     draw(plot2)
     draw(plot3)
   end do
 end do
     
 tres = True
 tres@txFontHeightF = 0.015
 tres@txAngleF      = 0.
 gsn_text_ndc(wks,"CESM-LE",.31,.99,tres)
 gsn_text_ndc(wks,"GFDL-CM3-LE",.64,.99,tres)

 tres@txFontHeightF = 0.013
 tres@txAngleF      = 90.
 xpos = .1
 ypos = (/.89,.72,.547,.375,.205/)
 do r = 0,nregs-1
   gsn_text_ndc(wks,regtitle(r),xpos,ypos(r),tres)
 end do
 
 lbres = True
 lbres@vpWidthF           = 0.06           ; labelbar width
 lbres@vpHeightF          = 0.12           ; labelbar height
 lbres@lbBoxMajorExtentF  = 0.14          ; puts space between color boxes
 lbres@lbMonoFillPattern  = True          ; Solid fill pattern
 lbres@lbLabelFontHeightF = 0.014         ; font height. default is small
 lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
 lbres@lbJustification    = "CenterLeft"  ; left justify boxes
 lbres@lbPerimOn          = False

 labels = (/"Memory","ENSO","Memory+ENSO"/)
 xpos := (/.29,.42,.54/)
 ypos := 0.48

 do i = 0,2
   lbres@lbFillColors = colors(i)
   gsn_labelbar_ndc(wks,1,labels(i),xpos(i),ypos,lbres)
 end do
      
; drawNDCGrid(wks)

 frame(wks)
 system("convert -trim "+plotname+"."+plottype+" "+plotname+"."+plottype)
