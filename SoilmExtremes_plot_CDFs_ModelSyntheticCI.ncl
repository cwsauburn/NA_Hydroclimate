;***************************************************************************************
; Plots cumulative distributions of wet/dry events in CESM-LE and GFDLCM3-LE

; Severity of events in all periods is scaled by the stddev of event severity in 1970-99.
;***************************************************************************************

 url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"

; dircesm = "/Projects/SERDP_regimeshifts/data/CESM-LE/lnd/postprocessed/regions_clean/"
; dirgfdl = "/Projects/SERDP_regimeshifts/data/GFDL-CM3-LE/SOILM/"
 ;dir = "/Projects/SERDP_regimeshifts/data/"
 dir = url
 dirmod = (/"CESM-LE/lnd/postprocessed/regions_clean/","GFDL-CM3-LE/SOILM/"/)
 dirext = (/"WetCDF_rndmSM0_v3/","DroughtCDF_rndmSM0_v3/"/)

 varmodnames = (/"H2OSOI_CESM-LE_","SOILM_GFDL-CM3-LE_"/)
 extfilename = (/"1940-2089_30yr_PluvialCDFs","1940-2089_30yr_DroughtCDFs"/)
 
 models = (/"CESM","GFDL-CM3"/)
 exttype = (/"Pluvials","Droughts"/)

 ; Regions to evaluate/plot
 region = (/"Southwest","Southeast","SouthernUSPlains","NorthernUSPlains","CanadianPlains"/)
 region_title = (/"Southwest","Southeast","South US Plains","North US Plains","Canadian Plains"/)
 nregs = dimsizes(region)
 
 ; Periods of analyses
 yrstr = (/1940,1970,2000,2030,2060/)
 yrend = (/1969,1999,2029,2059,2089/)
 ctryr = (/1955,1985,2015,2045,2075/)
 periodlabels = (/"1940-69","1970-99","2000-29","2030-59","2060-89"/)
 ; pick a period to compare against 1970-99
 p = 4

 ; Percentile rank IDs for given # of iterations
 ; (will be used to find median and 95% confidence interval in synthetic data)
 nit = 1000
 med =  toint(round(50./100.*(nit+1),0))
 pile025 = toint(round(2.5/100.*(nit+1),0))
 pile975 = toint(round(97.5/100.*(nit+1),0))

 ; Percentile array for interpolation of synthetic CDFs
 ; (will be used to find median and 95% confidence interval across each percentile level)
 pctiles = tofloat(ispan(0,995,5))/10. ; = (/0,.5,1.,1.5,...,98.5,99.,99.5/)
 npcts = dimsizes(pctiles)
 pct025 = toint(2.5/100.*nit)
 pct500 = toint(50./100.*nit)
 pct975 = toint(97.5/100.*nit)


 P1_temp = new((/npcts, 7/), "float", -9999.0)
 P1_temp(:, 0) = pctiles


;***************************************************************************************

 do r = 0,nregs-1  ; loop through regions
   print_table([/region(r)/],"%s")

   do m = 0,1  ; loop through models
     print_table([/models(m)/],"%s")

     do e = 0,1  ; loop through extremes type (0 = Pluvials, 1 = Droughts)
       print_table([/exttype(e)/],"%s")

       ; read in CDFs for 1970-99 events from model's soil moisture
       ; compute stddev of severity of the events for this period and use it as scaling reference
       ; sevsigma = stddev(cdf_severity)
   
       statfile = dir+dirmod(m)+dirext(e)+varmodnames(m)+extfilename(e)+"_detrendedAnoms_"+region(r)+"_19701999.nc"
       sf = addfile(statfile,"r")
       probax := sf->prob_axis
       sevcdf := sf->cdf_severity
       durcdf := sf->cdf_duration
   
       sevsigma = stddev(sevcdf)
       sevcdf = sevcdf / sevsigma

       ; read in CDFs for 1970-99 events from synthetic soil moisture, scale severity by sevsigma and
       ; compute the confidence interval of full CDF at each percentile level

       statfile = dir+dirmod(m)+dirext(e)+varmodnames(m)+extfilename(e)+"_SyntheticDetrendedAnoms_"+region(r)+"_19701999.nc"
       sf = addfile(statfile,"r")
       nevs := sf->nevents_synt
       probaxsynt := sf->prob_axis_synt
       sevcdfsynt := sf->cdf_severity_synt
       sevcdfsynt = sevcdfsynt / sevsigma
       durcdfsynt := sf->cdf_duration_synt

       sevcdfsyntinterp = new((/npcts,nit/),"float")
       durcdfsyntinterp = sevcdfsyntinterp
       do i = 0,nit-1
         ne = nevs(i)
         probaxisit := probaxsynt(i,0:ne-1)
         durcdfit := durcdfsynt(i,0:ne-1)
         sevcdfit := sevcdfsynt(i,0:ne-1)

         sevsinterp = linint1(probaxisit,sevcdfit,False,pctiles,0)
         sevsinterp(npcts-1) = sevcdfit(ne-1)
         sevcdfsyntinterp(:,i) = sevsinterp
         dursinterp = linint1(probaxisit,durcdfit,False,pctiles,0)
         dursinterp = round(dursinterp,0)
         dursinterp(npcts-1) = durcdfit(ne-1)
         durcdfsyntinterp(:,i) = dursinterp
       end do
    
       ; get the median and 95% confidence interval of CDFs at each percentile level
       cdfsevsyntlohi = new((/2,npcts/),"float")
       cdfsevsyntmed = new((/npcts/),"float")
       cdfdursyntlohi = new((/2,npcts/),"float")
       cdfdursyntmed = new((/npcts/),"float")
       do pct = 0,npcts-1
         samplesev = sevcdfsyntinterp(pct,:)
         qsort(samplesev)
         cdfsevsyntlohi(0,pct) = samplesev(pct025)
         cdfsevsyntlohi(1,pct) = samplesev(pct975)
         cdfsevsyntmed(pct) = samplesev(pct500)
         sampledur = durcdfsyntinterp(pct,:)
         qsort(sampledur)
         cdfdursyntlohi(0,pct) = sampledur(pct025)
         cdfdursyntlohi(1,pct) = sampledur(pct975)
         cdfdursyntmed(pct) = sampledur(pct500)
       end do

       ; read in CDFs for events in the future 30-year period, detrended and with trend,
       ; from both model and synthetic soil moisture, and scale severity by sevsigma.
       ; For synthetic events, compute the confidence interval of full CDF at each percentile level.

       ; detrended
       statfile = dir+dirmod(m)+dirext(e)+varmodnames(m)+extfilename(e)+"_detrendedAnoms_"+region(r)+"_"+yrstr(p)+yrend(p)+".nc"
       sf = addfile(statfile,"r")
       futprobax_dtr := sf->prob_axis
       futsevcdf_dtr := sf->cdf_severity
       futsevcdf_dtr = futsevcdf_dtr / sevsigma
       futdurcdf_dtr := sf->cdf_duration

       ; detrended, synthetic
       statfile = dir+dirmod(m)+dirext(e)+varmodnames(m)+extfilename(e)+"_SyntheticDetrendedAnoms_"+region(r)+"_"+yrstr(p)+yrend(p)+".nc"
       sf = addfile(statfile,"r")
       futne_dtr := sf->nevents_synt
       futprobaxsynt_dtr := sf->prob_axis_synt
       futsevcdfsynt_dtr := sf->cdf_severity_synt
       futsevcdfsynt_dtr = futsevcdfsynt_dtr / sevsigma
       futdurcdfsynt_dtr := sf->cdf_duration_synt

       futsevcdfsyntinterp_dtr = new((/npcts,nit/),"float")
       futdurcdfsyntinterp_dtr = futsevcdfsyntinterp_dtr
       do i = 0,nit-1
         ne = futne_dtr(i)
         probaxisit := futprobaxsynt_dtr(i,0:ne-1)
         sevcdfit := futsevcdfsynt_dtr(i,0:ne-1)
         durcdfit := futdurcdfsynt_dtr(i,0:ne-1)

         sevsinterp = linint1(probaxisit,sevcdfit,False,pctiles,0)
         sevsinterp(npcts-1) = sevcdfit(ne-1)
         futsevcdfsyntinterp_dtr(:,i) = sevsinterp
         dursinterp = linint1(probaxisit,durcdfit,False,pctiles,0)
         dursinterp = round(dursinterp,0)
         dursinterp(npcts-1) = durcdfit(ne-1)
         futdurcdfsyntinterp_dtr(:,i) = dursinterp
       end do
       
       ; get the median and 95% confidence interval of CDFs at each percentile level
       futcdfsevsyntlohi_dtr = new((/2,npcts/),"float")
       futcdfsevsyntmed_dtr = new((/npcts/),"float")
       futcdfdursyntlohi_dtr = new((/2,npcts/),"float")
       futcdfdursyntmed_dtr = new((/npcts/),"float")
       do pct = 0,npcts-1
         samplesev = futsevcdfsyntinterp_dtr(pct,:)
         qsort(samplesev)
         futcdfsevsyntlohi_dtr(0,pct) = samplesev(pct025)
         futcdfsevsyntlohi_dtr(1,pct) = samplesev(pct975)
         futcdfsevsyntmed_dtr(pct) = samplesev(pct500)
         sampledur = futdurcdfsyntinterp_dtr(pct,:)
         qsort(sampledur)
         futcdfdursyntlohi_dtr(0,pct) = sampledur(pct025)
         futcdfdursyntlohi_dtr(1,pct) = sampledur(pct975)
         futcdfdursyntmed_dtr(pct) = sampledur(pct500)
       end do

       ; trended
       statfile = dir+dirmod(m)+dirext(e)+varmodnames(m)+extfilename(e)+"_trendedAnoms_"+region(r)+"_"+yrstr(p)+yrend(p)+".nc"
       sf = addfile(statfile,"r")
       futprobax_tr := sf->prob_axis_tanom
       futsevcdf_tr := sf->cdf_severity_tanom
       futsevcdf_tr = futsevcdf_tr / sevsigma
       futdurcdf_tr := sf->cdf_duration_tanom

       ; trended, synthetic
       statfile = dir+dirmod(m)+dirext(e)+varmodnames(m)+extfilename(e)+"_SyntheticTrendedAnoms_"+region(r)+"_"+yrstr(p)+yrend(p)+".nc"
       sf = addfile(statfile,"r")
       futne_tr := sf->nevents_synt_tanom
       futprobaxsynt_tr := sf->prob_axis_synt_tanom
       futsevcdfsynt_tr := sf->cdf_severity_synt_tanom
       futsevcdfsynt_tr = futsevcdfsynt_tr / sevsigma
       futdurcdfsynt_tr := sf->cdf_duration_synt_tanom

       futsevcdfsyntinterp_tr = new((/npcts,nit/),"float")
       futdurcdfsyntinterp_tr = futsevcdfsyntinterp_tr
       do i = 0,nit-1
         ne = futne_tr(i)
         probaxisit := futprobaxsynt_tr(i,0:ne-1)
         sevcdfit := futsevcdfsynt_tr(i,0:ne-1)
         durcdfit := futdurcdfsynt_tr(i,0:ne-1)

         sevsinterp = linint1(probaxisit,sevcdfit,False,pctiles,0)
         sevsinterp(npcts-1) = sevcdfit(ne-1)
         futsevcdfsyntinterp_tr(:,i) = sevsinterp
         dursinterp = linint1(probaxisit,durcdfit,False,pctiles,0)
         dursinterp = round(dursinterp,0)
         dursinterp(npcts-1) = durcdfit(ne-1)
         futdurcdfsyntinterp_tr(:,i) = dursinterp
       end do

       ; get the median and 95% confidence interval of CDFs at each percentile level
       futcdfsevsyntlohi_tr = new((/2,npcts/),"float")
       futcdfsevsyntmed_tr = new((/npcts/),"float")
       futcdfdursyntlohi_tr = new((/2,npcts/),"float")
       futcdfdursyntmed_tr = new((/npcts/),"float")
       do pct = 0,npcts-1
         samplesev = futsevcdfsyntinterp_tr(pct,:)
         qsort(samplesev)
         futcdfsevsyntlohi_tr(0,pct) = samplesev(pct025)
         futcdfsevsyntlohi_tr(1,pct) = samplesev(pct975)
         futcdfsevsyntmed_tr(pct) = samplesev(pct500)
         sampledur = futdurcdfsyntinterp_tr(pct,:)
         qsort(sampledur)
         futcdfdursyntlohi_tr(0,pct) = sampledur(pct025)
         futcdfdursyntlohi_tr(1,pct) = sampledur(pct975)
         futcdfdursyntmed_tr(pct) = sampledur(pct500)
       end do

;***************************************************************************************

       ;dirout = "/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/plots/CDFs/"
       dirout = "/glade/work/kumar34/hydroclimate_ncc/"
       typPlt = "x11"
       typPlt@wkWidth  = 3000
       typPlt@wkHeight = 3000

       plotname = dirout+"SoilmExtremes_CDFs_1970-99_vs_"+periodlabels(p)+"_"+region(r)+"_"+models(m)+"_"+exttype(e)
       system("rm -f "+plotname+"."+typPlt)
       wks  = gsn_open_wks (typPlt,plotname)

       res                          = True
       res@gsnDraw                  = False
       res@gsnFrame                 = False

       res@vpHeightF                = 0.4
       res@vpWidthF                 = 0.4
       res@vpXF                     = 0.1
       res@vpYF                     = 0.8

       res@trXMinF                  = 0
       res@trXMaxF                  = 20 ; maxsev,maxdur
       res@trYMinF                  = 50
       res@trYMaxF                  = 100

       res@tiYAxisFontHeightF       = .035
      ; res@tiXAxisFontHeightF       = .035
       res@tiYAxisString            = " "

       res@tmXMajorGrid                = True          ; implement x grid 
       res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
       res@tmXMajorGridLineDashPattern = 2             ; select short dash lines
       res@tmYMajorGrid                = True          ; implement y grid
       res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
       res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
       res@tmYLLabelsOn                = False

       resM = res
       resM@xyMarkLineModes          = "Markers"
       resM@xyMarkers                = 12
       resM@xyMarkerSizeF            = 0.008
       resM@xyMarkerThicknessF       = 5
       starcolor = (/"dodgerblue4","darkgoldenrod2","firebrick3"/)

       resSci = res
       resSci@xyLineColor            = -1   ; transparent
       resSci@gsnXYFillOpacities     = .5
       shadecolor = (/"deepskyblue","gold","tomato"/)

       resSm = res
       resSm@xyMarkLineModes         = "Lines"
       resSm@xyLineThicknessF        = 15
       linecolor = (/"deepskyblue3","gold3","tomato3"/)

       ;======= Plot #1 :: Severity =======
       ; 1970-99, Synthetic CI
       resSci@gsnXYFillColors = shadecolor(0)
       resSci@tmYLLabelsOn    = True
       plot = gsn_csm_xy(wks,cdfsevsyntlohi,pctiles,resSci)
       P1_temp(:, 2) = cdfsevsyntlohi(0, :)
       P1_temp(:, 3) = cdfsevsyntlohi(1, :)


       ; 1970-99, Synthetic EM
       resSm@xyLineColor = linecolor(0)
       plotSm = gsn_csm_xy(wks,cdfsevsyntmed,pctiles,resSm)
       overlay(plot,plotSm)
       P1_temp(:, 1) = cdfsevsyntmed

       ; 1970-99, model CDF
       resM@xyMarkerColor = starcolor(0)
       plotm = gsn_csm_xy(wks,sevcdf,probax,resM)
       overlay(plot,plotm)

       P3_temp = new((/dimsizes(probax), 2/), "float", -9999.0)
       P3_temp(:, 0) = probax
       P3_temp(:, 1) = sevcdf



       ; future period detrended, Synthetic CI
       resSci@gsnXYFillColors = shadecolor(1)
       resSci@tmYLLabelsOn    = False
       plotf1 = gsn_csm_xy(wks,futcdfsevsyntlohi_dtr,pctiles,resSci)
       overlay(plot,plotf1)

       ; future period detrended, Synthetic EM
       resSm@xyLineColor = linecolor(1)
       plotfSm1 = gsn_csm_xy(wks,futcdfsevsyntmed_dtr,pctiles,resSm)
       overlay(plot,plotfSm1)

       ; future period detrended, model CDF
       resM@xyMarkerColor = starcolor(1)
       plotfm1 = gsn_csm_xy(wks,futsevcdf_dtr,futprobax_dtr,resM)
       overlay(plot,plotfm1)

       ; future period w/ trend, Synthetic CI
       resSci@gsnXYFillColors = shadecolor(2)
       plotf2 = gsn_csm_xy(wks,futcdfsevsyntlohi_tr,pctiles,resSci)
       overlay(plot,plotf2)

       P1_temp(:, 5) = futcdfsevsyntlohi_tr(0, :)
       P1_temp(:, 6) = futcdfsevsyntlohi_tr(1, :) 

       ; future period w/ trend, Synthetic EM
       resSm@xyLineColor = linecolor(2)
       plotfSm2 = gsn_csm_xy(wks,futcdfsevsyntmed_tr,pctiles,resSm)
       overlay(plot,plotfSm2)
       P1_temp(:, 4) = futcdfsevsyntmed_tr 

       ; future period w/ trend, model CDF
       resM@xyMarkerColor = starcolor(2)
       plotfm2 = gsn_csm_xy(wks,futsevcdf_tr,futprobax_tr,resM)
       overlay(plot,plotfm2)
       
       x_m = dimsizes(futprobax_tr)

       P2_temp = new((/x_m, 2/), "float", -9999.0) 

       P2_temp(:, 0) = futprobax_tr
       P2_temp(:, 1) = futsevcdf_tr
 
       fmtx = "7f15.8"
       opt1 = True
       opt1@fout = plotname+"severity_hist_fut_sythetic.txt"
       write_matrix(P1_temp, fmtx, opt1)
                

       fmtx = "2f15.8"
       opt1 = True
       opt1@fout = plotname+"severity_model_fut.txt"
       write_matrix(P2_temp, fmtx, opt1)


       opt1@fout = plotname+"severity_model_hist.txt"
       write_matrix(P3_temp, fmtx, opt1)


       delete(P2_temp)     
       delete(P3_temp)

       draw(plot)

       ;======= Plot #2 :: Duration =======
       
       resSci@vpXF                  = 0.59
       resSm@vpXF                   = resSci@vpXF
       resM@vpXF                    = resSci@vpXF

       ; 1970-99, Synthetic CI
       resSci@gsnXYFillColors = shadecolor(0)
       resSci@tmYLLabelsOn    = True
       plot = gsn_csm_xy(wks,cdfdursyntlohi,pctiles,resSci)

       ; 1970-99, Synthetic EM
       resSm@xyLineColor = linecolor(0)
       plotSm = gsn_csm_xy(wks,cdfdursyntmed,pctiles,resSm)
       overlay(plot,plotSm)

       ; 1970-99, model CDF
       resM@xyMarkerColor = starcolor(0)
       plotm = gsn_csm_xy(wks,durcdf,probax,resM)
       overlay(plot,plotm)

       ; future period detrended, Synthetic CI
       resSci@gsnXYFillColors = shadecolor(1)
       resSci@tmYLLabelsOn    = False
       plotf1 = gsn_csm_xy(wks,futcdfdursyntlohi_dtr,pctiles,resSci)
       overlay(plot,plotf1)

       ; future period detrended, Synthetic EM
       resSm@xyLineColor = linecolor(1)
       plotfSm1 = gsn_csm_xy(wks,futcdfdursyntmed_dtr,pctiles,resSm)
       overlay(plot,plotfSm1)

       ; future period detrended, model CDF
       resM@xyMarkerColor = starcolor(1)
       plotfm1 = gsn_csm_xy(wks,futdurcdf_dtr,futprobax_dtr,resM)
       overlay(plot,plotfm1)

       ; future period w/ trend, Synthetic CI
       resSci@gsnXYFillColors = shadecolor(2)
       plotf2 = gsn_csm_xy(wks,futcdfdursyntlohi_tr,pctiles,resSci)
       overlay(plot,plotf2)

       ; future period w/ trend, Synthetic EM
       resSm@xyLineColor = linecolor(2)
       plotfSm2 = gsn_csm_xy(wks,futcdfdursyntmed_tr,pctiles,resSm)
       overlay(plot,plotfSm2)

       ; future period w/ trend, model CDF
       resM@xyMarkerColor = starcolor(2)
       plotfm2 = gsn_csm_xy(wks,futdurcdf_tr,futprobax_tr,resM)
       overlay(plot,plotfm2)

       draw(plot)

       delete([/res,resM,resSci,resSm/])
 
;       drawNDCGrid(wks)

       frame(wks)  
    
       system("convert -trim "+plotname+"."+typPlt+" "+plotname+"."+typPlt)

    end do  ; end of extremes type loop
  end do  ; end of models loop
end do  ; end of regions loop
