;***************************************************************************************
; Plots probability of exceedance of wet/dry events in CESM-LE and GFDLCM3-LE

; Takes the 95th percentile threshold event during each model's 1970-99 climatology
; (which has a 5% probability of exceedance), finds the equivalent severity/duration
; of that event in a future CDF of events for both trended and detrended soil moisture
; anomalies, and computes the probability of exceedance of such event in the future CDFs.

; Severity of events in all periods is scaled by the stddev of the model's event severity in 1970-99.
;***************************************************************************************

 url = "http://www.esrl.noaa.gov/psd/thredds/dodsC//Projects/SERDP_regimeshifts/data/"

 ;dirin = (/"/Projects/SERDP_regimeshifts/data/CESM-LE/lnd/postprocessed/regions_clean/CDFs/", \
 ;         "/Projects/SERDP_regimeshifts/data/GFDL-CM3-LE/SOILM/CDFs/"/)
 
 dirin = (/url+"CESM-LE/lnd/postprocessed/regions_clean/CDFs/", url+"GFDL-CM3-LE/SOILM/CDFs/"/)
 ;exit 
 ;dir = "/Projects/SERDP_regimeshifts/data/"
 dir = url
 dirmod = (/"CESM-LE/lnd/postprocessed/regions_clean/","GFDL-CM3-LE/SOILM/"/)

 dirext = (/"WetCDF_rndmSM0_v3/","DroughtCDF_rndmSM0_v3/"/)
 varmodnames = (/"H2OSOI_CESM-LE_","SOILM_GFDL-CM3-LE_"/)
 extfilename = (/"1940-2089_30yr_PluvialCDFs","1940-2089_30yr_DroughtCDFs"/)

 gcmvarname = (/"CESM-LE_H2OSOI_","GFDLCM3-LE_SOILM_"/)

 EvntType = (/"Pluvials","Droughts"/)
 datasrc = (/"Model","Synthetic"/)
 datatrnd = (/"Dtrend","wTrend"/)

 ; Pick the two regions to evaluate/plot
 ;region = (/"Southwest","Southeast"/) ; (/"CanadianPlains"/) ; (/"SouthernUSPlains","NorthernUSPlains"/) ; 
 ;region_title = (/"Southwest","Southeast"/) ; (/"Canadian Plains"/) ; (/"South US Plains","North US Plains"/) ; 

  ;region = (/"SouthernUSPlains","NorthernUSPlains"/) ;
  ;region_title = (/"South US Plains","North US Plains"/) ;


  region = (/"CanadianPlains"/) ; (/"SouthernUSPlains","NorthernUSPlains"/) ;
  region_title = (/"Canadian Plains"/) ; (/"South US Plains","North US Plains"/) ;


 ; Periods of analyses
 yrstr = (/1940,1970,2000,2030,2060/)
 yrend = (/1969,1999,2029,2059,2089/)
 ctryr = (/1955,1985,2015,2045,2075/)
 periodlabels = (/"1940-69","1970-99","2000-29","2030-59","2060-89"/)
 ; pick a period to compare against 1970-99
 p = 4

 ; Percentile rank IDs for given # of iterations
 ; (will be used to find median and 95% confidence interval in model and synthetic data)
 nit = 1000
 med =  toint(round(50./100.*(nit+1),0))
 pile025 = toint(round(2.5/100.*(nit+1),0))
 pile975 = toint(round(97.5/100.*(nit+1),0))

;***************************************************************************************

 ;dirout = "/Users/cdewes/Documents/SERDP_ClimateRegimeShifts/Analyses/plots/"
 dirout = "/glade/work/kumar34/hydroclimate_ncc/"
 typPlt = "png"
 typPlt@wkWidth  = 4200
 typPlt@wkHeight = 4200

 plotname = dirout+"SoilmExtremes_ProbExceed_1970-99_vs_"+periodlabels(p)+"_"+region(0)+"_v24"
; plotname = dirout+"SoilmExtremes_ProbExceed_1970-99_vs_"+periodlabels(p)+"_"+region(0)+"_v2"
 system("rm -f "+plotname+"."+typPlt)
 wks  = gsn_open_wks (typPlt,plotname)


;***************************************************************************************
 do r = 0,0  ; loop through regions (one panel for each)
   print_table([/region(r)/],"%s")

   ; Settings and resources for main plots
   res           = True
   res@gsnDraw   = False
   res@gsnFrame  = False

   res@vpHeightF = 0.2
   res@vpWidthF  = 0.3
   
   res@trXMinF   = 0
   res@trXMaxF   = 12
   res@trYMinF   = .4
   res@trYMaxF   = 90
   
   res@tmXBMode = "Explicit"
   xtickpos = (/2.5,3.5,4.5,7.5,8.5,9.5/)
;   xtickpos = (/1.25,2,2.75,4.75,5.5,6.25/)
;   xtickpos = (/1.25,2,2.75,5.25,6,6.75/)
   res@tmXBValues = xtickpos
   res@tmXBLabels = (/"H","F~B~dt","F~B~tr","H","F~B~dt","F~B~tr"/)
   res@tmXBLabelFontHeightF = 0.009
   res@tmXBLabelDeltaF = -.7
   res@tmXTOn = False
   res@tmYLMode = "Explicit"
   ytickpos = (/1,2,5,10,15,20,30,50,70,90/)
   res@tmYLValues = ytickpos
   res@tmYLLabels = ""+ytickpos
   res@xyYStyle = "Log"
   res@tmYLLabelFontHeightF = 0.009
   res@tmYLLabelDeltaF = -.7
   res@tmYMajorGrid                = True          ; implement y grid
   res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
   res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
   res@tmYROn = False

   res@tiYAxisFontHeightF = .015

   ; add "floating" bar chart resources to res
   res@gsnXYBarChart = True
   barwidth = .6
   res@gsnXYBarChartBarWidth = barwidth
   res@gsnXYBarChartOutlineThicknessF = 0.
   res@gsnYRefLineColor = "transparent"
   barcolor = (/"deepskyblue","gold","tomato","deepskyblue","gold","tomato"/)

   ; horizontal bar resources
   ebres = True
   ebres@gsMarkerIndex = NhlNewMarker(wks,">",34,0.,0.,2.,1.,0.)
   ebres@gsLineThicknessF = 10.

   ; "star" marker resources
   smres = True
   smres@gsMarkerIndex = NhlNewMarker(wks,"z",35,0.,0.,1.,1.8,0.)
   smres@gsLineThicknessF = 3.
   starcolor = (/"dodgerblue4","darkgoldenrod2","firebrick3","dodgerblue4","darkgoldenrod2","firebrick3"/)

   ; overlay "star" marker resources (for black outline)
   smres2 = True
   smres2@gsMarkerIndex = NhlNewMarker(wks,"]",19,0.,0.,1.,1.8,0.)
   smres2@gsLineThicknessF = 4.
   smres2@gsMarkerColor = "black"

   ; error bar resources
   ebres2 = True
   ebres2@gsMarkerIndex = NhlNewMarker(wks,">",34,0.,0.,1.5,.7,0.)
   ebres2@gsLineThicknessF = 5.
   
;***************************************************************************************

   do et = 0,1    ; type of event (pluvials or droughts)
     print_table([/EvntType(et)/],"%s")

     ; create arrays that will store the probability of exceedance values for: 
                 ; 1) CESM historical; 2) CESM future detrended; 3) CESM future w/ trend
                 ; 4) GFDL historical; 5) GFDL future detrended; 6) GFDL future w/ trend
     ; these will be used for plotting!

     sev_model := new((/6/),"float")   ; median of model probability of exceedance
     sev_modelhi := sev_model          ; 95% CI high value
     sev_modello := sev_model          ; 95% CI low value

     sev_synt := sev_model             ; median of synthetic probability of exceedance
     sev_synthi := sev_model           ; 95% CI high value
     sev_syntlo := sev_model           ; 95% CI low value

     dur_model := sev_model            ; same as above but for duration
     dur_modelhi := sev_model
     dur_modello := sev_model

     dur_synt := sev_model
     dur_synthi := sev_model
     dur_syntlo := sev_model

     do gcm = 0,1    ; CESM or GFDL
       print_table([/gcmvarname(gcm)/],"%s")

       ; read in model CDF (not the bootstrapped CDFs) and find the 95th-percentile threshold for 1970-99 events
       ; compute stddev of severity of the events for this period and use it as scaling reference
       ; sevsigma = stddev(cdf_severity)
       cdffile = dir+dirmod(gcm)+dirext(et)+varmodnames(gcm)+extfilename(et)+"_detrendedAnoms_"+region(r)+"_19701999.nc"
       cf = addfile(cdffile,"r")
       probax := cf->prob_axis
       sevcdf := cf->cdf_severity
       durcdf := cf->cdf_duration
       delete(cf)

       ;print(probax)
       ; print(durcdf)
       ;exit


       sevsigma = stddev(sevcdf)
       sevcdf = sevcdf / sevsigma

       prob95id = closest_val(95.,probax)
       ;prob95id = closest_val(95.5,probax)
       probexceed = 100. - probax(prob95id)
       sev95hist = sevcdf(prob95id)
       print_table([/"Severity of the 1970-99 95th-%ile "+EvntType(et)+" event is "+sprintf("%4.1f",sev95hist)+" standard deviations"/],"%s")
       dur95hist = durcdf(prob95id)
       print_table([/"Duration of the 1970-99 95th-%ile "+EvntType(et)+" event is "+sprintf("%4.1f",dur95hist)+" years"/],"%s")
       
       ;exit
       ; store in plotting variables
       if(gcm.eq.0) then
         sev_model(0) = probexceed
         dur_model(0) = probexceed
       else
         sev_model(3) = probexceed
         dur_model(3) = probexceed
       end if

       ; read in model (bootstrapped)/synthetic CDFs for 1970-99 events, scale severity by sevsigma and
       ; compute the probability of exceedance of the 95th-%ile threshold referenced above [sev95hist and dur95hist]

       do sms = 0,1    ; source of soilm anoms (model or synthetic)
         print_table([/datasrc(sms)/],"%s")

         cdffile = dir+dirmod(gcm)+"CDFs/"+gcmvarname(gcm)+region(r)+"_"+EvntType(et)+"CDFs_"+datasrc(sms)+"Anoms_Dtrend_19701999.nc"
         cf = addfile(cdffile,"r")
         nevents := cf->nevents
         probax := cf->prob_axis
         sevcdf := cf->cdf_severity
         sevcdf = sevcdf / sevsigma
         durcdf := cf->cdf_duration
         delete(cf)

         sev95probexc = new((/nit/),"float")
         dur95probexc = new((/2,nit/),"float",default_fillvalue("float"))
         do i = 0,nit-1
            ne = nevents(i)
            idsev95pile = closest_val(sev95hist,sevcdf(i,0:ne-1))
            sev95probexc(i) = 100. - probax(i,idsev95pile)
            iddur95pile := ind(durcdf(i,0:ne-1).eq.dur95hist)
            ;print(iddur95pile)
            ;iddur95pile := closest_val(dur95hist,durcdf(i,0:ne-1))
            ;print(iddur95pile)
            ;exit
           if(any(ismissing(iddur95pile))) then
             continue
           else
             probexcdur95hi = 100. - probax(i,iddur95pile(0))
             probexcdur95lo = 100. - probax(i,iddur95pile(dimsizes(iddur95pile)-1))
             dur95probexc(:,i) = (/probexcdur95lo,probexcdur95hi/)
           end if
         end do
       
         ; get the median and 95% confidence interval of probexceeds for 1970-99 synthetic events
         qsort(sev95probexc)
         sevprobexcmed = sev95probexc(med)
         sevprobexchi = sev95probexc(pile975)
         sevprobexclo = sev95probexc(pile025)
    
         dur95probexc1d = ndtooned(dur95probexc)
         qsort(dur95probexc1d)
         miss = num(ismissing(dur95probexc1d))
         nvals = dimsizes(dur95probexc1d) - miss
         md =  toint(round(50./100.*(nvals+1),0))
         p025 = toint(round(2.5/100.*(nvals+1),0))
         p975 = toint(round(97.5/100.*(nvals+1),0))
         durprobexcmed = dur95probexc1d(md)
         durprobexchi = dur95probexc1d(p975)
         durprobexclo = dur95probexc1d(p025)
       
         ; store in plotting variables
         if(gcm.eq.0 .and. sms.eq.0) then
           sev_modelhi(0) = sevprobexchi
           sev_modello(0) = sevprobexclo
           dur_modelhi(0) = durprobexchi
           dur_modello(0) = durprobexclo
         end if
         if(gcm.eq.0 .and. sms.eq.1) then
           sev_synt(0) = sevprobexcmed
           sev_synthi(0) = sevprobexchi
           sev_syntlo(0) = sevprobexclo
           dur_synt(0) = durprobexcmed
           dur_synthi(0) = durprobexchi
           dur_syntlo(0) = durprobexclo
         end if
         if(gcm.eq.1 .and. sms.eq.0) then
           sev_modelhi(3) = sevprobexchi
           sev_modello(3) = sevprobexclo
           dur_modelhi(3) = durprobexchi
           dur_modello(3) = durprobexclo
         end if
         if(gcm.eq.1 .and. sms.eq.1) then
           sev_synt(3) = sevprobexcmed
           sev_synthi(3) = sevprobexchi
           sev_syntlo(3) = sevprobexclo
           dur_synt(3) = durprobexcmed
           dur_synthi(3) = durprobexchi
           dur_syntlo(3) = durprobexclo
         end if
       end do
       
       ; read in model and synthetic CDFs for the future period's events (detrended and with trend),
       ; scale severity by sevsigma, and compute the probability of exceedance of the historical 95th-%ile event

       do sms = 0,1    ; source of soilm anoms (model or synthetic)
         print_table([/datasrc(sms)/],"%s")

         do tr = 0,1     ; detrended or with trend
           print_table([/datatrnd(tr)/],"%s")

           cdffile = dir+dirmod(gcm)+"CDFs/"+gcmvarname(gcm)+region(r)+"_"+EvntType(et)+"CDFs_"+datasrc(sms)+"Anoms_"+datatrnd(tr)+"_"+yrstr(p)+yrend(p)+".nc"
           cf = addfile(cdffile,"r")
           nevents := cf->nevents
           probax := cf->prob_axis
           sevcdf := cf->cdf_severity
           sevcdf = sevcdf / sevsigma
           durcdf := cf->cdf_duration
           delete(cf)

           sev95probexc = new((/nit/),"float")
           dur95probexc = new((/2,nit/),"float",default_fillvalue("float"))
           do i = 0,nit-1
             ne = nevents(i)
             idsev95pile = closest_val(sev95hist,sevcdf(i,0:ne-1))
             sev95probexc(i) = 100. - probax(i,idsev95pile)
             ;iddur95pile := ind(durcdf(i,0:ne-1).eq.dur95hist)
             iddur95pile := closest_val(dur95hist,durcdf(i,0:ne-1))
             if(any(ismissing(iddur95pile))) then
               continue
             else
               probexcdur95hi = 100. - probax(i,iddur95pile(0))
               probexcdur95lo = 100. - probax(i,iddur95pile(dimsizes(iddur95pile)-1))
               dur95probexc(:,i) = (/probexcdur95lo,probexcdur95hi/)
             end if
           end do

           ; get the median and 95% confidence interval of probexceeds for future model/synthetic events
           qsort(sev95probexc)
           sevprobexcmed = sev95probexc(med)
           sevprobexchi = sev95probexc(pile975)
           sevprobexclo = sev95probexc(pile025)

           dur95probexc1d = ndtooned(dur95probexc)
           qsort(dur95probexc1d)
           miss = num(ismissing(dur95probexc1d))
           nvals = dimsizes(dur95probexc1d) - miss
           if(nvals.lt.1) then
           ;if(nvals.lt.1) then
           ;if(nvals.le.100) then
             print(nvals)
             durprobexcmed = default_fillvalue("float")
             durprobexchi = default_fillvalue("float")
             durprobexclo = default_fillvalue("float")
             print_table([/"sample of future events too small to estimate new probability of exceedance"/],"%s")
             ;exit
           else
             md =  toint(round(50./100.*(nvals+1),0))
             p025 = toint(round(2.5/100.*(nvals+1),0))
             p975 = toint(round(97.5/100.*(nvals+1),0))
             durprobexcmed = dur95probexc1d(md)
             durprobexchi = dur95probexc1d(p975)
             durprobexclo = dur95probexc1d(p025)
           end if
           
           ; store in plotting variables
           if(gcm.eq.0 .and. sms.eq.0 .and. tr.eq.0) then
             sev_model(1) = sevprobexcmed
             sev_modelhi(1) = sevprobexchi
             sev_modello(1) = sevprobexclo
             dur_model(1) = durprobexcmed
             dur_modelhi(1) = durprobexchi
             dur_modello(1) = durprobexclo
           end if
           if(gcm.eq.0 .and. sms.eq.0 .and. tr.eq.1) then
             sev_model(2) = sevprobexcmed
             sev_modelhi(2) = sevprobexchi
             sev_modello(2) = sevprobexclo
             dur_model(2) = durprobexcmed
             dur_modelhi(2) = durprobexchi
             dur_modello(2) = durprobexclo
           end if
           if(gcm.eq.0 .and. sms.eq.1 .and. tr.eq.0) then
             sev_synt(1) = sevprobexcmed
             sev_synthi(1) = sevprobexchi
             sev_syntlo(1) = sevprobexclo
             dur_synt(1) = durprobexcmed
             dur_synthi(1) = durprobexchi
             dur_syntlo(1) = durprobexclo
           end if
           if(gcm.eq.0 .and. sms.eq.1 .and. tr.eq.1) then
             sev_synt(2) = sevprobexcmed
             sev_synthi(2) = sevprobexchi
             sev_syntlo(2) = sevprobexclo
             dur_synt(2) = durprobexcmed
             dur_synthi(2) = durprobexchi
             dur_syntlo(2) = durprobexclo
           end if
           if(gcm.eq.1 .and. sms.eq.0 .and. tr.eq.0) then
             sev_model(4) = sevprobexcmed
             sev_modelhi(4) = sevprobexchi
             sev_modello(4) = sevprobexclo
             dur_model(4) = durprobexcmed
             dur_modelhi(4) = durprobexchi
             dur_modello(4) = durprobexclo
           end if
           if(gcm.eq.1 .and. sms.eq.0 .and. tr.eq.1) then
             sev_model(5) = sevprobexcmed
             sev_modelhi(5) = sevprobexchi
             sev_modello(5) = sevprobexclo
             dur_model(5) = durprobexcmed
             dur_modelhi(5) = durprobexchi
             dur_modello(5) = durprobexclo
           end if
           if(gcm.eq.1 .and. sms.eq.1 .and. tr.eq.0) then
             sev_synt(4) = sevprobexcmed
             sev_synthi(4) = sevprobexchi
             sev_syntlo(4) = sevprobexclo
             dur_synt(4) = durprobexcmed
             dur_synthi(4) = durprobexchi
             dur_syntlo(4) = durprobexclo
           end if
           if(gcm.eq.1 .and. sms.eq.1 .and. tr.eq.1) then
             sev_synt(5) = sevprobexcmed
             sev_synthi(5) = sevprobexchi
             sev_syntlo(5) = sevprobexclo
             dur_synt(5) = durprobexcmed
             dur_synthi(5) = durprobexchi
             dur_syntlo(5) = durprobexclo
           end if
         end do
       end do
     end do
     
;***************************************************************************************
     ;=====================
     ; SEVERITY PLOT
     if(r.eq.0 .and. et.eq.0) then
       res@vpYF = 0.97
       res@tmXBOn = False
       res@trYReverse = False 
     end if
     if(r.eq.0 .and. et.eq.1) then
       res@vpYF = 0.77
       res@tmXBOn = True
       res@trYReverse = True 
     end if
     if(r.eq.1 .and. et.eq.0) then
       res@vpYF = 0.44
       res@tmXBOn = False
       res@trYReverse = False 
     end if
     if(r.eq.1 .and. et.eq.1) then
       res@vpYF = 0.24
       res@tmXBOn = True
       res@trYReverse = True 
     end if
     res@vpXF = 0.2
     res@tmYLOn = True
     res@tmYLLabelsOn = True
     res@tmYMajorGrid = True

     ; plot "floating" bar charts for synthetic 95% confidence range,
     ; marking the median value with a black horizontal bar
     plot = new((/2,dimsizes(xtickpos)/),graphic)
     plot_med = new((/dimsizes(xtickpos)/),graphic)
     do i = 0,dimsizes(xtickpos)-1
       if(ismissing(sev_synthi(i)).or.ismissing(sev_syntlo(i))) then
         continue
       else
       res@gsnXYBarChartColors = barcolor(i)
       res@gsnYRefLine = sev_synt(i)
       plot(0,i) = gsn_csm_xy(wks,xtickpos(i),sev_synthi(i),res)
       ;print(sev_synthi(i))
       res@tmXBOn = False
       res@tmYLOn = False
       res@tmYMajorGrid = False
       plot(1,i) = gsn_csm_xy(wks,xtickpos(i),sev_syntlo(i),res)
       
       plot_med(i) = gsn_add_polymarker(wks,plot(1,i),xtickpos(i),sev_synt(i),ebres)
       end if
     end do
   
     ; plot model values as "stars" with "error bars" to designate the 95% confidence range
     modelplot = new((/dimsizes(xtickpos)/),graphic)
     modelplot2 = modelplot
     cihi = modelplot
     cilo = modelplot
     vertbar = modelplot
     do i = 0,dimsizes(xtickpos)-1
       cihi(i) = gsn_add_polymarker(wks,plot(1,i),xtickpos(i),sev_modelhi(i),ebres2)
       cilo(i) = gsn_add_polymarker(wks,plot(1,i),xtickpos(i),sev_modello(i),ebres2)
       vertbar(i) = gsn_add_polyline(wks,plot(1,i),(/xtickpos(i),xtickpos(i)/),(/sev_modello(i),sev_modelhi(i)/),ebres2)
       smres@gsMarkerColor = starcolor(i)
       modelplot(i) = gsn_add_polymarker(wks,plot(1,i),xtickpos(i),sev_model(i),smres)
       modelplot2(i) = gsn_add_polymarker(wks,plot(1,i),xtickpos(i),sev_model(i),smres2)
       print(sev_model(i))
       print(sev_modelhi(i))
       print(sev_modello(i))
     end do

     draw(plot)

     ;=====================
     ; DURATION PLOT
   
     if(r.eq.0 .and. et.eq.0) then
       res@vpYF = 0.97
       res@tmXBOn = False
       res@trYReverse = False 
     end if
     if(r.eq.0 .and. et.eq.1) then
       res@vpYF = 0.77
       res@tmXBOn = True
       res@trYReverse = True 
     end if
     if(r.eq.1 .and. et.eq.0) then
       res@vpYF = 0.44
       res@tmXBOn = False
       res@trYReverse = False 
     end if
     if(r.eq.1 .and. et.eq.1) then
       res@vpYF = 0.24
       res@tmXBOn = True
       res@trYReverse = True 
     end if
     res@vpXF = 0.5
     res@tmYLOn = True
     res@tmYLLabelsOn = False
     res@tmYMajorGrid = True
     res@tmYLMajorOutwardLengthF = .005
   
     ; plot "floating" bar charts for synthetic 95% confidence range,
     ; marking the median value with a black horizontal bar
     plot := new((/2,dimsizes(xtickpos)/),graphic)
     plot_med = new((/dimsizes(xtickpos)/),graphic)
     do i = 0,dimsizes(xtickpos)-1
       if(ismissing(dur_synthi(i)).or.ismissing(dur_syntlo(i))) then
         continue
       else
       res@gsnXYBarChartColors = barcolor(i)
       res@gsnYRefLine = dur_synt(i)
       plot(0,i) = gsn_csm_xy(wks,xtickpos(i),dur_synthi(i),res)
       res@tmYLOn = False
       res@tmYMajorGrid = False
       plot(1,i) = gsn_csm_xy(wks,xtickpos(i),dur_syntlo(i),res)
       plot_med(i) = gsn_add_polymarker(wks,plot(1,i),xtickpos(i),dur_synt(i),ebres)
       end if
     end do
   
     ; plot model values as "stars" with "error bars" to designate the 95% confidence range
     modelplot = new((/dimsizes(xtickpos)/),graphic)
     modelplot2 = modelplot
     cihi = modelplot
     cilo = modelplot
     vertbar = modelplot
     do i = 0,dimsizes(xtickpos)-1
       cihi(i) = gsn_add_polymarker(wks,plot(1,i),xtickpos(i),dur_modelhi(i),ebres2)
       cilo(i) = gsn_add_polymarker(wks,plot(1,i),xtickpos(i),dur_modello(i),ebres2)
       if(.not.ismissing(dur_modelhi(i)) .and. .not.ismissing(dur_modello(i))) then
         vertbar(i) = gsn_add_polyline(wks,plot(1,i),(/xtickpos(i),xtickpos(i)/),(/dur_modello(i),dur_modelhi(i)/),ebres2)
       end if
       smres@gsMarkerColor = starcolor(i)
       modelplot(i) = gsn_add_polymarker(wks,plot(1,i),xtickpos(i),dur_model(i),smres)
       modelplot2(i) = gsn_add_polymarker(wks,plot(1,i),xtickpos(i),dur_model(i),smres2)
       print(starcolor(i))
       print(dur_model(i))
       print(dur_model(i))
       print(dur_modelhi(i))
       print(dur_modello(i))
   

     end do
   
     draw(plot)

   end do
 end do
 
   ;======= Plot Annotations =======

   txres = True

   txres@txFontHeightF = 0.011
   ypos = .54 ; top panel
   gsn_text_ndc(wks,"CESM-LE",.288,ypos,txres)
   gsn_text_ndc(wks,"GFDL-CM3-LE",.41,ypos,txres)
   gsn_text_ndc(wks,"CESM-LE",.588,ypos,txres)
   gsn_text_ndc(wks,"GFDL-CM3-LE",.71,ypos,txres)
   ypos = .01 ; bottom panel
   gsn_text_ndc(wks,"CESM-LE",.288,ypos,txres)
   gsn_text_ndc(wks,"GFDL-CM3-LE",.41,ypos,txres)
   gsn_text_ndc(wks,"CESM-LE",.588,ypos,txres)
   gsn_text_ndc(wks,"GFDL-CM3-LE",.71,ypos,txres)

   txres@txFontHeightF = 0.018
   ypos = .985
   gsn_text_ndc(wks,"Severity",.35,ypos,txres)
   gsn_text_ndc(wks,"Duration",.65,ypos,txres)

   txres@txAngleF = 90.
   gsn_text_ndc(wks,region_title(0),.09,.77,txres) ; top panel
   gsn_text_ndc(wks,region_title(1),.09,.24,txres) ; bottom panel
   txres@txFontHeightF = 0.015
   gsn_text_ndc(wks,"Probability of Exceedance (%)",.125,.77,txres) ; top panel
   gsn_text_ndc(wks,"Probability of Exceedance (%)",.125,.24,txres) ; bottom panel

   txres@txFontHeightF = 0.015
   gsn_text_ndc(wks,"wet extremes",.155,.90,txres) ; top panel
   gsn_text_ndc(wks,"dry extremes",.155,.64,txres)
   gsn_text_ndc(wks,"wet extremes",.155,.37,txres) ; bottom panel
   gsn_text_ndc(wks,"dry extremes",.155,.11,txres)

   ; draw arrows indicating direction of probability scale
   polyres = True
   polyres@gsMarkerColor = "black"
   uparrow = NhlNewMarker(wks,"-",34,0.,0.,2.,2,0.)
   polyres@gsMarkerIndex = uparrow
   gsn_polymarker_ndc(wks,.155,.80,polyres) ; top panel
   gsn_polymarker_ndc(wks,.155,.27,polyres) ; bottom panel
   downarrow = NhlNewMarker(wks,"/",34,0.,0.,2.,2,0.)
   polyres@gsMarkerIndex = downarrow
   gsn_polymarker_ndc(wks,.155,.74,polyres) ; top panel
   gsn_polymarker_ndc(wks,.155,.21,polyres) ; bottom panel
   
   ;======= Plot Legend =======
   
   lglabels = (/"historical, 1970-99 clim","future detrended, "+periodlabels(p)+" clim", \
                "future w/ trend, "+periodlabels(p)+" w/ 1970-99 clim", \
                "Model, historical","Model, future detrended","Model, future w/ trend", \
                "Synthetic, historical","Synthetic, future detrended","Synthetic, future w/ trend"/)
   sqr_mid = NhlNewMarker(wks,"y",35,0.,0.,1.,2,0.)
   str_mid = NhlNewMarker(wks,"z",35,0.,0.,1.,2,0.)
   H_mid = NhlNewMarker(wks,"H",21,0.,0.,1.,1.5,0.)
   Fdt_mid = NhlNewMarker(wks,"F dt",21,0.,0.,1.,1.5,0.)
   Ftr_mid = NhlNewMarker(wks,"F tr",21,0.,0.,1.,1.5,0.)
   mkids = (/H_mid,Fdt_mid,Ftr_mid,str_mid,str_mid,str_mid,sqr_mid,sqr_mid,sqr_mid/)
   mkcolors = (/"black","black","black",starcolor(0),starcolor(1),starcolor(2),barcolor(0),barcolor(1),barcolor(2)/)
   xpos := (/.10,.10,.10,.40,.40,.40,.60,.60,.60/)
   ypos := (/.51,.49,.47,.51,.49,.47,.51,.49,.47/)

   txres@txFontHeightF = 0.01
   txres@txAngleF = 0.
   txres@txJust = "BottomLeft"
   do m = 0,dimsizes(lglabels)-1
     polyres@gsMarkerIndex = mkids(m)
     polyres@gsMarkerColor = mkcolors(m)
     gsn_polymarker_ndc(wks,xpos(m),ypos(m),polyres)
     gsn_text_ndc(wks,lglabels(m),xpos(m)+.02,ypos(m)-.0075,txres)
   end do

   ; draw a box around legend items
   xboxcoords = (/.08,.80,.80,.08,.08/)
   yboxcoords = (/.52,.52,.46,.46,.52/)
   lnres = True
   lnres@gsLineColor = "black"
   lnres@gsLineThicknessF = 3.
   gsn_polyline_ndc(wks,xboxcoords,yboxcoords,lnres)
   
   ;===========================

;   drawNDCGrid(wks)

   frame(wks)      
   system("convert -trim "+plotname+"."+typPlt+" "+plotname+"."+typPlt)

