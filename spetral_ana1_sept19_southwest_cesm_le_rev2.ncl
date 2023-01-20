;*****pec_1.ncl
;
; Concepts illustrated:
;   - Calculating and plotting spectra
;************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
;************************************************

;************************************************
; variable and file handling
;************************************************

 dirIn = "/glade/work/kumar34/hydroclimate_ncc/"

 nyrp = 30
 yrstr = (/1940,1970,2000,2030,2060/)
 yrend = yrstr + nyrp - 1

 region = (/"Southwest","Southeast","SouthernUSPlains","NorthernUSPlains","CanadianPlains"/)
 nrgs = dimsizes(region)

 model = (/"H2OSOI_CESM-LE", "SOILM_GFDLCM3-LE"/) 

 regN = 0
 modN = 0
 
 f1 = addfile(dirIn+model(modN)+"_1940-2089_30yr_TimeSeries_"+region(regN)+".nc", "r")
 print(f1)

 var_mod = f1->cesm_soilm(:, :, :)    ; cesm_soilm(periods, ens, time) 
 
 var_syn = f1->cesm_synt_soilm(:, :, :, :)   ; cesm_synt_soilm(periods, iter, ens, time)  
 
;************************

   d   = 0
   sm  = 1         ; periodogram
   pct = 0.10

  ;************************************************
  ; calculate mean spectrum spectrum and lag1 auto cor
  ;************************************************
  
  ; loop over each segment of length ntim
  
  
  
 ;  spcavg = new(5, floattoint(ntim/2), typeof(var_mod)) 
  
   ntim = dimsizes(var_mod&time) 

   nseg = dimsizes(var_mod&ens)
  
   dlen = ntim   

   spcavg = new((/5, floattoint(ntim/2)/), typeof(var_mod), var_mod@_FillValue)
  

   spcavg = 0.0
   
   r1zsum = spcavg(:, 0) 



  do i = 0, 4
      x = var_mod(i, :, :)
      
      do n = 0, nseg-1
         dof = specx_anal(x(n,:),d,sm,pct)         
         temp1 = dof@spcx
         temp2 = temp1*100/(dim_variance_n(x(n, :), 0)*dlen)
         spcavg(i, :) = spcavg(i, :) + temp2
         r1 = dof@xlag1
         r1zsum(i) = r1zsum(i)  + 0.5*(log((1+r1)/(1-r1))) ; sum the Fischer Z
         
      end do
   
  end do

   spcavg = spcavg/nseg
   r1z = r1zsum
   r1 :=  r1z

   do i = 0, 4
  
   r1z(i)  = r1zsum(i)/nseg                           ; average r1z
   r1(i)   = (exp(2*r1z(i))-1)/(exp(2*r1z(i))+1)            ; transform back

  end do

  print(sum(spcavg)/5)

  ;************************************************
  ; plotting
  ;************************************************
    wks  = gsn_open_wks("x11","spec")              ; Opens a ps file
  
    res = True
    ;res@tiMainString = "Mean Spectra: "+nseg+" segments, dof="+df                   ; title
    ;res@tiXAxisString = "Frequency (cycles/month)"  ; xaxis
     res@tiXAxisString = "Frequency (cycles/"+dlen+"year)"  ; xaxis
    res@tiYAxisString = "Variance"                  ; yaxis
  
    ;splt = specx_ci(df, 0.05, 0.95)                 ; confidence interval
  
    ;plot = gsn_csm_xy(wks, df@frq, splt ,res)      ; create plot


   res@gsnDraw = False
   res@gsnFrame = False

    res@vpHeightF= 0.4                    ; change aspect ratio of plot
    res@vpWidthF = 0.6

    res@trYMinF  = 2.0                   ; min value on y-axis
    res@trYMaxF  = 18.0                   ; max value on y-axis

    res@trXMinF  = 1.0                   ; min value on y-axiss
    res@trXMaxF  = 10.0                   ; max value on y-axis



     cmap = read_colormap_file("MPL_rainbow")
     idc = (/2, 20, 40, 93, 127/)
     clevs = cmap(idc,:)

     res@xyDashPattern  = 0                         ; Make curves all solid

     res@xyLineThicknesses = (/  4.0,   4.0, 4.0, 4.0, 4.0/)          ; make second line thicker
     res@xyLineColors      = clevs

     plot = gsn_csm_xy(wks, dof@frq*dlen, spcavg, res)      ; create plot


    df      = 2.0*nseg                           ; deg of freedom
                                                ; all segments
    df@spcx = spcavg(0, :)                             ; assign the mean spc
    df@frq  = dof@frq
    df@xlag1= r1(0)                            ; assign mean lag-1
    
    splt = specx_ci(df, 0.05, 0.95)                 ; confidence interval

    plot1 = gsn_csm_xy(wks, df@frq*dlen, splt ,res)      ; create plot

    overlay(plot, plot1)

    draw(plot)
    frame(wks)





