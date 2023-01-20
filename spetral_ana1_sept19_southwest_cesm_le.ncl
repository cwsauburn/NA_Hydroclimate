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

 varP_mod = f1->cesm_soilm(1, :, :)    ; cesm_soilm(periods, ens, time) 
 varF_mod = f1->cesm_soilm(4, :, :)  
 
 varP_syn = f1->cesm_synt_soilm(1, :, :, :)   ; cesm_synt_soilm(periods, iter, ens, time)  
 varF_syn = f1->cesm_synt_soilm(4, :, :, :)

;************************

   d   = 0
   sm  = 1         ; periodogram
   pct = 0.10

  ;************************************************
  ; calculate mean spectrum spectrum and lag1 auto cor
  ;************************************************
  
  ; loop over each segment of length ntim
  
   x = varF_mod
   ntim = dimsizes(varP_mod&time) 

   x1d_var = dim_variance_n(ndtooned(x), 0)

   x_var = dim_variance_n(x, 1)


   dlen = ntim

   nseg = dimsizes(varP_mod&ens)

   spcavg = new ( floattoint(ntim/2), typeof(x))
   spcavg = 0.0
  
   r1zsum = 0.0

   n = 0
   dof1 = specx_anal(x(n,:),d,sm,pct)
   printVarSummary(dof1) 

   ;exit
 
   do n=0,nseg-1
      dof    = specx_anal(x(n,:),d,sm,pct)      ; current segment spc
      ;printVarSummary(dof)
      ;exit
      temp1 = dof@spcx
      print("i am printing")
      print(sum(temp1))
      temp1 = temp1 *100.0/ (x_var(n)*(dlen-1))
      ;temp2 = sum(temp1)
      ;temp3 = temp1*100.0/temp2
      print(sum(temp1))
      ;spcavg = spcavg + dof@spcx                ; sum spc of each segment
      spcavg = spcavg + temp1
      r1     = dof@xlag1                        ; extract segment lag-1
      r1zsum = r1zsum  + 0.5*(log((1+r1)/(1-r1))) ; sum the Fischer Z
   end do
  
   r1z  = r1zsum/nseg                           ; average r1z
   r1   = (exp(2*r1z)-1)/(exp(2*r1z)+1)            ; transform back
    
   ;spcavg = spcavg /x1d_var

   print(sum(spcavg))                                            ; this is the mean r1
   spcavg  = spcavg/nseg
   print(sum(spcavg))                        ; average spectrum
   printVarSummary(spcavg)
  
  ;************************************************
  ; Assign mean spectrum to data object
  ;************************************************
  
   df      = 2.0*nseg                           ; deg of freedom
                                                ; all segments
   df@spcx = spcavg                             ; assign the mean spc
   df@frq  = dof@frq
   df@xlag1= r1                            ; assign mean lag-1
  
  ;************************************************
  ; plotting
  ;************************************************
    wks  = gsn_open_wks("x11","spec")              ; Opens a ps file
  
    res = True
    res@tiMainString = "Mean Spectra: "+nseg+" segments, dof="+df                   ; title
    res@tiXAxisString = "Frequency (cycles/"+dlen+"year)"  ; xaxis
    res@tiYAxisString = "Variance"                  ; yaxis
  
    splt = specx_ci(df, 0.05, 0.95)                 ; confidence interval
  
    plot = gsn_csm_xy(wks, df@frq*dlen, splt ,res)      ; create plot

    ;plot = gsn_csm_xy(wks, df@frq*dlen, spcavg, res)      ; create plot





