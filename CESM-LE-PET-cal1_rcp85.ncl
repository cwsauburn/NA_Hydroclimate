load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;******README for revision*************************
;***In this reviesion for mod 9 only 1st ensemble is used*************
;*************done reading with GLDAS data****************
;******************* adding GPCP precip data ***************************

 fm = addfile("/glade/p/uola0001/prj_data/surf_data/CLM4.0CN.h0.fixed.nc", "r")
 lm = fm->landmask
 area = fm->area
 landfrac = fm->landfrac
 area@_FillValue = 1.0*10^36
 landfrac@_FillValue = 1.0*10^36
 arwt = area
 arwt = area*landfrac/100
 lat = fm->lat
 lon = fm->lon

;****************************************************************************

 nEns = 40
 
 ens_ind = new((/nEns/), "integer", -99)
 
 ens_ind(0:34) = ispan(1, 35, 1)
 ens_ind(35:39) = ispan(101, 105, 1)
 
 varTrend = new((/2, nEns, dimsizes(lat), dimsizes(lon)/), "float", -9999.0)
 
 varTrend!0 = "MKZ_P_TSA"
 varTrend!1 = "ensemble"
 varTrend!2 = "lat"
 varTrend!3 = "lon"
 
 varTrend&lat = lat
 varTrend&lon = lon

 ;*******************************************************************************
 
  inDirT = "/glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/TREFHT/"
  inDirP = "/glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/PSL/"
  inDirU = "/glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/U10/"
  inDirRH = "/glade/p/cesmLE/CESM-CAM5-BGC-LE/lnd/proc/tseries/monthly/RH2M/"
  inDirS = "/glade/p/cesmLE/CESM-CAM5-BGC-LE/lnd/proc/tseries/monthly/FSA/"
  inDirL = "/glade/p/cesmLE/CESM-CAM5-BGC-LE/lnd/proc/tseries/monthly/FIRA/"
  inDirSM = "/glade/p/cesmLE/CESM-CAM5-BGC-LE/lnd/proc/tseries/monthly/H2OSOI/"
   
  run = (/"B20TRC5CNBDRD","BRCP85C5CNBDRD"/)
 
  outDir = "/glade/scratch/kumar34/CESM-LE"
   
  nEns = 40
  ens_ind = new((/nEns/), "integer", -99)
  ens_ind(0:34) = ispan(1, 35, 1)
  ens_ind(35:39) = ispan(101, 105, 1)
  
  
   
  do i = 31, 39
      file1 = systemfunc("ls "+inDirT+"b.e11."+run(1)+".f09_g16."+sprinti("%0.3i", ens_ind(i))+".cam.h0.TREFHT.*.nc")
      ft1 := addfile(file1(0),"r")
      date1 := ft1->date
      sti1 = ind(date1.eq.20060201)
      endi1 = ind(date1.eq.20810101)
      TREFHT := ft1->TREFHT(sti1:endi1, :, :)
     
      ;****unit converstion to degree C
      TREFHT = TREFHT - 273.15      ; in degree C
     
      ES = TREFHT
     
      ES = 611.0 * exp((17.27*TREFHT)/(237.3+TREFHT))   ; N/m2
     
      grdES = ES 
     
      grdES = (/(4098.0*ES/((237.3+TREFHT)^2))/)          ; unit N/m2 per degree C
     
      file2 = systemfunc("ls "+inDirRH+"b.e11."+run(1)+".f09_g16."+sprinti("%0.3i", ens_ind(i))+".clm2.h0.RH2M.*.nc")
      ft2 := addfile(file2(0),"r")
      date2 := ft2->mcdate
      sti2 = ind(date2.eq.20060201)
      endi2 = ind(date2.eq.20810101)
      RH2M := ft2->RH2M(sti2:endi2, :, :)
      
      VP = RH2M
      VP = (/(RH2M * ES /100.0)/)  
      
      VPD = RH2M
      VPD = (/ES - VP/)     
      
      file3 = systemfunc("ls "+inDirP+"b.e11."+run(1)+".f09_g16."+sprinti("%0.3i", ens_ind(i))+".cam.h0.PSL.*.nc")
      ft3 := addfile(file3(0),"r")
      date3 := ft3->date
      sti3 = ind(date3.eq.20060201)
      endi3 = ind(date3.eq.20810101)
      PSL := ft3->PSL(sti3:endi3, :, :)
      
      
      CP = 1.00464 * 10^3    ; J kg-1 K-1
      
      LV = 2.501 * 10^6     ; J kg-1
      
      PSY = PSL
      
      PSY = (/ (PSL * (CP/(0.622*LV)))/)    ; N m-2 K-1
      
      conv_fact = (24.0 * 60.0* 60.0 / LV)   ; watts/m2 to mm/day
      
      ;************
      
      file4 = systemfunc("ls "+inDirS+"b.e11."+run(1)+".f09_g16."+sprinti("%0.3i", ens_ind(i))+".clm2.h0.FSA.*.nc")
      ft4 := addfile(file4(0),"r")
      date4 := ft4->mcdate
      sti4 = ind(date4.eq.20060201)
      endi4 = ind(date4.eq.20810101)
      FSA := ft4->FSA(sti4:endi4, :, :)
     
      file5 = systemfunc("ls "+inDirL+"b.e11."+run(1)+".f09_g16."+sprinti("%0.3i", ens_ind(i))+".clm2.h0.FIRA.*.nc")
      ft5 := addfile(file5(0),"r")
      date5 := ft5->mcdate
      sti5 = ind(date5.eq.20060201)
      endi5 = ind(date5.eq.20810101)
      FIRA := ft5->FIRA(sti5:endi5, :, :)
      
      
      RN = FSA
      RN = FSA - FIRA    ; FIRA is postive upward (see page 63 in CLM4.5 tech doc)
     
     ;***********************************************************************************************************
     file6 = systemfunc("ls "+inDirU+"b.e11."+run(1)+".f09_g16."+sprinti("%0.3i", ens_ind(i))+".cam.h0.U10.*.nc")
     ft6 := addfile(file6(0),"r")
     date6 := ft6->date
     sti6 = ind(date6.eq.20060201)
     endi6 = ind(date6.eq.20810101)
     U10 := ft6->U10(sti6:endi6, :, :)
      
     ;*************************************************************************************************************
     
     U2 = U10
     U2 = (/U10 * (log(128.0)/log(661.3)) /)
     
     ;************************************************************************************************************
     
     
     PET = RN
     PET@long_name = "Potential ET"
     PET@units     = "mm/day"
     
     t1 = grdES/(grdES+PSY)
     
     
     t2 = PSY/(grdES+PSY)
     Udum = U2
     Udum = 1.0
     t3 = 6430 * (Udum + 0.536 * U2)
     
     
     PET = (/(t1*RN*conv_fact)+(t2*t3*VPD/LV)/)           ; PE = (grdES/(grdES+PSY))* RNET + (PSY/(grdES+PSY)) * (6430 (1+0.536U2)*VPD/LV)  ...Donohue et al., 2012
     
     system("/bin/rm -f "+outDir+"/b.e11."+run(1)+".f09_g16."+sprinti("%0.3i", ens_ind(i))+".clm2.h0.PET.200601-208012.nc")
     fout1 = addfile(outDir+"/b.e11."+run(1)+".f09_g16."+sprinti("%0.3i", ens_ind(i))+".clm2.h0.PET.200601-208012.nc", "c")    
     fout1->PET = PET
    
     
 end do

 ; date = array_append_record(date1(sti1:endi1),date2(sti2:endi2),0)
 
 ; time1 = ft1->time(sti1:endi1)
 ; time2 = ft2->time(sti2:endi2)
 ; time2 = time_to_newtime(time2,time1@units)
 ; time = array_append_record(time1,time2,0)

 ;TREFHT_full&time = time
 ;TREFHT_full!0 = "ens"
 ;TREFHT_full&ens = ens_ind

 ; printVarSummary(TREFHT_full) 
 
 
