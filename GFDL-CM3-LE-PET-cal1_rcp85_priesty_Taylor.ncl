load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;******README for revision*************************
;***In this reviesion for mod 9 only 1st ensemble is used*************
;*************done reading with GLDAS data****************
;******************* adding GPCP precip data ***************************

 ;;fm = addfile("/glade/p/uola0001/prj_data/surf_data/CLM4.0CN.h0.fixed.nc", "r")
 ;;lm = fm->landmask
 ;;area = fm->area
 ;;landfrac = fm->landfrac
 ;;area@_FillValue = 1.0*10^36
 ;;landfrac@_FillValue = 1.0*10^36
 ;;arwt = area
 ;;arwt = area*landfrac/100
 ;;lat = fm->lat
 ;;lon = fm->lon

;************************************
;****************************************

 fm = addfile("/glade/collections/cdg/data/CLIVAR_LE/gfdl_cm3_lens/Amon/tas/tas_Amon_GFDL-CM3_historical_rcp85_r1i1p1_192001-210012.nc", "r")
 lat = fm->lat
 lon = fm->lon




;************************************
 nEns = 20
 
 ens_ind = new((/nEns/), "integer", -99)
 
 ens_ind(0:19) = ispan(1, 20, 1)
 
   
  inDirT = "/glade/collections/cdg/data/CLIVAR_LE/gfdl_cm3_lens/Amon/tas/"
  inDirP = "/glade/collections/cdg/data/CLIVAR_LE/gfdl_cm3_lens/Amon/psl/"
  inDirU = "/glade/collections/cdg/data/CLIVAR_LE/gfdl_cm3_lens/Amon/ua/"
  inDirV = "/glade/collections/cdg/data/CLIVAR_LE/gfdl_cm3_lens/Amon/va/"
  inDirRH = "/glade/collections/cdg/data/CLIVAR_LE/gfdl_cm3_lens/Amon/hus/"
  inDirS_D = "/glade/collections/cdg/data/CLIVAR_LE/gfdl_cm3_lens/Amon/rsds/"
  inDirS_U = "/glade/collections/cdg/data/CLIVAR_LE/gfdl_cm3_lens/Amon/rsus/"
  inDirL   = "/glade/collections/cdg/data/CLIVAR_LE/gfdl_cm3_lens/Amon/flns/"


    
  outDir = "/glade/work/kumar34/GFDL-CM3LE/PET/"







  nEns = 20
  ens_ind = new((/nEns/), "integer", -99)
  ens_ind(0:19) = ispan(1, 20, 1)
  ;ens_ind(35:39) = ispan(101, 105, 1)
  
 ;i = 0 
   
  do i = 0, 19
      ft1 = addfile(inDirT+"tas_Amon_GFDL-CM3_historical_rcp85_r"+ens_ind(i)+"i1p1_192001-210012.nc", "r")
      TREFHT := ft1->tas
    
      printVarSummary(TREFHT)

      ;exit
 
      ;****unit converstion to degree C
      TREFHT = TREFHT - 273.15      ; in degree C
     
      ES = TREFHT
     
      ES = 611.0 * exp((17.27*TREFHT)/(237.3+TREFHT))   ; N/m2
     
      grdES = ES 
     
      grdES = (/(4098.0*ES/((237.3+TREFHT)^2))/)          ; unit N/m2 per degree C
     
      ft2 = addfile(inDirRH+"hus_Amon_GFDL-CM3_historical_rcp85_r"+ens_ind(i)+"i1p1_192001-210012.nc", "r")

      hus := ft2->hus(:, 0, :, :)
      printVarSummary(hus)
      print(avg(hus))
     
      ft3 = addfile(inDirP+"psl_Amon_GFDL-CM3_historical_rcp85_r"+ens_ind(i)+"i1p1_192001-210012.nc", "r")
      PSL = ft3->psl
      PSL = (/PSL * 10^2/)    ; hPa to Pa
      printVarSummary(PSL)


      hus_sat = hus
      hus_sat := (/(0.622 * ES / PSL)/)
      RH2M = hus
      RH2M = hus * (/(100.0 / hus_sat)/)
      print(avg(RH2M))



      VP = RH2M
      VP = (/(RH2M * ES /100.0)/)  
      
      VPD = RH2M
      VPD = (/ES - VP/)     
      
      
      CP = 1.00464 * 10^3    ; J kg-1 K-1
      
      LV = 2.501 * 10^6     ; J kg-1
      
      PSY = PSL
      
      PSY = (/ (PSL * (CP/(0.622*LV)))/)    ; N m-2 K-1
      
      conv_fact = (24.0 * 60.0* 60.0 / LV)   ; watts/m2 to mm/day
      
      ;************
      
      
      ft4 = addfile(inDirS_D+"rsds_Amon_GFDL-CM3_historical_rcp85_r"+ens_ind(i)+"i1p1_192001-210012.nc", "r")
      rsds = ft4->rsds


      ft5 = addfile(inDirS_U+"rsus_Amon_GFDL-CM3_historical_rcp85_r"+ens_ind(i)+"i1p1_192001-210012.nc", "r")
      rsus = ft5->rsus

      ft6 = addfile(inDirL+"flns_Amon_GFDL-CM3_historical_rcp85_r"+ens_ind(i)+"i1p1_192001-210012.nc", "r")
      flns = ft6->flns


      
      RN = rsds
      RN = rsds - rsus + flns    ; flns is postive downward (see page 63 in CLM4.5 tech doc)
     
      
      ft7 = addfile(inDirU+"ua_Amon_GFDL-CM3_historical_rcp85_r"+ens_ind(i)+"i1p1_192001-210012.nc", "r")
      ua = ft7->ua(:, 0, :, :)

      ft8 = addfile(inDirV+"va_Amon_GFDL-CM3_historical_rcp85_r"+ens_ind(i)+"i1p1_192001-210012.nc", "r")
      va = ft8->va(:, 0, :, :)

      U10 = ua
      U10 = (/(ua^2 + va^2)^0.5/)
   


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
     
     
     ;PET = (/(t1*RN*conv_fact)+(t2*t3*VPD/LV)/)           ; PE = (grdES/(grdES+PSY))* RNET + (PSY/(grdES+PSY)) * (6430 (1+0.536U2)*VPD/LV)  ...Donohue et al., 2012
     
     PET = (/t1 * RN*conv_fact * 1.3/)                          ; Dewes et al., 2017 


     system("/bin/rm -f "+outDir+"/PET_pt_Amon_GFDL-CM3_historical_rcp85_r"+ens_ind(i)+"i1p1_192001-210012.nc")
     fout1 = addfile(outDir+"/PET_pt_Amon_GFDL-CM3_historical_rcp85_r"+ens_ind(i)+"i1p1_192001-210012.nc", "c")    
     fout1->PET = PET
    
     
 end do

