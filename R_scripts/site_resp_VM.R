site_resp_VM <- function(period,PGAr,vs30,Z_LCG){
# Function that estimates the site response for the Valley of Mexico conditioned
# on PGAr, Vs30, and Z_LCG.
# Ref: Contreras V, Stewart JP, Mayoral JM, Pérez-Campos X. (2025). Valley of 
# Mexico Site Response from Non-Reference Site Approach, Earthquake Spectra
# (accepted)
#
# INPUT
# period (s): Oscillator period (PGV = -1, PGA = 0, or period between 0.01s-10s)
# PGAr (g)  : Peak Ground Acceleration at reference rock site (Vs30 = 760m/s)
#             (originally computed using NGA-Sub Parker et al. 2022 GMM)
# vs30 (m/s): Time-averaged shear wave velocity in the upper 30 m of a site
# Z_LCG (m) : Sediment depth (measured to the Lower Coarse Grained layer)
#
# OUTPUT
# Fs_VM     : Total site response (nonlinear combined model) in natural log units
# Amp       : Total site amplification factor (nonlinear combined model)
# Fs_NL     : Nonlinear part of the site response in natural log units
# Fs_L      : Linear part of the site response in natural log units
# Fs_Z      : Sediment depth adjustment in natural log units
# dZ (m)    : Differential depth
#
# Victor Contreras, Universidad Diego Portales
# 2024/12/14 (latest version)  

# Linear response
  coeffs <- read.csv("VM_site_response_model_coeffs.csv", header =  TRUE, skip = 0)
  IM <- coeffs$IM
  s1 <- interp1(IM,coeffs$s_1vm,period)
  s2 <- interp1(IM,coeffs$s_2vm,period)
  V1 <- interp1(IM,coeffs$V_1vm,period)
  f_VM <- interp1(IM,coeffs$f_VM,period)
  Vref <- 760
  if(vs30 <= V1){
    Fs_L <- s1*log(vs30/V1) + s2*log(V1/Vref) + f_VM
  }else if(vs30 > V1){
    Fs_L <- s2*log(vs30/Vref) + f_VM
  }

# Nonlinear response
  f3 <- 0.01
  if(vs30>250 & vs30<=760){ # Zone I
    f2 <- interp1(IM,coeffs$f_2ZI,period)
  }else if(vs30>140 & vs30<=250){ # Zone II
    f2 <- interp1(IM,coeffs$f_2ZII,period)
  }else if(vs30>90 & vs30<=140){ # Zone IIIa
    f2 <- interp1(IM,coeffs$f_2ZIIIa,period)
  }else if(vs30>=60 & vs30<=90){  # Zone IIIb
    f2 <- interp1(IM,coeffs$f_2ZIIIb,period)
  }else if(vs30>=40 & vs30<60){ # Zones IIIc-IIId
    f2 <- mean(c(interp1(IM,coeffs$f_2ZIIIc,period),interp1(IM,coeffs$f_2ZIIId,period)))
  }
  Fs_NL <- f2*log((PGAr + f3)/f3)
  
# Sediment depth correction
# Differential depth
  if(Z_LCG==0 | vs30>=250){
    dZ <- 0
  }else{
    dZ <-  Z_LCG - exp(8.684 - 1.22*log(vs30))
  }
# d1, d2 smoothing             IM       d1      d2
  data_smoothing <-   matrix(c(0.65, -0.0155, -0.0117,
                              0.90, -0.0140, -0.0100,
                              1.10, -0.0130, -0.0085,
                              1.30, -0.0125, -0.0080,
                              1.75, -0.0130, -0.0100,
                              2.25, -0.0360, -0.0460,
                              2.75, -0.0400, -0.0470,
                              3.50, -0.0045, -0.0100,
                              4.40, -0.0000,  0.0305),nrow=9,ncol=3,byrow=T)
  IM_sm <- c(IM,data_smoothing[,1])
  d1_sm <- c(coeffs$d_1,data_smoothing[,2])
  d2_sm <- c(coeffs$d_2,data_smoothing[,3])
  IM_sm <- sort(IM_sm,index.return=T)
  ind <- IM_sm$ix
  IM_sm <- IM_sm$x
  d1_sm <- d1_sm[ind];
  d2_sm <- d2_sm[ind];
  d1 <- interp1(IM_sm,d1_sm,period);
  d2 <- interp1(IM_sm,d2_sm,period);
  
  if(dZ < (-5)){
    Fs_Z <- d1*(dZ + 5)
  }else if(dZ > 5){
    Fs_Z <- d2*(dZ - 5)
  }else{
    Fs_Z <- 0
  }

# Total site response and amplification factor
  Fs_VM <- Fs_NL + Fs_L + Fs_Z
  Amp <- exp(Fs_VM)
  return(c(Fs_VM,Amp,Fs_NL,Fs_L,Fs_Z,dZ))
}