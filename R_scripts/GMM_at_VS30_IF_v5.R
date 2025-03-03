#Grace Parker
#Modified February 26, 2020, to expand comments
#Modified March 25, 2020, to take coefficients from "Table_E1_Interface_Coefficients.csv"
#Modified December 21, 2024 to incorporate CAM updates concerning backarc
#anelastic attenuation and site response in the Valley of Mexico (VM) (Victor Contreras)

# Input Parameters --------------------------------------------------------

#Event type: 0 == interface, 1 == slab

#region corresponds to options in the DatabaseRegion column of the flatfile, plus global. Must be a string. If no matches, default will be global model:
  # "global", "Alaska", "Cascadia", "CAM", "Japan", "SA" or "Taiwan"

#Saturation Region corresponds to regions defined by R. Archuleta and C. Ji:
  # "global", "Aleutian","Alaska","Cascadia","CAM_S", "CAM_N", "Japan_Pac","Japan_Phi","SA_N","SA_S", "Taiwan_W","Taiwan_E"

# Rrup is number in kilometers

# period can be: (-1,0,0.01,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,2.5,3,4,5,7.5,10) 
  # where -1 == PGV and 0 == PGA

#VS30 in units m/s

#Z2.5 in units m. Only used if DatabaseRegion == "Japan" or "Cascadia". Can also specify "default" to get no basin term

#basin is only used if DatabaseRegion == "Cascadia". Value can be 0, 1, or 2, where 0 == having an estimate of Z2.5 outside mapped basin, 1 == Seattle basin, and 0 == other mapped basin (Tacoma, Everett, Georgia, etc.)

# CAM only, Dec/24 updates:
# Rrup_b is the backarc portion of the rupture distance in km
# File "CAM_backarc_attenuation_coeffs.csv" with coefficients must be in the active directory

# VM_flag is a flag indicating if the site is located in the VM (1 = inside VM, 0 = outside VM) 

# Z_LCG is the sediment depth in m (measured to the Lower Coarse Grained layer in VM)

# Other pertinent information ---------------------------------------------
#Coefficient files must be in the active working directory
# "GMM_at_VS30_IF_v5.R" calls function "GMM_at_760_IF_v5.R" to compute PGAr in the nonlinear site term. This function must be in the R environment else an error will occur.
# The output is the desired median model prediction in LN units
# Take the exponential to get PGA, PSA in g or the PGV in cm/s

GMM_at_VS30_IF_v5 <- function(event.type,region,saturation.region,Rrup,M,period,VS30,Z2.5,basin,Rrup_b,VM_flag,Z_LCG){
 
   if(event.type == 1 | event.type == 5){
    stop("This function is only for IF")
  }
  
# Import Master Coefficient Table -----------------------------------------

  coefficients <- read.csv("Table_S1_Interface_Coefficients.csv", header =  T)
 
  #Isolate desired period
   coefficients.T <- subset(coefficients,Period..s. == period)

  #Define Mb
  if(saturation.region == "global"){
    Mb <- 7.9
  }else{
    Interface.saturation.regions <- data.frame(SBZ = c("Aleutian","Alaska","-999","Cascadia","CAM_S", "CAM_N", "Japan_Pac","Japan_Phi", "New_Zealand_N","New_Zealand_S","SA_N","SA_S", "Taiwan_W","Taiwan_E"),
                                               Mb = c(8,8.6,NA,7.7,7.4,7.4,8.5,7.7,NA,NA,8.5,8.6,7.1,7.1))
    Mb <- as.numeric(subset(Interface.saturation.regions, SBZ == as.character(saturation.region), select = Mb))
    if(is.na(Mb)){
      Mb <- 7.9
    }
  }

# Constant ----------------------------------------------------------------
  
  if(region == "global"){
    c0 <- as.numeric(subset(coefficients.T, select = Global_c0))
  }else{
    c0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(saturation.region,"_c0",sep=""))))))
  }

# Path Term ---------------------------------------------------------------
  h <- 10^(-0.82 + 0.252*M)
  Rref <- sqrt(1 + h^2)
  R <- sqrt(Rrup^2 + h^2)
  LogR <- log(R)
  R_Rref <- log(R/Rref)
  
  #Need  to isolate regional anelastic coefficient, a0
  if(region == "global" | region == "Cascadia"){
    a0 <- as.numeric(subset(coefficients.T, select = c(Global_a0)))
  }else{
    a0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_a0",sep=""))))))
  }
  
  if(is.na(a0)){
    a0 <- as.numeric(subset(coefficients.T, select = c(Global_a0)))
  }

  Fp <- as.numeric(coefficients.T$c1*LogR + (coefficients.T$b4*M)*R_Rref + a0*R)
  

# Magnitude Scaling -------------------------------------------------------

  #compute M-scaling term
  func1 <- function(M,Mb){ifelse(M <= Mb,M-Mb,0)}
  func2 <- function(M,Mb){ifelse(M > Mb,M-Mb,0)}
  func3 <- function(M,Mb){ifelse(M <= Mb, (M-Mb)^2, 0)}
  Fm <- coefficients.T$c4*func1(M,Mb) + coefficients.T$c6*func2(M,Mb) + coefficients.T$c5*func3(M,Mb)

# Linear Site Amplification ----------------------------------------------

  #Site Coefficients
  V1 <- coefficients.T$V1..m.s.
  V2 <- coefficients.T$V2..m.s.
  Vref <- coefficients.T$Vref..m.s.
  
  if(region == "global"| region == "CAM"){
    s2 <- coefficients.T$Global_s2
    s1 <- s2
  }else if(region =="Taiwan" | region == "Japan"){
    s2 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_s2",sep=""))))))
    s1 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_s1",sep=""))))))
  }else{
    s2 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_s2",sep=""))))))
    s1 <- s2
  }
  
  #Compute linear site term
    if(VS30 <= V1){
      Flin <- s1*log(VS30/V1) + s2*log(V1/Vref)
    }else if(VS30 <= V2){
      Flin <- s2*log(VS30/Vref)
    }else{
      Flin <- s2*log(V2/Vref)
    }
 
# Nonlinear Site Term -----------------------------------------------------
  
  PGAr <- exp(GMM_at_760_IF_v5(event.type,region,saturation.region,Rrup,M,0,Rrup_b))
  f3 <- 0.05
  Vb <- 200
  Vref.Fnl <- 760
  
#  if(period >= 3){ #GAP 1/15/25: hard period cutoff at 3s was changed at late stage of publication, with tapering enforced in fnl coefficients instead for a smoother spectra
#    Fnl = 0;
 # }else{
    f2 <- coefficients.T$f4*(exp(coefficients.T$f5*(min(VS30,Vref.Fnl)-Vb))-exp(coefficients.T$f5*(Vref.Fnl-Vb)))
    Fnl = 0 + f2*log((PGAr+f3)/f3)
 # }

# Basin Term --------------------------------------------------------------
    #GAP 1/15/25: updated Cascadia to reflect version published in EQS
    #import coefficients
    
    if(Z2.5 == "default" | Z2.5 <= 0 | (region != "Japan" & region != "Cascadia")){
      Fb = 0
    }else{
      if(region == "Cascadia"){
        theta0 <- 3.75
        theta1 <- -0.74
        vmu <- 500
        vsig <- 0.42
        e1 <- coefficients.T$C_e1
        e2 <-  coefficients.T$C_e2 
        e3 <- coefficients.T$C_e3 
        
      }else if(region == "Japan"){
        theta0 <- 3.05
        theta1 <- -0.8
        vmu <- 500
        vsig <-0.33
        e3 <- coefficients.T$J_e3
        e2 <- coefficients.T$J_e2
        e1 <- coefficients.T$J_e1
      }
      
      Z2.5.pred <- 10^(theta0 + theta1*(1 + erf((log10(VS30) - log10(vmu))/(vsig*sqrt(2)))))
      delZ2.5 <- log(Z2.5) - log(Z2.5.pred)
      
      if(delZ2.5 <= (e1/e3)){
        Fb <- e1
      }else if(delZ2.5 >= (e2/e3)){
        Fb <- e2
      }else{
        Fb <- e3*delZ2.5
      }
    }
  
  # CAM updates concerning backarc anelastic attenuation and site response in 
  # the Valley of Mexico, VM (Dec/2024)
  if(region == "CAM"){
    Fp <- Fp + CAM_backarc_attenuation(period,Rrup_b,event.type)
    
    if(VM_flag == 1){
      site_response <- site_resp_VM(period,PGAr,VS30,Z_LCG)
      Fnl <- site_response[3]
      Flin <- site_response[4]
      Fb <- site_response[5]
    }
  }
  
  
# Add it all up! ----------------------------------------------------------
  
  mu <- c0 + Fp + Fnl + Fb + Flin + Fm
  return(c(mu))
}
