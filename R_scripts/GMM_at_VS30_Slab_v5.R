#Grace Parker
#Modified February 26 to expand comments
#Modified March 25, 2020 to call consolidated coefficient table
#Modified December 21, 2024 to incorporate CAM updates concerning backarc
#anelastic attenuation and site response in the Valley of Mexico (VM) (Victor Contreras)

# Input Parameters --------------------------------------------------------

#Event type: 0 == interface, 1 == slab

#region corresponds to options in the DatabaseRegion column of the flatfile, plus global. Must be a string. If no matches, default will be global model:
# "global", "Alaska", "Cascadia", "CAM", "Japan", "SA" or "Taiwan"

#Saturation Region corresponds to regions defined by R. Archuleta and C. Ji:
# "global", "Aleutian","Alaska","Cascadia","Central_America_S", "Central_America_N", "Japan_Pac","Japan_Phi","South_America_N","South_America_S", "Taiwan_W","Taiwan_E"

# Rrup is number in kilometers

#Hypocentral depth in km. To use Ztor value to estimate hypocentral depth, see Ch. XXX of Parker et al. PEER report

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
# "GMM_at_VS30_Slab_v5.R" calls function "GMM_at_760_Slab_v5.R" to compute PGAr in the nonlinear site term. This function must be in the R environment else an error will occur.
# The output is the desired median model prediction in LN units
# Take the exponential to get PGA, PSA in g or  PGV in cm/s


#Function to compute GMM predictions at various VS30s for slab
GMM_at_VS30_Slab_v5 <- function(event.type,region,saturation.region,Rrup,M,hypocentral.depth,period,VS30, Z2.5,basin,Rrup_b,VM_flag,Z_LCG){
  
  if(event.type == 0){
    stop("This function is only for slab")
  }

# Import Master Coefficient Table -----------------------------------------
  
  coefficients <- read.csv("Table_S2_Slab_Coefficients.csv", header =  T)
  
  #Isolate desired period
  coefficients.T <- subset(coefficients,Period == period)
  
  #Define mb based on Archuleta and Ji (2019)
  if(saturation.region == "global"){
    Mb <- 7.6
  }else{
    slab.saturation.regions <- data.frame(SBZ = c("Aleutian","Alaska","-999","Cascadia","CAM_S", "CAM_N", "Japan_Pac","Japan_Phi", "New_Zealand_N","New_Zealand_S","SA_N","SA_S", "Taiwan_W","Taiwan_E"),
                                          Mb = c(7.98,7.2,NA,7.2,7.6,7.4,7.65,7.55,7.6,7.4,7.3,7.25,7.7,7.7))
    Mb <- as.numeric(subset(slab.saturation.regions, SBZ == as.character(saturation.region), select = Mb))
  }
  
# Constant ----------------------------------------------------------------
  
  #Isolate constant
  if(region == "global" | region == "Cascadia"){
    c0 <- as.numeric(subset(coefficients.T,select = Global_c0))
  }else if(region == "Alaska" | region == "SA"){
    c0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(saturation.region,"_c0",sep=""))))))
  }else{
    c0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_c0",sep=""))))))
  }
  
# Path Term ---------------------------------------------------------------  
 
  #near-source saturation
  if(M <= Mb){
    m <- (log10(35)-log10(3.12))/(Mb-4)
    h <- 10^(m*(M-Mb) + log10(35))
  }else{
    h <- 35
  }
  
  Rref <- sqrt(1 + h^2)
  R <- sqrt(Rrup^2 + h^2)
  LogR <- log(R)
  R_Rref <- log(R/Rref)
  
  #Need  to isolate regional anelastic coefficient, a0
  if(region == "global"){
    a0 <- coefficients.T$Global_a0
  }else{
    a0 <- as.numeric(subset(coefficients.T, select = c(eval(as.name(paste(region,"_a0",sep=""))))))
  }
  if(is.na(a0)){
    a0 <- coefficients.T$Global_a0
  }
  
  Fp <- as.numeric(coefficients.T$c1*LogR + (coefficients.T$b4*M)*R_Rref + a0*R)
  
# Magnitude Scaling ------------------------------------------------------- 

  func1 <- function(M,Mb){ifelse(M <= Mb,M-Mb,0)}
  func2 <- function(M,Mb){ifelse(M > Mb,M-Mb,0)}
  func3 <- function(M,Mb){ifelse(M <= Mb, (M-Mb)^2, 0)}
  
  Fm <- coefficients.T$c4*func1(M,Mb) + coefficients.T$c6*func2(M,Mb) + coefficients.T$c5*func3(M,Mb)
  

# Source Depth Scaling ----------------------------------------------------

  if(hypocentral.depth >= coefficients.T$db..km.){
    Fd <-coefficients.T$d
  }else if(hypocentral.depth <= 20){
    Fd <- coefficients.T$m*(20-coefficients.T$db..km.) + coefficients.T$d
  }else{
    Fd <- coefficients.T$m*(hypocentral.depth - coefficients.T$db..km.) + coefficients.T$d
  }
  
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
  
  PGAr <- exp(GMM_at_760_Slab_v5(event.type,region,saturation.region,Rrup,M,hypocentral.depth,0,Rrup_b))
  f3 <- 0.05
  Vb <- 200
  Vref.Fnl <- 760
  
  #if(period >= 3){  #GAP 1/15/25: hard period cutoff at 3s was changed at late stage of publication, with tapering enforced in fnl coefficients instead for a smoother spectra
   # Fnl = 0;
  #}else{
    f2 <- coefficients.T$f4*(exp(coefficients.T$f5*(min(VS30,Vref.Fnl)-Vb))-exp(coefficients.T$f5*(Vref.Fnl-Vb)))
    Fnl = 0 + f2*log((PGAr+f3)/f3)
  #}
  
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
  
   mu <- c0 + Fp + Fm + Fd + Fnl + Flin + Fb
  return(c(mu)) 
}


