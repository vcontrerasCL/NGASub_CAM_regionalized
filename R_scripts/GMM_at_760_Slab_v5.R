#Grace Parker
#Modified February 26 to expand comments
#Modified March 25, 2020 to call coefficients from master table
#Modified December 16, 2024 to incorporate CAM updates concerning backarc
#anelastic attenuation (Victor Contreras)

# Input Parameters --------------------------------------------------------

#Event type: 0 == interface, 1 == slab

#region corresponds to options in the DatabaseRegion column of the flatfile, plus global. Must be a string. If no matches, default will be global model:
# "global", "Alaska", "Cascadia", "CAM", "Japan", "SA" or "Taiwan"

#Saturation Region corresponds to regions defined by R. Archuleta and C. Ji:
# "global", "Aleutian","Alaska","Cascadia","Central_America_S", "Central_America_N", "Japan_Pac","Japan_Phi","South_America_N","South_America_S", "Taiwan_W","Taiwan_E"

# Rrup is number in kilometers

#Hypocentral depth in km. To use Ztor value to estimate hypocentral depth, see Ch. 4.3.3 of Parker et al. PEER report

# period can be: (-1,0,0.01,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1,1.5,2,2.5,3,4,5,7.5,10) 
# where -1 == PGV and 0 == PGA

# Rrup_b is the backarc portion of the rupture distance in km (CAM only, Dec/24 updates)
# File "CAM_backarc_attenuation_coeffs.csv" with coefficients must be in the active directory

# Other pertinent information ---------------------------------------------
#Coefficient files must be in the active working directory
# This function has no site term. Can only estimate ground motion at the reference condition VS30 = 760m/s
# The output is the desired median model prediction in LN units
# Take the exponential to get PGA, PSA in g or the PGV in cm/s

#Function to compute GMM predictions at 760m/s for slab and interface
GMM_at_760_Slab_v5 <- function(event.type,region,saturation.region,Rrup,M,hypocentral.depth,period,Rrup_b){
  
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
  
  # CAM updates concerning backarc anelastic attenuation (Dec/2024)
  if(region == "CAM"){
    Fp <- Fp + CAM_backarc_attenuation(period,Rrup_b,event.type)
  }
  
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
  

   mu <- c0 + Fp + Fm + Fd
  return(mu) 
}
