CAM_backarc_attenuation <- function(period,Rrup_b,event_flag){
# Function that estimates the backarc anelastic attenuation in Central
# America & Mexico (CAM) conditioned on T, Rrup_b, and event_flag.
# Ref: Contreras V, Stewart JP, Mayoral JM, Pérez-Campos X. (2025).
# Regionalization of Global Subduction Ground Motion Model for Application in
# Central America and Mexico (accepted)
#
# INPUT
# period (s)    : Oscillator period (PGV = -1, PGA = 0, T between 0.01s - 10s)
# Rrup_b (km)   : Backarc portion of the rupture distance
# event_flag    : 0 = interface, 1 = intraslab
#
# OUTPUT
# Fp_BA         : Backarc anelastic attenuation (in natural log units)
#
# Victor Contreras, Universidad Diego Portales
# 2024/12/14 (latest version)

# Backarc anelastic attenuation term
coeffs <- read.csv('CAM_backarc_attenuation_coeffs.csv', header =  TRUE, skip = 0)
IM <- coeffs$IM
if (event_flag == 0){
  a0b <- interp1(IM,coeffs$a0b_interface,period)
}else if(event_flag == 1){
  a0b <- interp1(IM,coeffs$a0b_intraslab,period)
}
Fp_BA <- a0b*Rrup_b
return(Fp_BA)
}