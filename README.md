# NGASub_CAM_regionalized
This implementation is an updated version of the Parker et al. (2022) Ground Motion Model, incorporating modifications specific to the Central America and Mexico (CAM) region. These updates address anelastic attenuation at backarc sites and site response in the Valley of Mexico (VM), including Mexico City. The users can compute PGA, PGV and 5%-damped elastic response spectra for individual scenarios. The backarc attenuation model is conditioned on the portion of the rupture distance in the backarc (Rrup,b) and earthquake type (interface or intraslab). The VM site response model is conditioned on peak ground acceleration at a reference rock site (PGAr), time-average shear wave velocity in the upper 30 m (VS30), and sediment depth measured to the Lower Coarse Grained layer (zLCG). The foundation of this revised model is detailed in the following two papers:

1. Contreras V, Stewart JP, Mayoral JM, Pérez-Campos X. (2025). Regionalization of Global Subduction Ground Motion Model for Application in Central America and Mexico, Earthquake Spectra, accepted.
2. Contreras V, Stewart JP, Mayoral JM, Pérez-Campos X. (2025). Valley of Mexico Site Response from Non-Reference Site Approach, Earthquake Spectra, https://doi.org/10.1177/87552930251316816.

NOTE: The model was not modified for any region outside of CAM.

The foundation of the original model is detailed in the following paper:
Parker G.A., Stewart J.P., Boore D.M., Atkinson G.M., Hassani B. (2022). NGA-subduction global ground motion models with regional adjustment factors. Earthquake Spectra, 38(1), 456-493.
