# Sensing Aided OTFS Channel Estimation for Massive MIMO Systems

This is a Python code package related to the following article:
S. Jiang and A. Alkhateeb, "[Sensing Aided OTFS Channel Estimation for Massive MIMO Systems](https://arxiv.org/abs/2209.11321)," in IEEE ICC Workshops, 2023

# Abstract of the Article
Orthogonal time frequency space (OTFS) modulation has the potential to enable robust communications in highlymobile scenarios. Estimating the channels for OTFS systems, however, is associated with high pilot signaling overhead that scales with the maximum delay and Doppler spreads. This becomes particularly challenging for massive MIMO systems where the overhead also scales with the number of antennas. An important observation however is that the delay, Doppler, and angle of departure/arrival information are directly related to the distance, velocity, and direction information of the mobile user and the various scatterers in the environment. With this motivation, we propose to leverage radar sensing to obtain this information about the mobile users and scatterers in the environment and leverage it to aid the OTFS channel estimation in massive MIMO systems. As one approach to realize this vision, this paper formulates the OTFS channel estimation problem in massive MIMO systems as a sparse recovery problem and utilizes the radar sensing information to determine the support (locations of the non-zero delay-Doppler taps). The proposed radar sensing aided sparse recovery algorithm is evaluated based on an accurate 3D raytracing framework with co-existing radar and communication data. The results show that the developed sensing-aided solution consistently outperforms the standard sparse recovery algorithms (that do not leverage radar sensing data) and leads to a significant reduction in the pilot overhead, which highlights a promising direction for OTFS based massive MIMO systems.# License and Referencing

# License and Referencing
This code package is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/). 
If you in any way use this code for research that results in publications, please cite our original article:
> S. Jiang and A. Alkhateeb, "Sensing Aided OTFS Channel Estimation for Massive MIMO Systems," in IEEE ICC Workshops, 2023
