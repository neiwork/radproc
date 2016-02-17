#pragma once

#include <fparticle\particle.h>



/* pairBH estimates the pair injection due to the channel 
for pair production of photohadronic interaction; Bethe Heitler*/ 
double pairBH(double E, Particle& particle, Particle& proton, fun1 tpf);  
