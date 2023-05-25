Nk = 400 #number of k point mesh
NT = 20 #number of section for temperature
boltz = 1.380649 * 10^(-23) 
ele_charge = 1.602176634*10^(-19)
fac = 10^6 * boltz / ele_charge
tau =  1

function calc_mu(filling)
    
end 

for imu = 0:5
    filling = 1.7 + 0..05 * imu 

    for iT = 1:NT 
        T = iT/NT * 0.3
        mu = calc_mu(filling)

    end
end  

