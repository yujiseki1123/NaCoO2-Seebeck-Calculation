#Constants 
Nk = 400 #Number of k point mesh
NT = 20 #Number of section for temperature
boltz = 1.380649 * 10^(-23) #Boltzmann constant
electron_charge = 1.602176634*10^(-19)　#Electron charge
fac = 10^6 * boltz / electron_charge #Factor constant
tau =  1 #Relaxation time 
Nhop = 18 #Number of hopping coefficients
tpi = 2 * pi 

#Array of position and its hopping coefficients
pos = Array[]
HR = []

#Reading position and hopping parameter from posHR.txt file 
open("posHR.txt", "r") do f
    for line in eachline(f)
        line_num = split(line)
        pos_ele = []
        num_1 = parse(Float64,line_num[1])
        num_2 = parse(Float64,line_num[2])
        num_3 = parse(Float64,line_num[3])
        push!(pos_ele,num_1)
        push!(pos_ele,num_2)
        push!(pos,pos_ele)
        push!(HR,num_3)
    end
end 

#Reciprocal lattice vector 
b1 = (0.5, 0.5 * sqrt(3)) 
b2 = (-0.5,0.5 * sqrt(3))

#Function to calculate energy of wavelength k = (kx, ky)
function calc_Ek(kx,ky)
    Ek = 0 
    for i = 1 : Nhop
        theta = tpi * (kx * pos[i][1] + ky*pos[i][2])
        Ek = Ek + HR[i] * cos(theta)
    end 
    return Ek 
end 

#Function to calculate the dispersion 
function calc_fk(Ek,mu,T)
    r1 = (Ek - mu) / T 
    fk = 0 
    if r1 > 600
        fk = 0 
    elseif r1 < -600
        fk = 1 
    else 
        fk = 1 / (exp(r1) + 1)
    end 
    return fk 
end 

#Function to calulate the velocity for k_x direction
function calc_vk_x(kx,ky)
    vk_x = 0 
    for i = 1 : Nhop
        theta = tpi * (kx * pos[i][1] + ky * pos[i][2])
        vk_x = vk_x - HR[i] * sin(r1) * pos[i][1]
    end 
    return vk_x 
end 

#Function to calulate the velocity for k_y direction 
function calc_vk_xy(kx,ky)
    vk_y = 0 
    for i = 1 : Nhop
        theta = tpi * (kx * pos[i][1] + ky * pos[i][2])
        vk_y = vk_y - HR[i] * sin(r1) * pos[i][2]
    end 
    return vk_x 
end 

#Function to calculate the derivitive of fk
function calc_fk_deriv(Ek,mu,T,fk) 
    r1 = (Ek - mu)/T 
    fk_deriv = 0 
    if -600 <= r1 <= 600
        fk_deriv = - exp(r1) / T / (exp(r1) + 1)^2 
    end 
    return fk_deriv 
end 

#Function to calculate mu by binary search 
function calc_mu(filling,T)
    mu_max = 7 
    mu_min = 0 
    mu_mid = (mu_max + mu_min) / 2 
    mu = 0 

    for i2 = 1 : 1000
        i1 = 0 
        fmin = 0 
        fmid = 0 
        fmax = 0 
        for k2 = -Nk : (Nk-1) 
            for k1 = -Nk : (Nk-1)
                i1 = i1 + 1 
                r1 = (k1 / Nk) * 0.5 
                r2 = (k2 / Nk) * 0.5 
                kx = r1 * b1[1] + r2 * b2[1]
                ky = r1 * b1[2] + r2 * b2[2]
                Ek = calc_Ek(kx,ky)
                fmin = fmin + 2 * calc_fk(Ek,mu_min,T) 
                fmid = fmid + 2 * calc_fk(Ek,mu_mid,T) 
                fmax = fmax + 2 * calc_fk(Ek,mu_max,T)
            end 
        end
        fmin = fmin / i1 
        fmid = fmid / i1 
        fmax = fmax / i1

        if fmid < (filling + 10^(-8)) 
            mu_min = mu_mid 
            mu_mid = (mu_min + mu_max)/2 
        elseif fmid > (filling + 10^(-8)) 
            mu_max = mu_mid 
            mu_mid = (mu_min + mu_max)/2
        else 
            mu = mu_mid 
            break 
        end 
    end 

    return mu 
end 

#Start Calulation of Seebeck Coefficients
for imu = 0:5
    filling = 1.7 + 0.05 * imu 

    for iT = 1:NT 
        T = iT/NT * 0.3
        mu = calc_mu(filling,T)
        i1 = 0 
        dos = 0 
        fmid = 0 
        K0mat = zeros(2,2)
        K1mat = zeros(2,2)
        for k2 = -Nk : Nk-1
            for k1 = -Nk : Nk-1 
                i1 = i1 + 1 
                r1 = (k1/Nk) * 0.5
                r2 = (k2/Nk) * 0.5
                kx = r1 * b1[1] + r2 * b2[1]
                ky = r1 * b1[2] + r2 * b2[2]
                Ek = calc_Ek(kx,ky)
                print(Ek)
                vk_x = calc_vk_x(kx,ky)
                vk_y = calc_vk_y(kx,ky)
                fk = calc_fk(Ek,mu,T)
                fk_deriv = calc_fk_deriv(Ek,mu,T,fk)   
                dos = dos - fk_deriv * 2 

                #Update K0mat 
                K0mat[1,1] += 2 * tau * vk_x * vk_x * (-fk_deriv)
                K0mat[2,1] += 2 * tau * vk_y * vk_x * (-fk_deriv)
                K0mat[1,2] += 2 * tau * vk_x * vk_y * (-fk_deriv)
                K0mat[2,2] += 2 * tau * vk_y * vk_y * (-fk_deriv)

                #Update K1mat 
                K1mat[1,1] += 2 * tau * vk_x * vk_x * (-fk_deriv)*(Ek-mu)
                K1mat[2,1] += 2 * tau * vk_y * vk_x * (-fk_deriv)*(Ek-mu)
                K1mat[1,2] += 2 * tau * vk_x * vk_y * (-fk_deriv)*(Ek-mu)
                K1mat[2,2] += 2 * tau * vk_y * vk_y * (-fk_deriv)*(Ek-mu)
            end 
        end 

        fmid = fmid / i1 
        dos = dos / i1 
        K0mat ./ i1 
        K1mat ./ i1 

        Smat = inv(K0mat) * K1mat
        Smat ./ -(T*fac)
    end
end  

