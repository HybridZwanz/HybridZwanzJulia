include("../../LJChain.jl")
include("../../CGStats.jl")
# include("LJChain.jl")
# include("CGStats.jl")

using Main.LJChain,Main.CGStats,Statistics,DelimitedFiles

C = Chain();
S = Simulator(dt=1e-3, Stochdt=1e-3, FullSteps=Int(1e7), ApproxSteps=Int(1e7));

################################################################################

# Q_FGD = readdlm("/home/mspino/MarcoMasterProject/Julia_Code/Store_Trajectory/Outcomes/Q_FGD_1e7(B1.096)_CORRETTA.csv", Float64);
Q_FGD = readdlm("/home/mspino/MarcoMasterProject/Julia_Code/Store_Trajectory/Outcomes/Q_FGD_1e7.csv", Float64);
# Q_FGD = readdlm("Store_Trajectory/Outcomes/Q_FGD_1e5.csv", Float64);

lm = size(Q_FGD,1);
dQ1 = zeros(lm, C.NBeads);
for i in 1:lm
    dQ1[i, :] = circshift(Q_FGD[i, :], -1) - Q_FGD[i, :];
    dQ1[i, C.NBeads] += C.Period;
end

N = size(dQ1, 1);
l = N - 1;

msd1 = zeros(l);

n = 1;
sumtot = 0;
while n < N
    if n < 11
        for i in 1:(N-n)
            global sumtot += sum((dQ1[i+n,:] - dQ1[i,:]) .^ 2) / C.NBeads;
        end
        msd1[n] = sumtot / (N-n);
        global sumtot = 0;
    else
        n = Int(trunc(1.1 * n));
        if n < l
            for i in 1:(N-n)
                global sumtot += sum((dQ1[i+n,:] - dQ1[i,:]) .^ 2) / C.NBeads;
            end
            msd1[n] = sumtot / (N-n);
            global sumtot = 0;
        end
    end
    global n += 1;
end

v = findall(msd1 .!= 0);
msd1 = msd1[v];

# open("Data_Analysis/MSD/Outcomes/INT_msd_FGD_1e7(B1.096)_CORRETTA.csv", "w+") do x
#     writedlm(x, [ S.dt*100*(v.) msd1 ]);
# end

open("Data_Analysis/MSD/Outcomes/INT_msd_FGD_1e7.csv", "w+") do x
    writedlm(x, [ S.dt*100*v msd1 ]);
end

# open("Data_Analysis/MSD/Outcomes/INT_msd_FGD_1e5.csv", "w+") do x
#     writedlm(x, [ S.dt*100*v msd1 ]);
# end


################################################################################

# Q_DCGD = readdlm("/home/mspino/MarcoMasterProject/Julia_Code/Store_Trajectory/Outcomes/QApprox_DCGD_1e7(B1.096)_CORRETTA.csv", Float64);
Q_DCGD = readdlm("/home/mspino/MarcoMasterProject/Julia_Code/Store_Trajectory/Outcomes/QApprox_DCGD_1e7.csv", Float64);
# Q_DCGD = readdlm("Store_Trajectory/Outcomes/QApprox_DCGD_1e5.csv", Float64);

lm = size(Q_DCGD,1);
dQ2 = zeros(lm, C.NBeads);
for i in 1:lm
    dQ2[i, :] = circshift(Q_DCGD[i, :], -1) - Q_DCGD[i, :];
    dQ2[i, C.NBeads] += C.Period;
end

N = size(dQ2, 1);
l = N - 1;

msd2 = zeros(l);

n = 1;
sumtot = 0;
while n < N
    if n < 11
        for i in 1:(N-n)
            global sumtot += sum((dQ2[i+n,:] - dQ2[i,:]) .^ 2) / C.NBeads;
        end
        msd2[n] = sumtot / (N-n);
        global sumtot = 0;
    else
        n = Int(trunc(1.1 * n));
        if n < l
            for i in 1:(N-n)
                global sumtot += sum((dQ2[i+n,:] - dQ2[i,:]) .^ 2) / C.NBeads;
            end
            msd2[n] = sumtot / (N-n);
            global sumtot = 0;
        end
    end
    global n += 1;
end

v = findall(msd2 .!= 0);
msd2 = msd2[v];

open("Data_Analysis/MSD/Outcomes/INT_msd_DCGD_1e7.csv", "w+") do x
    writedlm(x, [ S.dt*100*v msd2 ]);
end


# open("Data_Analysis/MSD/Outcomes/INT_msd_DCGD_1e5.csv", "w+") do x
#     writedlm(x, [ S.dt*100*v msd2 ]);
# end


################################################################################

# Q_MMZD = readdlm("/home/mspino/MarcoMasterProject/Julia_Code/Store_Trajectory/Outcomes/QApprox_MMZD_1e7(B1.096)_CORRETTA.csv", Float64);
Q_MMZD = readdlm("/home/mspino/MarcoMasterProject/Julia_Code/Store_Trajectory/Outcomes/QApprox_MMZD_1e7(B1.096).csv", Float64); #NB: è calcolata con gli stessi parametri di quella sopra
# Q_MMZD = readdlm("Store_Trajectory/Outcomes/QApprox_MMZD_1e5(B1.096).csv", Float64);

lm = size(Q_MMZD,1);
dQ3 = zeros(lm, C.NBeads);
for i in 1:lm
    dQ3[i, :] = circshift(Q_MMZD[i, :], -1) - Q_MMZD[i, :];
    dQ3[i, C.NBeads] += C.Period;
end

N = size(dQ3, 1);
l = N - 1;

msd3 = zeros(l);

n = 1;
sumtot = 0;
while n < N
    if n < 11
        for i in 1:(N-n)
            global sumtot += sum((dQ3[i+n,:] - dQ3[i,:]) .^ 2) / C.NBeads;
        end
        msd3[n] = sumtot / (N-n);
        global sumtot = 0;
    else
        n = Int(trunc(1.1 * n));
        if n < l
            for i in 1:(N-n)
                global sumtot += sum((dQ3[i+n,:] - dQ3[i,:]) .^ 2) / C.NBeads;
            end
            msd3[n] = sumtot / (N-n);
            global sumtot = 0;
        end
    end
    global n += 1;
end

v = findall(msd3 .!= 0);
msd3 = msd3[v];

open("Data_Analysis/MSD/Outcomes/INT_msd_MMZD_1e7(B1.096).csv", "w+") do x
    writedlm(x, [ S.Stochdt*100*v msd3 ]);
end

# open("Data_Analysis/MSD/Outcomes/INT_msd_MMZD_1e5(B1.096).csv", "w+") do x
#     writedlm(x, [ S.Stochdt*100*v msd3 ]);
# end
