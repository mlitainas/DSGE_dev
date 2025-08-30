using MacroModelling
import StatsPlots

@model RBCsims begin

    uc[0] = θ / c[0]

    un[0] = -ψ / (1-n[0])

    uc[0] = β *  uc[1] * (r[1] + (1 - δ))  

    k[1] = i[0] +(1 - δ)* k[0] 

    i[0] = y[0] - c[0]

    y[0] = a[0] * k[-1]^α *  n[0]^(1- α) 

    r[0] = α * y[0] / k[-1] 

    w[0] = (1 - α) * y[0] / n[0]
    
    # w[0] = -un[0]/uc[0]
    # n[0] = w[0] / c[0] -1  

    n[0] = - ψ * c[0]/ θ * w[0] + 1 

    log(a[0]) = ρ*log(a[-1])  + std_z * e[x]

end;

@parameters RBCsims begin
    
    std_z = 0.01
    θ = 0.5
    ρ = 0.95
    σ = 1
    δ = 0.02
    α = 0.33
    β = 0.95
    ψ = 1

end;

get_solution(RBCsims)

# get_steady_state(RBC)
 plot_irf(RBCsims)