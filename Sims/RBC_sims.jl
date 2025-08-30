using MacroModelling
import StatsPlots


@model RBCsims begin

    1/c[0] = β * 1/c[1] * (r[1] + (1−δ))  

    k[1] = i[0] +(1−δ)* k[0] 

    y[0] = a[0] * k[0]^α *  n[0]^(1−α) 

    r[0] = α * a[0] * k[0]^(α−1) * n[0]^(1−α) 

    w[0] = (1−α)*a[0] * k[0]^α * n[0]^(−α) 
    
    θ * n[0]^χ  = 1/c[0] * w[0] 

    y[0] = c[0] +i[0] 

    log(a[0]) = ρ*log(a[-1])  + std_z * e[x]

end;

@parameters RBCsims begin
    
    std_z = 0.01
    χ=0.5
    θ = 0.5
    ρ = 0.95
    σ = 1
    δ = 0.02
    α = 0.33
    β = 0.95

end;

get_solution(RBCsims)

# get_steady_state(RBC)
# plot_irf(RBCsims)