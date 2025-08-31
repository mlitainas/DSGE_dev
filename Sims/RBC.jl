using MacroModelling
import StatsPlots

@model RBCsims begin

    uc[0] = θ / c[0]

    # un[0] = -ψ / (1-n[0])

    n[0] = - ψ * c[0]/ θ * w[0] + 1 

    uc[0] = β *  uc[1] * (r[1] + (1 - δ))  

    k[1] = i[0] +(1 - δ)* k[0] 

    i[0] = y[0] - c[0]

    y[0] = a[0] * k[-1]^α *  n[0]^(1- α) 

    r[0] = α * y[0] / k[-1] 

    w[0] = (1 - α) * y[0] / n[0]
    
    # w[0] = -un[0]/uc[0]
    # n[0] = w[0] / c[0] -1  
    log(a[0]) = ρ*log(a[-1])  + std_z * e[x]

end;

@parameters RBCsims begin
    
    std_z = 0.01
    θ = 1
    ρ = 0.9
    σ = 1
    δ = 0.025
    α = 0.35
    β = 0.95
    ψ = 1.6

end;

# Policy function
get_solution(RBCsims)

plot_solution(RBCsims,:k)

# get_steady_state
get_steady_state(RBC)

# plot irf
plot_irf(RBCsims)

plot_irf(RBCsims, parameters = :α => 0.5, variables = [:y,:i,:r,:w,:c], shocks = :e)

get_irf(RBCsims)

plot_simulations(RBCsims, variables = [:y,:i,:r,:w,:c])

get_autocorrelation(RBC)

shock_series = KeyedArray(zeros(2,12), Shocks = [:e], Periods = 1:12)