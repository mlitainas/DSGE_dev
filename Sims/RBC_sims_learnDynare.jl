using MacroModelling
import StatsPlots

@model SimsDynare begin

    uc[0] =  c[0]^(-σ)

    uc[0] = β *  uc[1] * (r[1] + (1 - δ))  

    k[0] = i[0] +(1 - δ)* k[-1] 

    i[0] = y[0] - c[0]

    y[0] = a[0] * k[-1]^α 

    uc[0] = β *  uc[1] * (1 + rb[0])  

    r[0] = α * y[0] / k[-1] 

    w[0] = (1 - α) * y[0] 
    
    log(a[0]) = ρ*log(a[-1])  + std_z * e[x]

end;

@parameters SimsDynare begin
    
    std_z = 0.01
    ρ = 0.95
    σ = 1
    δ = 0.025
    α = 1/3
    β = 0.99


end;


get_solution(SimsDynare)

get_steady_state(SimsDynare)

plot_irf(SimsDynare)

# plot_simulations(RBCsims)
