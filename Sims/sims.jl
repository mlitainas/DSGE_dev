using MacroModelling
import StatsPlots

@model RBC begin

    1/c[0]  =  β * 1/c[1] * ( α * a[1] * k[0]^(α-1)+(1-δ))

    k[0] = a[0]* k[-1]^(α) - c[0] + (1-δ) * k[-1]
    
    a[0] = ρ * a[-1] + std_z * eps_z[x]

    # a[0] = log(z[0]) 
end;


# @model RBC begin
#     c[0]^(-1) = β  * 1/c[1] * (α * exp(z[1]) * k[-1]^(α - 1) + (1 - δ))
#     k[0] = (1 - δ) * k[-1] + q[0] - c[0]
#     q[0] = exp(z[0]) * k[-1]^α
#     z[0] = ρ * z[-1] + std_z * eps_z[x]
# end;


@parameters RBC begin
    std_z = 0.01
    ρ = 0.9
    δ = 0.02
    α = 0.33
    β = 0.95
end;

# get_solution(RBC)

# get_steady_state(RBC)
plot_irf(RBC)