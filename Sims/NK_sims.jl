using MacroModelling
import StatsPlots
using Plots

@model NKsims begin

    uc[0] =  c[0]^(-σ)

    uc[0] = β *  uc[1] * (1 + i[0]) / pi[1] 

    # 1
    lamda[0] = β *  (c[0] /c[-1]) ^(-σ)

    # 2
    θ * n[0]^χ = w[0] * uc[0]

    # 3
    m[0] = ψ * c[0]^σ * ((1+i[0])/i[0])

    (1-ϕ)* pisharp[0]^(1-ϵ) + ϕ * pi[0]^(ϵ -1) - 1

    pisharp[0] = ϵ/(ϵ-1) * x1[0] / x2[0]

    x1[0] = mc[0]* y[0] + ϕ * lamda[1] * pi[1]^ϵ * x1[1]

    x2[0] = y[0] + ϕ * lamda[1] * pi[1]^(ϵ-1) * x2[1]

    mc[0] = w[0] / a[0]

    y[0] =  c[0]

    a[0] * n[0] = y[0] * up[0]  

    up[0] = (1-ϕ) * pisharp[0]^(-ϵ) + ϕ* pi[0]^ϵ * up[-1]
       
    log(a[0]) = rhoa * log(a[-1])  + sda * ea[x]
    
    gm[0] = (1 - rhom)*gm[ss] * rhom * gm[-1] + sdm * em[x]

    gm[0] = log(m[0]/m[-1]) + log(pi[0]) 
end;



@parameters NKsims begin  
 ϕ = 2/3
 σ = 3 
 β = 0.99
 χ = 1
 θ = 1
 ϵ = 11
 ψ = 1
 gm = 0
 rhom = 0.5
 rhoa = 0.9
 sdm = 0.01
 sda = 0.01
end;


get_solution(NKsims)
# get_steady_state(NKsims)
plot_irfs(NKsims)

# plot_solution(NKsims, :pi, variables = [:a,:])


irf = get_irfs(NKsims)
# # plot_simulations(RBCsims)
# plot IRF
horizon = 20

# 1 for productivity, 2 for monetary 
Shock = 2
irf_dt = irf[:,1:horizon,Shock]*100;

# variable names in the same order as slice rows
#vars = [:a, :c, :lamda, :mc, :pi, :w, :y, :m, :pisharp, :up]
vars = [:m, :y, :pi, :i, :mc, :gm];

row_index = findfirst(==( :gm ), vars);   # returns 2

irf_dt[row_index,:] =  cumsum(irf_dt[row_index,:] )

# build one small plot per variable
plots = [plot(1:horizon, irf_dt[i, :], title=string(v), legend=false) for (i, v) in enumerate(vars)];

# arrange them in a grid layout (e.g. 4 rows × 4 columns for 16 vars)
#plot(plots..., layout=(5,3))
plot(plots..., layout=(4,3), size=(800,800))

##############################################
#        TFP shock
##############################################
# 1 for productivity, 2 for monetary 
Shock = 2
irf_dt = irf[:,1:horizon,Shock]*100;

# variable names in the same order as slice rows
#vars = [:a, :c, :lamda, :mc, :pi, :w, :y, :m, :pisharp, :up]
vars = [:m, :y, :pi, :i, :mc, :gm];

row_index = findfirst(==( :gm ), vars);   # returns 2

irf_dt[row_index,:] =  cumsum(irf_dt[row_index,:] )

# build one small plot per variable
plots = [plot(1:horizon, irf_dt[i, :], title=string(v), legend=false) for (i, v) in enumerate(vars)];

# arrange them in a grid layout (e.g. 4 rows × 4 columns for 16 vars)
#plot(plots..., layout=(5,3))
plot(plots..., layout=(4,3), size=(800,800))
