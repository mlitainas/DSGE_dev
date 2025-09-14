# RBC model with governemnt Spendinf
# Notes from Nipsi-Landi 

using MacroModelling,DataFrames, StatsPlots
import StatsPlots

@model RBC_nl begin

    #  Households
    lambda[0] = c[0]^( -sigma)  

    lambda[0] = beta * lambda[1] / rr[0]

    lambda[0] = beta * lambda[1] * (rk[1] + (1 - delta) * q[1]) / q[0]  

    #  Labour demand
    h[0]^(phi) = w[0] * lambda[0]/kappaL 

    #  capital stock LoM 
    k[0] = (1-delta)*k[-1]+(1-kappaI/2*(i[0]/i[-1] - 1)^2) * i[0] 

    # Investment FOC
    -1 + q[0] * (1 - kappaI/2 *   (i[0]/i[-1] - 1)^2 - kappaI * (i[0]/i[-1] - 1 ) * i[0]/i[-1] ) + kappaI * beta * lambda[1] / lambda[0] *  q[1] * (i[1]/i[0] - 1 ) * (i[1]/i[0])^2      
    # -1 + q[0] *  (1 - kappaI/2 *   (i[0]/i[-1] - 1)^2 - kappaI * (i[0]/i[-1] - 1)  * i[0]/i[-1])  + kappaI * beta * lambda[1] / lambda[0] *  q[1] * (i[1]/i[0] - 1)  * (i[1]/i[0])^2

    #  Firms
    # Production function
    y[0] = a[0] *k[-1]^alpha * h[0]^(1-alpha)

    # labour demand
    w[0] = (1-alpha) * y[0]/h[0]

    # 11 capital demand - interest rate
    rk[0] = alpha * y[0] / k[-1] 

    # 10 Market clearing
    y[0] = c[0] + i[0] + g[0] 

    #  Shocks
    # 11 Technology
    log(a[0]) = (1-rhoa) * log(ass) + rhoa * log(a[-1]) + std_a * va[x]   
    # 12 Government spending
    log(g[0]) = (1-rhog) * log(gss) + rhog * log(g[-1]) + std_g * vg[x] + vg_news[-4]  
    # log(va_news[0]) =  rhog_news * log(va_news[-1]) + std_g_news * vg_news[x]  
    vg_news[0]  =  vg_news[x]  
end;

@parameters RBC_nl begin
    std_a = 0.01
    std_g = 0.01
    std_g_news = 0.01
    alpha=0.33;                   # elasticity of production wrt capital
    beta=0.99;                    # discount factor
 
    delta=0.025;                  # depreciation rate
    sigma=2;                      # relative risk aversion
    phi=1;                        # inverse of Frisch elasticity 
    
    ass = 0.1
    gss = 0.2
     # share of public spending in ss

    # Parameters not affecting the steady state
    kappaI=0;                 # investment adjustment cost. If 0, q is constant
    rhoa=0.9;                 # tfp persistence
    rhog=0.9                  #ß public spending persistence

    kappaL = 1
end;

# 
get_solution(RBC_nl)
# 
get_steady_state(RBC_nl)
# 
plot_irf(RBC_nl, variables = [:g, :y,:c,:w,:rk], shocks = :vg_news)
# plot_irf(RBC_nl;  variables = [:y, :c, :i], shock = "vg")
plot_irf(RBC_nl, variables = [:g, :y,:c,:w,:rk], shocks = :vg )

# 
sh_news = get_irfs(RBC_nl, variables = [:g, :y,:c,:w,:rk], shocks = :vg_news)
sh_g = get_irfs(RBC_nl, variables = [:g, :y,:c,:w,:rk], shocks = :vg)

# plot_irf(RBC_nl)

# Two subplots side by side
sh_g[:,:, end] 
