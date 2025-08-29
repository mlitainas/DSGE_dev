using MacroModelling
import StatsPlots

@model NKSOE begin
        
#  Household block
#  Households C foc 
lambda[0] = ( c[0] - kappaL *  h[0]^(1+phi) / (1+phi) )^(-sigma)

# Households euler equation FoC wrt Domestic bonds br
lamda[0] = beta * lambda[1] * r[0] /pii[1] 

# Households euler equation FoC wrt Forein bonds  bf
-kappaD * (d[0]-d[ss]) + 1 - ( beta * lambda[1] / lambda[0] * rz[0] * s[1]/ s[0])  

# Households euler equation FoC wrt capital 
q[0] = beta * ( lambda[1]/lambda[0] * ( rk[1] + (1-delta) * q[1] ) ) 

# Households euler equation FoC wrt labour 
w[0] = beta *  kappaL *  h[0]^(phi)       

# # Households FoC wrt investment 
- 1/q[0] + (1-kappaI/2 * ((i[0]/i[-1])-1)^2 - kappaI * (i[0]/i[-1]) * (i[0]/i[-1]-1) *i[0]/i[-1])+kappaI*beta*lambda[1]/lambda[0]*q[1]*(i[1]/i[0]-1)*(i[1]/i[0])^2    

# Low of motion for capital
k[0] = (1-delta)*k[-1] + (1-kappaI/2 * (i[0]/i[-1]-1)^2) * i[0]    

# Firms
# Production function
yH[0]= a[0]* k[-1]^alphaa * h[0]^(1-alphaa)

# labour demand
w[0] = (1-alphaa) * mc[0] * yH[0]/h[0] 

# capital demand
rk[0] = alphaa * mc[0] * yH[0] / k[-1]

# Philips curve
(piH[0]-pii[ss])*piH[0]=beta*(lambda[1]/lambda[0]*pH[1]*yH[1]/(pH[0]*yH[0])*piH[1]*(piH[1]-pi[ss]))+epsilon/kappaP*(mc[0]/pH[0]-(epsilon-1)/epsilon)

# Domestic goods prices
piH[0] = pH[0]/pH[-1] * pii[0]

# Market clearing
yH[0]=(1-gamma)*pH[0]^(-eta)*(c[0]+i[0])+g[0]+gammaz[0]*yz[0]*(tot[0])^eta+(kappaP/2*(piH[0]-pi[ss])^2)*yH[0] 

gdp[0]=c[0]+i[0]+pH[0]*g[0] - s[0]*d[0]+s[0]*rz[-1]*d[-1]+(kappaP/2*(piH[0]-pii[ss])^2)*gdp[0]+s[0]*kappaD/2*(d[0]-d[ss])^2


# Prices
-pH[0]^(1-eta) + 1 - gamma*s[0]^(1-eta) /(1-gamma)

# Policy
-r[0]/r[ss] + ( (pii[0] / pii[ss])^(phipi) * (gdp[0])^(phiy) * (De[0]/pii[ss])^phie)^(1-rhom)*(r[-1]/r[ss])^(rhom) * exp(vm[x]);

# Shocks
# Technology
log(a[0])  = (1-rhoa)*log(a[ss]) + rhoa * log(a[-1]) + va[x]
# Government spenidng 
log(g[0])  = (1-rhog)*log(g[ss]) + rhog * log(g[-1]) + vg[x]  
# Forein output
log(yz[0]) = rhoy * log(yz[-1]) + vy[x] 

# foreign monetary shock
rz[0]= (1-rhop) * 1/beta + rhop * rz[-1] + vp[x]

# Auxiliary variables
# rr[0]=r[0]/pii[1]
# s[0]/s[-1]=De[0]/pii[0]
# gdp[0]=pH[0]*yH[0]
# tb[0]=xp[0]-mp[0]
# xp[0]=pH[0]*gammaz[0]*yz[0]*(tot[0])^eta
# mp[0]=gamma*s[0]^(1-eta)*(c[0]+i[0])
# tot[0]=s[0]/pH[0]
# clog[0]=log(c[0]) 
# wlog[0]=log(w[0]) 
# hlog[0]=log(h[0]) 
# klog[0]=log(k[0]) 
# ilog[0]=log(i[0])
# dlog[0]=log(s[0]*d[0]/(4*gdp[0]))
# proof[0]=gdp[0]-(c[0]+i[0]+pH[0]*g[0]+tb[0]+(kappaP/2*(piH[0]-pii[ss])^2)*gdp[0])

end;



# model_summary(NKSOE)




# @parameters NKSOE begin
# # Structural Parameters (quarterly calibration)
# beta = 0.99                    # discount factor
# alphaa = 0.33                   # elasticity of production wrt capital
# epsilon = 6                    # elasticity of substitution btw differentiated goods
# delta = 0.025                  # depreciation rate
# sigma = 2                      # relative risk aversion
# phi = 1                        # inverse of Frisch elasticity
# g = 0.2                        # share of public spending in ss
# eta = 1.5                      # elasticity of intratemporal substitution
# gamma = 0.3                    # trade openness
# D = 0.25                       # external debt/GDP ratio


# ## Parameters not affecting the steady state
# phipi=1.5                # mp response to inflation
# phiy=0.125               # mp response to output
# phie=0                   # mp response to exchange rate
# kappaI=2.48              # investment adjustment cost (as in CEE). If 0, q is constant
# kappaD=0.001             # usually calibrated at a small value. If 0, the model has a unit root
# rhoa=0.9                 # tfp persistence
# rhog=0.9                 # public spending persistence
# rhop=0.9                 # foreign mp shock persistence
# rhoy=0.9                 # foreign demand shock persistence persistence
# rhom=0.8                 # monetary policy inertia
# calvo=0.66               # price rigidity in calvo framework
# #  adjusment cost coefficient to have the same linear Phillips Curve of the Calvo framework complicated
# kappaP=(epsilon-1)*calvo/(piss^2*(1-calvo)*(1-beta*calvo)) 


# end;


# # plot_irf(NKSOE)


# #  Steady State
# # pii = 1;                        % inflation targeting (quarterly calibration)
# # d = 4*D;                       % foreign debt
# # yH=1;                        % domestic output
# # pH=1;                        % domestic price
# # s=1;                         % real FX rate
# # h=1/3;                       % hours of work
# # rr=1/beta;                   % real interest rate
# # r=pii/beta;                   % nominal interest rate
# # q=1;                         % marginal value of investment (in terms of lambda)  
# # rk=1/beta-(1-delta);         % rental rate of capital
# # mc=(epsilon-1)/epsilon;      % real marginal costs
# # k=alpha*mc*yH/rk;            % capital
# # w=(1-alpha)*mc*yH/h;         % real wage
# # kappaL=w/(h^(phi));          % labor preference parameter
# # i=delta*k;                   % investment
# # c=yH-i-g-d*(1/beta-1);       % consumption
# # lambda=(c-kappaL/(1+phi)*h^(1+phi))^(-sigma);   % marginal utlity of consumption
# # a=yH/(k^(alpha)*h^(1-alpha)); % tfp
# # gammaz=1-(1-gamma)*(c+i)-g;   % foreign parameter


# # pi=1;  
# # piH=1;
# # d=dss;                       
# # yH=1;                        
# # pH=1;                        
# # s=1;                         
# # h=1/3;                       
# # rr=1/beta;                   
# # r=pi/beta;                   
# # q=1;                         
# # rk=1/beta-(1-delta);         
# # mc=(epsilon-1)/epsilon;      
# # k=alpha*mc*yH/rk;            
# # w=(1-alpha)*mc*yH/h;        
# # i=delta*k; 
# # g=gss;
# # c=yH-i-g-d*(1/beta-1);       
# # lambda=(c-kappaL/(1+phi)*h^(1+phi))^(-sigma);   
# # a=ass; 
# # rz=1/beta;
# # yz=1;
# # De=1;
# # gdp=pH*yH;
# # tot=s/pH;
# # xp=pH*gammaz*yz*(tot)^eta;
# # mp=gamma*s^(1-eta)*(c+i);
# # tb=xp-mp+s*kappaD/2*(d-dss)^2;
# # clog=log(c);
# # wlog=log(w);
# # hlog=log(h); klog=log(k); 
# # ilog=log(i); 
# # dlog=log(s*d/(4*gdp));
# # proof=0;