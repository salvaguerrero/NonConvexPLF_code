using PlotlyJS
using JuMP
using Gurobi
using Distributions
using PiecewiseLinearOpt


############################################  PiecewiseLinearOpt: HelloWorld Univariate ############################################

model = Model()
set_optimizer(model, Gurobi.Optimizer)
@variable(model, x)
d = 0:(pi/4):4pi
y_range = [sin(x) for x in d]
z = piecewiselinear(model, x, d, sin, method=:Incremental)
# = piecewiselinear(model, x, UnivariatePWLFunction(d,sin))
@objective(model, Max, z)

#Solve 
set_optimizer_attribute(model, "MIPGap", 0.1)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))

#Ploting
if termination_status(model) == OPTIMAL
    plot(d,y_range)
    scatter!([JuMP.value(x)],[JuMP.value(z)],label="min")    
end

############################################  PiecewiseLinearOpt: HelloWorld Bivariate ############################################

using Plots
using JuMP
using Gurobi
using Distributions
using PiecewiseLinearOpt
x_range = -4:1:4
y_range = -4:1:4
z_range = [ exp(x+y) for x in x_range, y in y_range]
model = Model()
set_optimizer(model, Gurobi.Optimizer)
@variable(model, x)
@variable(model, y)
@variable(model, z)

z = piecewiselinear(model, x, y, x_range, y_range, (u,v) -> exp(u+v))
@objective(model, Min, z)

#Solve 
set_optimizer_attribute(model, "MIPGap", 0.1)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))

println("[" ,JuMP.value(x),";",JuMP.value(y),";",JuMP.value(z),"]")

#Ploting
if termination_status(model) == OPTIMAL

    surface(x_range,y_range,z_range)
    scatter!([JuMP.value(x)],[JuMP.value(y)],[JuMP.value(z)],label="min")

end


############################################  PiecewiseLinearOpt: Highly non convex ############################################
f(x,y) = 3*(1-x)^2*exp(-(x^2) - (y+1)^2)  - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2)  - 1/3*exp(-(x+1)^2 - y^2)
g(x,y) = x + y + 3

x_range = -4:1:4
y_range = -4:1:4

sol = [[],[],[]]
model = Model()
set_optimizer(model, Gurobi.Optimizer)
@variable(model, x)
@variable(model, y)
@variable(model, z)

z = piecewiselinear(model, x, y, x_range, y_range, (u,v) -> f(u,v))
@constraint(model,z <= x + y + 5)

@objective(model, Max, z)

#Solve 
set_optimizer_attribute(model, "MIPGap", 0.1)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))

println("[" ,JuMP.value(x),";",JuMP.value(y),";",JuMP.value(z),"]")



#Ploting
if termination_status(model) == OPTIMAL
    push!(sol[1],JuMP.value(x))
    push!(sol[2],JuMP.value(y))
    push!(sol[3],JuMP.value(z))
end
@objective(model, Min, z)

#Solve 
set_optimizer_attribute(model, "MIPGap", 0.1)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))

#Ploting
if termination_status(model) == OPTIMAL
    push!(sol[1],JuMP.value(x))
    push!(sol[2],JuMP.value(y))
    push!(sol[3],JuMP.value(z))
end

surface(x_range, y_range,f)
surface!(x_range, y_range,g,alpha = 0.5)
scatter!(sol[1],sol[2],sol[3],label="min & max")




############################################  PiecewiseLinearOpt: Hydro Plant ############################################

using PlotlyJS
using JuMP
using Gurobi
using Distributions
using PiecewiseLinearOpt

#Data------------------------------------------>
#Generators
Opex  = [0,25]        #[€/MWh] Hydro, Gas
L_lim = [0,0]         #[MW] 
U_lim = [10000,10000] #[MW]

#Demand
Time = 1:1
T_num = length(Time)
dem = [200*sin(t*(2*pi)/24)+300 for t in Time] 

#Hydro Plant
#Water levels [#hm3]
W_lev_0 = 10
W_lev_f = 10
w_apo   = 20#aportaciones [m3/h]
w_lev_M = 100 
w_lev_m = 0 
#Turbined water [#hm3]
w_tur_M = 100  
w_tur_m = 0    
# cubic reservoir a x b [m]
a = 10        
b = 10
# efficiency
eff = 0.6    

model = Model()
set_optimizer(model, Gurobi.Optimizer)
@variable(model, g_gen   >= 0)
@variable(model, g_gen_h >= 0)
@variable(model, w_tur   >= 0)
@variable(model, w_lev_0 >= 0)
@variable(model, w_lev_f >= 0)
@variable(model, w_spi   >= 0)

@objective(model, Min, g_gen*Opex[2] )

@constraint(model,Demand, g_gen_h + g_gen == dem[1] ) 

#@constraint(model, PowerLimit[g in 1:2, t in Time], L_lim[g] <= g_gen[g,t] <= U_lim[g] )

@constraint(model,  w_lev_f == w_lev_0 - w_tur + w_apo - w_spi )
@constraint(model, w_lev_0 == W_lev_0 )
#@constraint(model, w_lev[T_num] == w_lev_f )

p(q,v) = 9.81*q*v#/a/b*eff

q_range = w_tur_m:1:w_tur_M
v_range = w_lev_m:1:w_lev_M
#t = 1
#g_gen[1,t] == piecewiselinear(model, w_tur[t], w_lev[t-1], q_range, v_range, (q,v) -> p(q,v))#,method=:ZigZag)
g_gen_h = piecewiselinear(model, w_tur, w_lev_f, q_range, v_range, p)#,method=:ZigZag)


#@constraint(model, WatTurb, g_gen[1,t] == 9.81*w_tur[t]*eff*0 )
#@constraint(model, WaterLimit[ t in Time], w_tur_m <= w_tur + w_spi <= w_tur_M )


set_optimizer_attribute(model, "MIPGap", 0.1)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))

JuMP.value(g_gen_h) 
JuMP.value(g_gen)
JuMP.value(w_tur)
JuMP.value(w_lev_0)

Gen = Dict()
Gen[1] = zeros(T_num)
Gen[2] = zeros(T_num)
Level  = zeros(T_num)
Turbined = zeros(T_num)
STMP = zeros(T_num)
Spilled = zeros(T_num)
Inflows = w_apo*ones(T_num)
if termination_status(model) == OPTIMAL
    for t in Time
        Level[t]    = round(JuMP.value(w_lev[t])   , digits = 2)
        Turbined[t] = round(JuMP.value(w_tur[t])   , digits = 2)
        Spilled[t]  = round(JuMP.value(w_spi[t])   , digits = 2)
        Gen[1][t]   = round(JuMP.value(g_gen_h) , digits = 2)
        Gen[2][t]   = round(JuMP.value(g_gen[2,t]) , digits = 2)
        STMP[t]     = 0#round(dual(Demand[t])        , digits = 2)
    end

    g1_power = scatter(;x=Time, y=Gen[1], mode="lines+markers",name="Hydro",stackgroup="one")
    g2_power = scatter(;x=Time, y=Gen[2], mode="lines+markers",name="Gas"  ,stackgroup="one")
    demand   = scatter(;x=Time, y=dem,    mode="markers"      ,name="Demand")
  
    water_lev   = scatter(;x=Time, y=Level,   mode="lines+markers",name="Water Level")

    water_inf   = scatter(;x=Time, y=Inflows,mode="lines+markers" ,name="Water Inflows")
    water_tur   = scatter(;x=Time, y=Turbined,mode="lines+markers",name="Turbined water")
    water_spill = scatter(;x=Time, y=Spilled,mode="lines+markers" ,name="Spilled water")

    price       = scatter(;x=Time, y=STMP,mode="lines+markers",name="STMP")

    p3 = [plot([g1_power,g2_power,demand],       Layout(title="Generated Power & Demand", yaxis_title="MW") ); 
          plot([water_lev],                      Layout(title="Reservoir", yaxis_title="hm3") ); 
          plot([water_spill,water_tur,water_inf],Layout(title="Water Inflows & Outflows", yaxis_title="hm3/h") ); 
          plot([price],                          Layout(title="Electrical Price", yaxis_title="€/MW", xaxis_title="hour of the day") );
        ]
end


