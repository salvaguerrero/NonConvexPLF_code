using PlotlyJS
using JuMP
using Gurobi
using PiecewiseLinearOpt
using Dates

include("my_lib.jl")

#################################### Small Hydro Case Study: Non Linear modeling of gravity
lib = "my_lib"

#Data------------------------------------------>
#Generators
Opex  = 25        #[€/MWh] Hydro, Gas
L_lim = 0	          #[MW] 
U_lim = 250			  #[MW]

#Demand
Time = 1:24
T_num = length(Time)
dem = [200*sin(t*(2*pi)/24)+300 for t in Time] 

#Hydro Plant
#Water levels [#hm3]
w_lev_0 = 10
w_lev_f = 10
w_apo   = 20#aportaciones [m3/h]
w_lev_M = 128 
w_lev_m = 0 
#Turbined water [#hm3]
w_tur_M = 128  
w_tur_m = 0    
# cubic reservoir a x b [m]
a = 10        
b = 10
# efficiency
eff = 0.6    

model = Model(Gurobi.Optimizer)
@variables(model, begin
    L_lim  <= g_gen[Time] <= U_lim
    w_tur[Time] >= 0
    w_lev[vcat(0, Time)] >= 0
    w_spi[Time] >= 0
end)
fix(w_lev[0],   w_lev_0; force = true)
fix(w_lev[end], w_lev_f; force = true)

p(q,v) = 9.81*q*v#/a/b*eff

q_range = w_tur_m:1:w_tur_M #0:128
v_range = w_lev_m:1:w_lev_M #0:128
if lib != "my_lib"
    g_gen_h = [
        piecewiselinear(model, w_tur[t] , w_lev[t], v_range, q_range, p,method=:ZigZag)#ZigZagInteger
        for t in Time
    ]
else
    @variable(model,g_gen_h[Time])
    for t in Time
	    myPiecewiseLinearOpt(model, w_tur[t] , w_lev[t], g_gen_h[t], v_range, q_range, p, "ZZB", string(t))
    end
end
set_lower_bound.(g_gen_h, L_lim)
set_upper_bound.(g_gen_h, U_lim)
@constraints(model, begin
    [t in Time], g_gen_h[t] + g_gen[t] == dem[t]
    [t in Time],  w_lev[t] == w_lev[t-1] + w_apo - (w_tur[t] + w_spi[t])
end) 

@objective(model, Min, Opex * sum(g_gen[t] for t in Time))

try
	rm("my_log_file.txt")
catch e
	println("Problem cleaning Gurobi Logs")
end
set_optimizer_attribute(model, "LogFile", "my_log_file.txt")
set_optimizer_attribute(model, "MIPGap", 0.01)
println("Running optimization process:")
t = DateTime(Dates.now())
optimize!(model)
println("    Execution status: ",termination_status(model))
println("Soving total time: ",DateTime(Dates.now())-t)

Gen = Dict()
Gen[1] = zeros(T_num)
Gen[2] = zeros(T_num)
Level  = zeros(T_num)
Turbined = zeros(T_num)
STMP = zeros(T_num)
Spilled = zeros(T_num)
Inflows = w_apo*ones(T_num)

Lin_e = zeros(T_num)

if termination_status(model) == OPTIMAL
    for t in Time
        Level[t]    = round(JuMP.value(w_lev[t])   , digits = 2)
        Turbined[t] = round(JuMP.value(w_tur[t])   , digits = 2)
        Spilled[t]  = round(JuMP.value(w_spi[t])   , digits = 2)
        Gen[1][t]   = round(JuMP.value(g_gen_h[t]) , digits = 2)
        Gen[2][t]   = round(JuMP.value(g_gen[t]) , digits = 2)
        STMP[t]     = 0#round(dual(Demand[t])        , digits = 2)
        Lin_e[t]    = abs(9.81*Turbined[t]*Level[t] - Gen[1][t])
    end
    println( "max linearization absolute erro: ", maximum(Lin_e))
   
    # Model Plotting
    g1_power = scatter(;x=Time, y=Gen[1], mode="lines+markers",name="Hydro",stackgroup="one")
    g2_power = scatter(;x=Time, y=Gen[2], mode="lines+markers",name="Gas"  ,stackgroup="one")
    demand   = scatter(;x=Time, y=dem,    mode="markers"      ,name="Demand")
  
    water_lev   = scatter(;x=Time, y=Level,   mode="lines+markers",name="Water Level")

    water_inf   = scatter(;x=Time, y=Inflows,mode="lines+markers" ,name="Water Inflows")
    water_tur   = scatter(;x=Time, y=Turbined,mode="lines+markers",name="Turbined water")
    water_spill = scatter(;x=Time, y=Spilled,mode="lines+markers" ,name="Spilled water")

    #price       = scatter(;x=Time, y=STMP,mode="lines+markers",name="STMP")

    p3 = [plot([g1_power,g2_power,demand],       Layout(title="Generated Power & Demand", yaxis_title="MW") ); 
          plot([water_lev],                      Layout(title="Reservoir", yaxis_title="hm3") ); 
          plot([water_spill,water_tur,water_inf],Layout(title="Water Inflows & Outflows", yaxis_title="hm3/h") ); 
          #plot([price],                          Layout(title="Electrical Price", yaxis_title="€/MW", xaxis_title="hour of the day") );
        ]
end


