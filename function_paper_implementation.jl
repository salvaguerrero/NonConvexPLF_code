
using Plots
using JuMP
using Gurobi
using Distributions
using PiecewiseLinearOpt

include("my_lib.jl")

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
