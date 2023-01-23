
using Plots
using JuMP
using Gurobi
using Distributions
using PiecewiseLinearOpt

#My verion of PiecewiseLinearOpt, it includes the function myPiecewiseLinearOpt, similar to piecewiselinear
include("my_lib.jl")

############################################  PiecewiseLinearOpt: HelloWorld Univariate ############################################

model = Model()
set_optimizer(model, Gurobi.Optimizer)
@variable(model, x)
d = 0:(pi/2):4pi
y_range = [sin(x) for x in d]
z = piecewiselinear(model, x, d, sin, method=:ZigZagInteger)

@objective(model, Max, z)


try
	rm("my_log_file.txt")
catch e
	println("Problem cleaning Gurobi Logs")
end

set_optimizer_attribute(model, "LogFile", "my_log_file.txt")
set_optimizer_attribute(model, "MIPGap", 0.01)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))

X = [JuMP.value(x)]
Z = [JuMP.value(z)]

@objective(model, Min, z)
optimize!(model)
X = push!(X,JuMP.value(x))
Z = push!(Z,JuMP.value(z))

#Ploting
if termination_status(model) == OPTIMAL
    plot(d,y_range)
    scatter!((X[1],Z[1]),label="Max")    
	scatter!((X[2],Z[2]),label="Min")    
end

############################################  PiecewiseLinearOpt: HelloWorld Bivariate ############################################




##################### PiecewiseLinearOpt and analytical


#f(x,y) = exp(x+y)
f(x,y) = 3*(1-x)^2*exp(-(x^2) - (y+1)^2)  - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2)  - 1/3*exp(-(x+1)^2 - y^2)

x_range = -4:1:4
y_range = -4:1:4
z_range = [ f(x,y) for x in x_range, y in y_range]

model = Model()
set_optimizer(model, Gurobi.Optimizer)
@variable(model, x)
@variable(model, y)
@variable(model, z)

z = piecewiselinear(model, x, y, x_range, y_range, (u,v) -> f(u,v))#,method=:ZigZagInteger)
#@constraint(model, x + y + 3 >= z)
@objective(model, Min, z)

#Solve 
set_optimizer_attribute(model, "MIPGap", 0.01)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))
println("")
println("piecewiselinear")
println("[" ,JuMP.value(x),";",JuMP.value(y),";",JuMP.value(z),"]")

println("analytical")
println("[" ,x_range[findmin(z_range)[2][1]],";",y_range[findmin(z_range)[2][2]],";",findmin(z_range)[1],"]")
#Ploting
if termination_status(model) == OPTIMAL
    surface(x_range,y_range,z_range)
    scatter!([JuMP.value(x)],[JuMP.value(y)],[JuMP.value(z)],label="min")
end


####################   myPiecewiseLinearOpt: ZZB 
model = Model()
set_optimizer(model, Gurobi.Optimizer)
@variable(model, x)
@variable(model, y)
@variable(model, z)
myPiecewiseLinearOpt(model, x, y, z, x_range, y_range, f, "ZZB", '1')
@constraint(model, x + y + 3 >= z)
@objective(model, Max, z)

#Solve 
set_optimizer_attribute(model, "MIPGap", 0.0001)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))
println("")

println("Mypiecewiselinear: ZZB")
println("[" ,JuMP.value(x),";",JuMP.value(y),";",JuMP.value(z),"]")


#####################   myPiecewiseLinearOpt: ZZI 
model = Model()
set_optimizer(model, Gurobi.Optimizer)
@variable(model, x)
@variable(model, y)
@variable(model, z)
myPiecewiseLinearOpt(model, x, y, z, x_range, y_range, f, "ZZI", '1')
@constraint(model, x + y + 3 >= z)
@objective(model, Max, z)

#Solve 
set_optimizer_attribute(model, "MIPGap", 0.0001)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))
println("")
println("Mypiecewiselinear: ZZI")
println("[" ,JuMP.value(x),";",JuMP.value(y),";",JuMP.value(z),"]")



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
set_optimizer_attribute(model, "MIPGap", 0.01)
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
set_optimizer_attribute(model, "MIPGap", 0.01)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))

#Ploting
if termination_status(model) == OPTIMAL
    push!(sol[1],JuMP.value(x))
    push!(sol[2],JuMP.value(y))
    push!(sol[3],JuMP.value(z))
end
surface(x_range,y_range,f)
surface(x_range,y_range,g,alpha=0.5)
scatter!(sol[1],sol[2],sol[3],label="min & max")
 

############################################  PiecewiseLinearOpt: 2 non linear functions ############################################

f(x,y) = 3*(1-x)^2*exp(-(x^2) - (y+1)^2)  - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2)  - 1/3*exp(-(x+1)^2 - y^2)
f(x,y) = x*y
g(x,y) = x*x + y

x_range = -4:1:4
y_range = -4:1:4

model = Model()
set_optimizer(model, Gurobi.Optimizer)
@variable(model, x)
@variable(model, y)
@variable(model, z)
@variable(model, x2)
@variable(model, y2)
@variable(model, z2)

z  = piecewiselinear(model, x, y, x_range, y_range, (u,v) -> f(u,v))
z2 = piecewiselinear(model, x2, y2, x_range, y_range, (u,v) -> g(u,v))

@objective(model, Max, z + z2)

#Solve 
set_optimizer_attribute(model, "MIPGap", 0.01)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))

println("[" ,JuMP.value(x),";",JuMP.value(y),";",JuMP.value(z),"]")
println("[" ,JuMP.value(x2),";",JuMP.value(y2),";",JuMP.value(z2),"]")


surface(x_range,y_range,f,alpha=1)
surface!(x_range,y_range,g,alpha=0.5)
scatter!([JuMP.value(x),JuMP.value(x2)],[JuMP.value(y),JuMP.value(y2)],[JuMP.value(z),JuMP.value(z2)],label="min & max")
