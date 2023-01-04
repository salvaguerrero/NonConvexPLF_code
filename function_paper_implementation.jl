
using Plots
using JuMP
using Gurobi
using Distributions
using PiecewiseLinearOpt

############################################  AUX Functions  ############################################
function congruent_modulo(a,b,n)
    #return true/false if a ≡ b (mod n)
    return mod(a,n)==mod(b,n)
    #     for k in -10:10
    #         if a == k*n+b 
    #             return true
    #         end
    #     end
    #     return false
end

function mesh(I,J)
	# input: I = -1:0.5:1
	#		 J = -1:0.5:1
	X = zeros( length(I),length(J) )
	Y = zeros( length(I),length(J) )
	ii = 1
	jj = 1
	for i in I
		for j in J
			Y[ii,jj] = i
			X[ii,jj] = j
			jj = jj + 1
		end
		ii = ii + 1
		jj = 1
	end
	return X,Y
end

function S_1(I,J)
    S1 = []
	for i in I
        for j in J
            if congruent_modulo(i,j,2) && congruent_modulo(i+j,2,4)
                push!(S1,(i,j))
            end
        end
    end
    return S1
end

function S_2(I,J)
    S2 = []
	for i in I
        for j in J
            if congruent_modulo(i,j,2) && congruent_modulo(i+j,0,4)
                push!(S2,(i,j))
            end
        end
    end
    return S2
end

function S_3(I,J)
    S3 = []
	for i in I
        for j in J
            if !congruent_modulo(i,j,2) && congruent_modulo(i+j,3,4)
                push!(S3,(i,j))
            end
        end
    end
    return S3
end

function S_4(I,J)
	S4 = []
    for i in I
        for j in J
            if !congruent_modulo(i,j,2) && congruent_modulo(i+j,1,4)
                push!(S4,(i,j))
            end
        end
    end
    return S4
end

function C_matrix(R)
    
    C = Dict()
    C[1] = [0;1]
    r = 1
    for r in 1:Int(R)
        C[r+1] = [ 
                   C[r]                              zeros(2^r,1);

                   C[r]+ones(2^r,1)*C[r][2^r,:]'      ones(2^r,1) ]
        # println("theoretical size: ",(2^(r+1),r+1)," Actual size ", size( C[r+1] ))
        # display(C[r])
    end

    return C
end

 
function myPiecewiseLinearOpt(model, x, y, z, x_range, y_range, f, model_type)

	I = 1:length(x_range)
	J = 1:length(y_range)

	m = length(I)
	n = length(J)

	d = length(I)-1 
	r = Int(log2(m-1))
	s = Int(log2(n-1))


	###########          Triangular Selection          ##################
	S1 = []
	S2 = []
	S3 = []
	S4 = []
	S1 = S_1(I,J)
	S2 = S_2(I,J)
	S3 = S_3(I,J)
	S4 = S_4(I,J)

	###########  Non linear function        ################## 
	Q,V = mesh(x_range,y_range)	
	x_hat = x_range
	y_hat = y_range 
	z_hat = f.(Q,V)

	###########  Zig-Zag Based Modeling        ################## 

	C = Dict()
	C = C_matrix(r)


	@variable(model, ø[i in 1:m,j in 1:n] >= 0)
	if model_type == "ZZI"
		@variable(model, ζ_k[k in 1:r], Int )       
		@variable(model, ζ_l[l in 1:s], Int )
	else
		@variable(model, ζ_k[k in 1:r], Bin )       
		@variable(model, ζ_l[l in 1:s], Bin )
	end
	#triangle selection variables
	@variable(model, z1, Bin )
	@variable(model, z2, Bin )



	#Binary Zig-Zag constraints
	@constraint(model, x == sum(x_hat[i]*ø[i,j]   for i in 1:m, j in 1:n) )   #13_a1
	@constraint(model, y == sum(y_hat[j]*ø[i,j]   for i in 1:m, j in 1:n) )   #13_a2
	@constraint(model, z == sum(z_hat[i,j]*ø[i,j] for i in 1:m, j in 1:n) )   #13_b
	@constraint(model, 1 == sum(ø[i,j]            for i in 1:m, j in 1:n) )   #13_c

	if model_type == "ZZI"

		#ZZI
		@constraint(model,tr_d_1[k in 1:r], sum( C[r][ (i == 1 ? 1   : i-1 ) ,k]*sum(ø[i,j] for j in 1:n) for i in 1:m) <= ζ_k[k])      #17_a
		@constraint(model,tr_d_2[k in 1:r], sum( C[r][ (i == m ? m-1 : i   ) ,k]*sum(ø[i,j] for j in 1:n) for i in 1:m) >= ζ_k[k])      #17_a

		@constraint(model,tr_e_1[l in 1:s], sum( C[s][ (j == 1 ? 1   : j-1 ) ,l]*sum(ø[i,j] for i in 1:m) for j in 1:n) <= ζ_l[l])      #17_b
		@constraint(model,tr_e_2[l in 1:s], sum( C[s][ (j == n ? n-1 : j   ) ,l]*sum(ø[i,j] for i in 1:m) for j in 1:n) >= ζ_l[l])      #17_b
	else

		# ZZB
		@constraint(model,tr_d_1[k in 1:r], sum( C[r][ (i == 1 ? 1   : i-1 ) ,k]*sum(ø[i,j] for j in 1:n) for i in 1:m) <= ζ_k[k] + sum(ζ_k[l]*2^(l-k-1) for l in (k+1):r) )      #13_d
		@constraint(model,tr_d_2[k in 1:r], sum( C[r][ (i == m ? m-1 : i   ) ,k]*sum(ø[i,j] for j in 1:n) for i in 1:m) >= ζ_k[k] + sum(ζ_k[l]*2^(l-k-1) for l in (k+1):r) )      #13_d

		@constraint(model,tr_e_1[l in 1:s], sum( C[s][ (j == 1 ? 1   : j-1 ) ,l]*sum(ø[i,j] for i in 1:m) for j in 1:n) <= ζ_l[l] + sum(ζ_l[ll]*2^(ll-l-1) for ll in (l+1):r))      #13_e
		@constraint(model,tr_e_2[l in 1:s], sum( C[s][ (j == n ? n-1 : j   ) ,l]*sum(ø[i,j] for i in 1:m) for j in 1:n) >= ζ_l[l] + sum(ζ_l[ll]*2^(ll-l-1) for ll in (l+1):r))      #13_e
	
	end


	#triangle selection constraints (Equations 15)
	@constraint(model, sum(ø[a[1],a[2]]   for a in S1) <=     z1  )
	@constraint(model, sum(ø[a[1],a[2]]   for a in S2) <= 1 - z1  )
	@constraint(model, sum(ø[a[1],a[2]]   for a in S3) <=     z2  )
	@constraint(model, sum(ø[a[1],a[2]]   for a in S4) <= 1 - z2  )


end


############################################  PiecewiseLinearOpt: Bivariate ############################################

x_range = -4:1:4
y_range = -4:1:4
f(x,y) = 3*(1-x)^2*exp(-(x^2) - (y+1)^2)  - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2)  - 1/3*exp(-(x+1)^2 - y^2)
f(x,y) = exp(x+y)
# be aware!!!!!!!!!!   length(x_range) == (1+2^w)       
Q,V = mesh(x_range,y_range)	
x_hat = x_range
y_hat = y_range 
z_hat = f.(Q,V)

model = Model(Gurobi.Optimizer)
@variable(model, x)
@variable(model, y)
@variable(model, z)
#myPiecewiseLinearOpt(model, x, y, z, x_range, y_range, f, "ZZI")
z = piecewiselinear(model, x, y, x_range, y_range, (u,v) -> f(u,v))#,method=:ZigZag)

@objective(model, Max, z)

# Min
@objective(model, Min, z)

set_optimizer_attribute(model, "MIPGap", 0.1)
#set_optimizer_attribute(model, "LogFile", "my_log_file.txt")
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))

println("Min:    [" ,JuMP.value(x),";",JuMP.value(y),";",JuMP.value(z),"]")
minimum(z_hat)



############################################  Hydro Plant ############################################

model = Model()
set_optimizer(model, Gurobi.Optimizer)
@variable(model, g_gen   >= 0)
@variable(model, g_gen_h >= 0)
@variable(model, w_tur   >= 0)
@variable(model, w_lev_0 >= 0)
@variable(model, w_lev_f >= 0)
@variable(model, w_spi   >= 0)

@objective(model, Min, g_gen*25 )

@constraint(model, g_gen_h + g_gen == 120 ) 

@constraint(model,  w_lev_f == w_lev_0 - w_tur + 20 - w_spi )
@constraint(model, w_lev_0 == 10 )
@constraint(model, g_gen_h <= 100 )

p(x,y) = 10*x*y

q_range = 0:1:128
v_range = 0:1:128
Q,V = mesh(q_range,v_range)	
z_hat = p.(Q,V)

#g_gen_h = piecewiselinear(model, w_tur, w_lev_f, q_range, v_range, p)
myPiecewiseLinearOpt(model, w_tur, w_lev_f, g_gen_h, q_range, v_range, p, "ZZI")
# be aware!!!!!!!!!!   length(x_range) == (1+2^w)       


set_optimizer_attribute(model, "MIPGap", 0.1)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))

JuMP.value(g_gen_h) 
JuMP.value(g_gen)
JuMP.value(w_tur)
JuMP.value(w_spi)
JuMP.value(w_lev_f)




