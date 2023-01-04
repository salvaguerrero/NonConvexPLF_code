
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


w = 3
I = 1:(1+2^w)
J = 1:(1+2^w)
x_range = -4:1:4
y_range = -4:1:4



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

# println("S1: ",S1)
# println("S2: ",S2)
# println("S3: ",S3)
# println("S4: ",S4)

###########  Non linear function        ################## 
# I = -1:0.5:1
# J = -1:0.5:1

Q,V = mesh(x_range,y_range)	
P(x,y) = 3*(1-x)^2*exp(-(x^2) - (y+1)^2)  - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2)  - 1/3*exp(-(x+1)^2 - y^2)


#Q->X->i->m
#v->Y->J->n
q_hat = x_range
v_hat = y_range
p_hat = P.(Q,V)

###########  Binary Zig-Zag Based Modeling        ################## 

C = Dict()
C = C_matrix(r)

#ZZB: 1
#ZZI: 2
model_type = 2

model = Model(Gurobi.Optimizer)
@variable(model, q)
@variable(model, v)
@variable(model, p)
@variable(model, ø[i in 1:m,j in 1:n] >= 0)
if model_type == 1
    @variable(model, ζ_k[k in 1:r], Bin )       
    @variable(model, ζ_l[l in 1:s], Bin )
else
    @variable(model, ζ_k[k in 1:r], Int )       
    @variable(model, ζ_l[l in 1:s], Int )
end
#triangle selection variables
@variable(model, z1, Bin )
@variable(model, z2, Bin )



#Binary Zig-Zag constraints
@constraint(model, q == sum(q_hat[i]*ø[i,j]   for i in 1:m, j in 1:n) )   #13_a1
@constraint(model, v == sum(v_hat[j]*ø[i,j]   for i in 1:m, j in 1:n) )   #13_a2
@constraint(model, p == sum(p_hat[i,j]*ø[i,j] for i in 1:m, j in 1:n) )   #13_b
@constraint(model, 1 == sum(ø[i,j]            for i in 1:m, j in 1:n) )   #13_c

if model_type == 1
    # ZZB
    @constraint(model,tr_d_1[k in 1:r], sum( C[r][ (i == 1 ? 1   : i-1 ) ,k]*sum(ø[i,j] for j in 1:n) for i in 1:m) <= ζ_k[k] + sum(ζ_k[l]*2^(l-k-1) for l in (k+1):r) )      #13_d
    @constraint(model,tr_d_2[k in 1:r], sum( C[r][ (i == m ? m-1 : i   ) ,k]*sum(ø[i,j] for j in 1:n) for i in 1:m) >= ζ_k[k] + sum(ζ_k[l]*2^(l-k-1) for l in (k+1):r) )      #13_d

    @constraint(model,tr_e_1[l in 1:s], sum( C[s][ (j == 1 ? 1   : j-1 ) ,l]*sum(ø[i,j] for i in 1:m) for j in 1:n) <= ζ_l[l] + sum(ζ_l[ll]*2^(ll-l-1) for ll in (l+1):r))      #13_e
    @constraint(model,tr_e_2[l in 1:s], sum( C[s][ (j == n ? n-1 : j   ) ,l]*sum(ø[i,j] for i in 1:m) for j in 1:n) >= ζ_l[l] + sum(ζ_l[ll]*2^(ll-l-1) for ll in (l+1):r))      #13_e
else
    # #ZZI
    @constraint(model,tr_d_1[k in 1:r], sum( C[r][ (i == 1 ? 1   : i-1 ) ,k]*sum(ø[i,j] for j in 1:n) for i in 1:m) <= ζ_k[k])      #17_a
    @constraint(model,tr_d_2[k in 1:r], sum( C[r][ (i == m ? m-1 : i   ) ,k]*sum(ø[i,j] for j in 1:n) for i in 1:m) >= ζ_k[k])      #17_a

    @constraint(model,tr_e_1[l in 1:s], sum( C[s][ (j == 1 ? 1   : j-1 ) ,l]*sum(ø[i,j] for i in 1:m) for j in 1:n) <= ζ_l[l])      #17_b
    @constraint(model,tr_e_2[l in 1:s], sum( C[s][ (j == n ? n-1 : j   ) ,l]*sum(ø[i,j] for i in 1:m) for j in 1:n) >= ζ_l[l])      #17_b
end


#triangle selection constraints (Equations 15)
@constraint(model, sum(ø[a[1],a[2]]   for a in S1) <=     z1  )
@constraint(model, sum(ø[a[1],a[2]]   for a in S2) <= 1 - z1  )
@constraint(model, sum(ø[a[1],a[2]]   for a in S3) <=     z2  )
@constraint(model, sum(ø[a[1],a[2]]   for a in S4) <= 1 - z2  )


#@constraint(model,p <= q + v + 5)


@objective(model, Min, p)

set_optimizer_attribute(model, "MIPGap", 0.1)
println("Running optimization process:")
optimize!(model)
println("    Execution status: ",termination_status(model))
println("Model min: " ,JuMP.value(p))
println("Analytical min: ",minimum(p_hat))
println("Analytical max: ",maximum(p_hat))
println("[" ,JuMP.value(p),";",JuMP.value(v),";",JuMP.value(q),"]")

