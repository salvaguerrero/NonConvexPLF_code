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

function myPiecewiseLinearOpt(model, x, y, z, x_range, y_range, f, model_type, id)

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


	ø = @variable(model, [i in 1:m,j in 1:n], base_name = "teta"*id)
	@constraint(model,[i in 1:m,j in 1:n],ø[i,j]  >= 0 )

	if model_type == "ZZI"
		ζ_k = @variable(model, [k in 1:r], Int, base_name = "ζ_k_"*id)       
		ζ_l = @variable(model, [l in 1:s], Int, base_name = "ζ_l_"*id)
	else
		ζ_k = @variable(model, [k in 1:r], binary = true, base_name = "ζ_k_"*id)       
		ζ_l = @variable(model, [l in 1:s], binary = true, base_name = "ζ_l_"*id)
	end
	#triangle selection variables
	z1 = @variable(model, binary = true, base_name = "z2:"*id)
	z2 = @variable(model, binary = true, base_name = "z1_"*id)

	#Binary Zig-Zag constraints
	@constraint(model, x == sum(x_hat[i]*ø[i,j]   for i in 1:m, j in 1:n) )   #13_a1
	@constraint(model, y == sum(y_hat[j]*ø[i,j]   for i in 1:m, j in 1:n) )   #13_a2
	@constraint(model, z == sum(z_hat[i,j]*ø[i,j] for i in 1:m, j in 1:n) )   #13_b
	@constraint(model, 1 == sum(ø[i,j]            for i in 1:m, j in 1:n) )   #13_c

	if model_type == "ZZI"

		#ZZI
		@constraint(model,[k in 1:r], sum( C[r][ (i == 1 ? 1   : i-1 ) ,k]*sum(ø[i,j] for j in 1:n) for i in 1:m) <= ζ_k[k])      #17_a
		@constraint(model,[k in 1:r], sum( C[r][ (i == m ? m-1 : i   ) ,k]*sum(ø[i,j] for j in 1:n) for i in 1:m) >= ζ_k[k])      #17_a

		@constraint(model,[l in 1:s], sum( C[s][ (j == 1 ? 1   : j-1 ) ,l]*sum(ø[i,j] for i in 1:m) for j in 1:n) <= ζ_l[l])      #17_b
		@constraint(model,[l in 1:s], sum( C[s][ (j == n ? n-1 : j   ) ,l]*sum(ø[i,j] for i in 1:m) for j in 1:n) >= ζ_l[l])      #17_b
	else

		# ZZB
		@constraint(model,[k in 1:r], sum( C[r][ (i == 1 ? 1   : i-1 ) ,k]*sum(ø[i,j] for j in 1:n) for i in 1:m) <= ζ_k[k] + sum(ζ_k[l]*2^(l-k-1) for l in (k+1):r) )      #13_d
		@constraint(model,[k in 1:r], sum( C[r][ (i == m ? m-1 : i   ) ,k]*sum(ø[i,j] for j in 1:n) for i in 1:m) >= ζ_k[k] + sum(ζ_k[l]*2^(l-k-1) for l in (k+1):r) )      #13_d

		@constraint(model,[l in 1:s], sum( C[s][ (j == 1 ? 1   : j-1 ) ,l]*sum(ø[i,j] for i in 1:m) for j in 1:n) <= ζ_l[l] + sum(ζ_l[ll]*2^(ll-l-1) for ll in (l+1):r))      #13_e
		@constraint(model,[l in 1:s], sum( C[s][ (j == n ? n-1 : j   ) ,l]*sum(ø[i,j] for i in 1:m) for j in 1:n) >= ζ_l[l] + sum(ζ_l[ll]*2^(ll-l-1) for ll in (l+1):r))      #13_e
	
	end


	#triangle selection constraints (Equations 15)
	@constraint(model, sum(ø[a[1],a[2]]   for a in S1) <=     z1  )
	@constraint(model, sum(ø[a[1],a[2]]   for a in S2) <= 1 - z1  )
	@constraint(model, sum(ø[a[1],a[2]]   for a in S3) <=     z2  )
	@constraint(model, sum(ø[a[1],a[2]]   for a in S4) <= 1 - z2  )


end
