using DataFrames, Random, Distributions

# need to generate Data's columns names
function DGP(coef_estimation::Vector{Float64}; noftime_perID=2, nofx=2, nofG=10, nofid_perG=4, intercept_exist::Bool = true)
	# nofx might include intercept
	# default require length of coef_estimation equal to 9
	nobs = nofG*nofid_perG*noftime_perID

	if intercept_exist
		(length(coef_estimation) != noftime_perID + nofx + 6 - 1) && throw("length of coef_estimation should equal to $( noftime_perID) + $(nofx) + 6 - 1, minus one to avoid dummy variable trap")
	else
		(length(coef_estimation) != noftime_perID + nofx + 6) && throw("length of coef_estimation should equal to $( noftime_perID) + $(nofx) + 6")
	end
	nofT = noftime_perID
	intercept_exist && (nofT = noftime_perID - 1 )

	β = coef_estimation[1:nofx]
	ϴ = coef_estimation[(nofx+1):(nofx+nofT)]

	# σᵤ₀, σᵤ_star, σ₍₀, σ₍_star, σw⁰, σw_star
	pos_beg = nofx + nofT + 1
	pos_end = nofx + nofT + 6
	σᵤ₀, σᵤ_star, σ₍₀, σ₍_star, σw⁰, σw_star = sqrt.(coef_estimation[pos_beg:pos_end])

    
	if intercept_exist
		X = hcat(randn(Float64, (nobs, nofx-1)), ones(nobs))
	else
		X = randn(Float64, (nobs, nofx))
	end

	
	
	T_vec = repeat(1:noftime_perID, outer = nofG*nofid_perG)
	# get dummy variables to generate Data Y
	T = Matrix{Int64}(undef, nobs, nofT)
	for t in 1:nofT
        T[:,t] = Float64.(T_vec .== t)  
    end 
    

	G = repeat(1:nofG, inner = nofid_perG*noftime_perID)
	id = repeat(1:(nofG*nofid_perG), inner = noftime_perID)

	# DGP for error term
	U = rand(Normal(0, σᵤ₀), nobs) .- rand(TruncatedNormal(0, σᵤ_star, 0, Inf), nobs)
	C = rand(Normal(0, σ₍₀), nobs) .- rand(TruncatedNormal(0, σ₍_star, 0, Inf), nobs)
	W = rand(Normal(0, σw⁰), nobs) .- rand(TruncatedNormal(0, σw_star, 0, Inf), nobs)

	Y = (X*β) .+ (T*ϴ) .+ C .+ W .+ U

	Data = DataFrame(hcat(Y, X, T_vec, G, id))

	
	x_names = ["x$(i)t" for i in 1:nofx]
	intercept_exist && (x_names[end] = "_cons")

	Names = vcat("yit", x_names, "Time", "Group", "id") 
	rename!(Data, Symbol.(Names[:]))

	return Data

	 
end

# println(DGP([1.0,1.0,2.2, 0.5,0.5,0.5,0.5,0.5,0.5]))
println(DGP([1.0,1.0,2.2, 0.43,0.5,0.1,0.7,0.9,0.6], nofG=100))



