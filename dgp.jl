# need to generate Data's columns names
# function DGP(coef_estimation::Vector{Float64}; nofT=2, nofx=2, nofG=10, nofF=4, intercept_exist::Bool = true)
	
function DGP(; nofG=10, nofF=4, nofT=2,
               nofx = 1, # at least 1 which is constant
	   		   b0 = 0.5,
			   sig2uo = 0.1, sig2us = 0.2, 
			   sig2co = 0.1, sig2cs = 0.2, 
			   sig2wo = 0.1, sig2ws = 0.2,
			   σₓₒ² = 1.0) 

	nofobs = nofG*nofF*nofT

	β = b0*ones(nofx) # vector, (nofx x 1)


    X = [randn(nofobs, nofx-1)*sqrt(σₓₒ²)   ones(nofobs,1)]  # (N x nofx)


	W0 = rand(Normal(0, sqrt(sig2wo)), nofG) .- rand(TruncatedNormal(0, sqrt(sig2ws), 0, Inf), nofG)  # group-specific
    W  = repeat(W0, inner = nofF*nofT)

	C0 = rand(Normal(0, sqrt(sig2co)), nofF*nofG) .- rand(TruncatedNormal(0, sqrt(sig2cs), 0, Inf), nofF*nofG)  # firm-specific
    C  = repeat(C0, inner = nofT)

	U = rand(Normal(0, sqrt(sig2uo)), nofobs) .- rand(TruncatedNormal(0, sqrt(sig2us), 0, Inf), nofobs)  # obs-specific
	
	Tid = repeat(1:nofT,        outer = nofG*nofF)  # time var
	Gid = repeat(1:nofG,        inner = nofF*nofT)    # group id
	Fid = repeat(1:(nofG*nofF), inner = nofT)  # firm id

	#=
	get dummy variables to generate Data Y
	T = Matrix{Int64}(undef, nobs, nofT)
	for t in 1:nofT
        T[:,t] = Float64.(Tid .== t)  
    end 

	Y = (X*β) .+ (T*ϴ) .+ C .+ W .+ U
	=#

	Y = (X*β) .+ U .+ C .+ W


	Data = DataFrame(hcat(Y, X, Gid, Fid, Tid), :auto)

	
	x_names      = ["x$(i)t" for i in 1:nofx]
    x_names[end] = "_cons"

	Names = vcat("yit", x_names, "Group", "Firm", "Time") 
	rename!(Data, Symbol.(Names[:]))

	return Data
	 
end



