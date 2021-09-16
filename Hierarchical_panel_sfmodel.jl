module Hierarchical_panel_sfmodel




export spec, init_vec, opt, 
       fit, predict,
       # likelihood functions 
       LL_T, 
       DGP, # Data generate process for test.jl
       # macros for spec() and init_vec(); 
       @is_intercept_exist, 
       @depvar, 
       @frontier, @timevar,
       @idvar, @Gvar,
       @sigma2_u_0, @σ²ᵤ₀, @sigma2_u_star, @σ²ᵤ_star, # for init_vec
       @sigma2_c_0, @σ²₍₀, @sigma2_c_star, @σ²₍_star, # for init_vec
       @sigma2_w_0, @σ²ω⁰, @sigma2_w_star, @σ²w_star, # for init_vec
       @eq,

       # functions for init_vec(); 
       is_intercept_exist,
       frontier, timevar,
       sigma2_u_0, σ²ᵤ₀, sigma2_u_star, σ²ᵤ_star, 
       sigma2_c_0, σ²₍₀, sigma2_c_star, σ²₍_star, 
       sigma2_w_0, σ²ω⁰, sigma2_w_star, σ²w_star,
       all_init,

       # functions for opt()
       warmstart_solver, warmstart_maxIT,
       main_solver, main_maxIT, tolerance, silent,
       # functions for sfmodel_fit
       useData

       
       

using Optim
using DataFrames
using NLSolversBase              # for hessian!
using StatsFuns                  # for normlogpdf(), normlogcdf()
using Statistics 
using StatsBase                  # for countmap()
using HypothesisTests            # for pvalue()
using LinearAlgebra              # extract diagnol and Matrix(I,...)
using Distributions              # for TDist, Normal
using DataStructures             # for OrderedDict
using PrettyTables
using Revise
using KahanSummation
using Kronecker                  # for loglikelyhood compuation
using Random                     # for loglikelyhood
using HaltonSequences

Random.seed!(0);                    # Fix random seed generator for reproducibility. Need Random package.

""" create  HaltonSequences  """

nofdraw = 100
randseq_tmp = HaltonPoint(4, start=20, length=nofdraw)  # dimension is 4(w_0, w_star, c_0, c_star)

randseq = zeros(nofdraw, 4)
for i in 1:nofdraw
    randseq[i,:] = randseq_tmp[i]
end

my_w_g_1 = quantile.(Normal(0,1), randseq[:,1])
my_w_g_2 =   quantile.(Normal(0,1), (0.5 .* randseq[:,2]) .+ 0.5)

my_c_i_1 = quantile.(Normal(0,1), randseq[:,3])
my_c_i_2 = quantile.(Normal(0,1), (0.5 .* randseq[:,4]) .+ 0.5)

################################################
##    include other files; order important    ##
################################################

include("macro.jl")
include("get_var.jl")
include("log_likelyhood.jl")
# include("getvar_simple.jl")
# include("log_likelyhood_simple.jl")
include("prediction.jl")

include("main.jl")
# include("main_simple.jl")
include("dgp.jl")



end # Hierarchical_panel_sfmodel module