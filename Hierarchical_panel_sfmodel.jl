module Hierarchical_panel_sfmodel


#! not specify σ²ᵤ₀, σ²ᵤ_star, σ²₍₀, σ²₍_star, σ²w⁰, σ²w_star since they are only parameterized by constant
export spec, init_vec, opt, 
       fit, predict,
       # likelihood functions 
       LL_T, 
       # macros for spec() and init_vec(); 
       @depvar, 
       @frontier, @timevar,
       @idvar, @Gvar,
       @sigma2_u_0, @σ²ᵤ₀, @sigma2_u_star, @σ²ᵤ_star, # for init_vec
       @sigma2_c_0, @σ²₍₀, @sigma2_c_star, @σ²₍_star, # for init_vec
       @sigma2_w_0, @σ²ω⁰, @sigma2_w_star, @σ²w_star, # for init_vec
       @eq,

       # functions for init_vec(); 
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
using Statistics                 #! not sure
using HypothesisTests            # for pvalue()
using LinearAlgebra              # extract diagnol and Matrix(I,...)
using Distributions              # for TDist, Normal
using DataStructures             # for OrderedDict
using PrettyTables
using Revise
using KahanSummation
using Kronecker                  # for loglikelyhood compuation
using Random                     # for loglikelyhood

# Random.seed!(0);                    # Fix random seed generator for reproducibility. Need Random package.
# already set random seed in model_loglikelyhood.jl, therefore not set random seed here

################################################
##    include other files; order important    ##
################################################

include("matrix_computation.jl")
include("macro.jl")
include("hierachy-model_loglikelyhood.jl")
include("getvar.jl")
include("prediction.jl")
include("main_function.jl")



end # Hierarchical_panel_sfmodel module