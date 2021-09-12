using CSV, DataFrames
using Statistics, StatsFuns
using LinearAlgebra,  NLSolversBase, Optim #,SparseArrays
using Revise, JLD2
using BenchmarkTools
using DataFrames
using Random
using Test
using StatFiles
using HaltonSequences
using Hierarchical_panel_sfmodel



# only intercept data & all_init initial value:

# df = CSV.read("sim_data1_only_intercept.csv", DataFrame; header=1, delim=",")
# df[!,:_cons] .= 1.0
# println(typeof(df))

# spec(@depvar(yit), @frontier(_cons), @timevar(time_id), @idvar(firm_id), @Gvar(group_id))
# println("spec function is ok")

# init_vec(all_init(1.0))
# opt( warmstart_solver(),    #* empty means no warmstart 
#      warmstart_maxIT(10),
# 	 main_solver(Newton()),      #* BFGS does not work
# 	 main_maxIT(2000), 
# 	 tolerance(1e-8)
# 	 #, silent(true)
#      )

# println("opt function is ok")
# res = ()

# @time res = fit(df)

# println("start to predict !!")
# println("estimation of frontier of first two observation:", predict(@eq(frontier), df)[1], predict(@eq(frontier), df)[1])
# println("estimation of log_σ²₍₀", predict(@eq(log_σ²₍₀), df))
# println("estimation of σ²₍₀", predict(@eq(σ²₍₀), df)) 

# --------------------------------------------------------------------------------------------- #
# df = CSV.read("sim_data1.csv", DataFrame; header=1, delim=",")
# df[!,:_cons] .= 1.0

# spec(@is_intercept_exist(true),@depvar(yit), @frontier(_cons,xit), @timevar(time_id), @idvar(firm_id), @Gvar(group_id))

# # @timebar() require user to set init_vec for every period, therefore user should know how many period span 
# # this data have two unique time , however we just give one intial elelment for @timevar in init_vec since we already have intercept
# init_vec(@frontier(0.5,1), @timevar(0.1), @σ²ᵤ₀(2), @σ²ᵤ_star(2), @σ²₍₀(0.75), @σ²₍_star(0.5), @σ²ω⁰(0.4), @σ²w_star(0.3))  
# opt( warmstart_solver(),    
#      warmstart_maxIT(10),
# 	 main_solver(Newton()),      
# 	 main_maxIT(2000), 
# 	 tolerance(1e-8)
# 	 #, silent(true)
#      )

# res = ()

# @time res = fit(df)

# println("start to predict !!")
# println("estimation of frontier of first two observation:", predict(@eq(frontier), df)[1], predict(@eq(frontier), df)[1])
# println("estimation of log_σ²₍₀", predict(@eq(log_σ²₍₀), df))
# println("estimation of σ²₍₀", predict(@eq(σ²₍₀), df)) 



# ------------------------------------------------------------------------- #

# df = CSV.read("sim_data1.csv", DataFrame; header=1, delim=",")


# spec(@is_intercept_exist(false), @depvar(yit), @frontier(xit), @timevar(time_id), @idvar(firm_id), @Gvar(group_id))
# # spec(@depvar(yit), @frontier(_cons,xit), @timevar(time_id), @idvar(firm_id), @Gvar(group_id))

# # @timebar() require user to set init_vec for every period, therefore user should know how many period span 
# init_vec(@frontier(0.5,1), @timevar(0.1))

# # init_vec(@frontier(0.5,1), @timevar(0.1, 0.1)) 

# opt( warmstart_solver(NelderMead()),    
#      warmstart_maxIT(10),
# 	 main_solver(NewtonTrustRegion()),      
# 	 main_maxIT(2000), 
# 	 tolerance(1e-8)
# 	 #, silent(true)
#      )
# res = ()

# @time res = fit(df)

# println("start to predict !!")
# println("estimation of frontier of first two observation:", predict(@eq(frontier), df)[1], predict(@eq(frontier), df)[1])
# println("estimation of log_σ²₍₀", predict(@eq(log_σ²₍₀), df))
# println("estimation of σ²₍₀", predict(@eq(σ²₍₀), df)) 

# ---------------------------------------------------------

df = DGP([1.0,1.0,2.2, 0.43,0.5,0.1,0.7,0.9,0.6], nofG=30, nofid_perG=5, noftime_perID=2, nofx=2)



spec(@is_intercept_exist(true), @depvar(yit), @frontier(_cons,x1t), @timevar(Time), @idvar(id), @Gvar(Group))
# spec(@depvar(yit), @frontier(_cons,xit), @timevar(time_id), @idvar(firm_id), @Gvar(group_id))

# @timebar() require user to set init_vec for every period, therefore user should know how many period span 
init_vec( @σ²ᵤ₀(2), @σ²ᵤ_star(2), @σ²₍₀(0.75), @σ²₍_star(0.5), @σ²ω⁰(0.4), @σ²w_star(0.3))

# init_vec(@frontier(0.5,1), @timevar(0.1, 0.1)) 

opt( warmstart_solver(NelderMead()),    
     warmstart_maxIT(10),
	 main_solver(NewtonTrustRegion()),      
	 main_maxIT(2000), 
	 tolerance(1e-8)
	 #, silent(true)
     )
res = ()

@time res = fit(df)


println("start to predict !!")
println("estimation of frontier of first two observation:", predict(@eq(frontier), df)[1], predict(@eq(frontier), df)[1])
println("estimation of log_σ²₍₀", predict(@eq(log_σ²₍₀), df))
println("estimation of σ²₍₀", predict(@eq(σ²₍₀), df)) 



