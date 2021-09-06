###################################################################
##  macro and functions for spec(), init_vec(), and fit()  ##
###################################################################


#? ---- macros for spec() and functions for init_vec() ------------
 
 macro depvar(arg::Vararg)
     return (:depvar, collect(:($(arg))))
 end
 
  # -------------------
 
 macro frontier(arg::Vararg)   # for spec()          
     return (:frontier, collect(:($(arg))))
 end 
 
 function frontier(arg::Vararg)    # for init_vec(0.1, 0.1)
     return (:frontier, collect(:($(arg))) )
 end
 
 
 function frontier(arg::Vector)    # for init_vec(b0)
     return (:frontier,  arg)
 end
 
 
 
  # -------------------
 
 macro sigma2_u_0(arg::Vararg)
     return (:σ²ᵤ₀, collect(:($(arg))))
 end

 
 
 function sigma2_u_0(arg::Vararg)
     return (:σ²ᵤ₀, collect(:($(arg))))
 end
 
 function sigma2_u_0(arg::Vector)
     return (:σ²ᵤ₀, arg)
 end

  # -------------------
 
 macro σ²ᵤ₀(arg::Vararg)
     return (:σ²ᵤ₀, collect(:($(arg))))
 end

 
 
 function σ²ᵤ₀(arg::Vararg)
     return (:σ²ᵤ₀, collect(:($(arg))))
 end
 
 function σ²ᵤ₀(arg::Vector)
     return (:σ²ᵤ₀, arg)
 end
 

  #-------------------
  macro sigma2_u_star(arg::Vararg)
     return ( :σ²ᵤ_star, collect(:($(arg))))
 end

 
 function sigma2_u_star(arg::Vararg)
     return ( :σ²ᵤ_star, collect(:($(arg))))
 end
 
 function sigma2_u_star(arg::Vector)
     return ( :σ²ᵤ_star, arg)
 end

  #-------------------
  macro σ²ᵤ_star(arg::Vararg)
     return ( :σ²ᵤ_star, collect(:($(arg))))
 end

 
 function σ²ᵤ_star(arg::Vararg)
     return ( :σ²ᵤ_star, collect(:($(arg))))
 end
 
 function σ²ᵤ_star(arg::Vector)
     return ( :σ²ᵤ_star, arg)
 end

  # -------------------
  macro sigma2_c_0(arg::Vararg)
     return (:σ²₍₀, collect(:($(arg))))
 end

    
 function sigma2_c_0(arg::Vararg)
     return (:σ₍₀², collect(:($(arg))))
 end
 
 function sigma2_c_0(arg::Vector)
     return (:σ²₍₀, arg)
 end

 # -------------------
  macro σ²₍₀(arg::Vararg)
     return (:σ²₍₀, collect(:($(arg))))
 end
 
 function σ²₍₀(arg::Vararg)
     return (:σ²₍₀, collect(:($(arg))))
 end
 
 function σ²₍₀(arg::Vector)
     return (:σ²₍₀, arg)
 end


  # -------------------
  macro sigma2_c_star(arg::Vararg)
     return (:σ²₍_star, collect(:($(arg))))
 end
 
 function sigma2_c_star(arg::Vararg)
     return (:σ²₍_star, collect(:($(arg))))
 end
 
 function sigma2_c_star(arg::Vector)
     return (:σ²₍_star, arg)
 end

   # -------------------
  macro σ²₍_star(arg::Vararg)
     return (:σ²₍_star, collect(:($(arg))))
 end
 
 function σ²₍_star(arg::Vararg)
     return (:σ²₍_star, collect(:($(arg))))
 end
 
 function σ²₍_star(arg::Vector)
     return (:σ²₍_star, arg)
 end
# -------------------
  macro sigma2_w_0(arg::Vararg)
     return (:σ²w⁰, collect(:($(arg))))
 end
 
 
 function sigma2_w_0(arg::Vararg)
     return (:σ²w⁰, collect(:($(arg))))
 end
 
 function sigma2_w_0(arg::Vector)
     return (:σ²w⁰, arg)
 end

 # -------------------
  macro σ²ω⁰(arg::Vararg)
     return (:σ²w⁰, collect(:($(arg))))
 end
 
 
 function σ²ω⁰(arg::Vararg)
     return (:σ²ω⁰, collect(:($(arg))))
 end
 
 function σ²ω⁰(arg::Vector)
     return (:σ²ω⁰, arg)
 end
 # -------------------

  macro sigma2_w_star(arg::Vararg)
     return (:σ²w_star, collect(:($(arg))))
 end
 
 function sigma2_w_star(arg::Vararg)
     return (:σ²w_star, collect(:($(arg))))
 end
 
 function sigma2_w_star(arg::Vector)
     return (:σ²w_star, arg)
 end

  # -------------------

  macro σ²w_star(arg::Vararg)
     return (:σ²w_star, collect(:($(arg))))
 end
 
 function σ²w_star(arg::Vararg)
     return (:σ²w_star, collect(:($(arg))))
 end
 
 function σ²w_star(arg::Vector)
     return (:σ²w_star, arg)
 end

  # -------------------


macro timevar(arg::Vararg)
    return (:timevar, collect(:($(arg))))
end

macro idvar(arg::Vararg)
    return (:idvar, collect(:($(arg))))
end


macro Gvar(arg::Vararg)
    return (:Gvar, collect(:($(arg))))
end

 # -------------------
 function all_init(arg::Float64)
     return (:all_init, arg)
 end
 

  function all_init(arg::Int64)
     return (:all_init, Float64(arg))
 end
 
 
 function warmstart_solver(arg=nothing)
     return (:warmstart_solver, arg)
 end
 
 
  # -------------------

 function warmstart_maxIT(arg=nothing)
     return (:warmstart_maxIT, arg)
 end
 

 function main_solver(arg=nothing)
     return (:main_solver, arg)
 end
 
 
 function main_maxIT(arg::Any=nothing)
     return (:main_maxIT, arg)
 end
 
 
  # -------------------
 
 function tolerance(arg::Float64=nothing)
     return (:tolerance, arg)
 end
 
   # --------------------------
 
 function silent(arg::Bool=true)
     return (:silent, arg)
 end
 

 function useData(D::DataFrame)  # doesn't work using macro (perhaps because of data), so...
     return D
 end
 
 #? ---- macros for predict() -------------

 macro eq(arg)   
    return [:($(arg))]
end 