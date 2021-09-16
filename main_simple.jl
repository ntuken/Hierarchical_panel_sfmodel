# global variable list:
# _dicM: specify dicttionary, construct it each time using spec(), and update it once you call spec()
# _dicINI: initial vector for optimization, construct it each time using spec(), and change it in init_vec()
# _dicOPT: maximization options dictionary, construct it in opt(), not necessary to update it each time


########################################################
####                spec()                  ####
########################################################

"""
    spec(<keyword arguments>)
Provide specifications of the hierarchical panel data stochastic frontier model.    
# Arguments
- `@depvar(::symbol)`: the dependent variable from a DataFrame.
- `@frontier(::symbol)`: a list of variables in the frontier function (not include timevar, specify timevar independently).
- `@timevar(::symbol)`: time variable from DataFrame
- `@idvar(::symbol)`: id(firm) variable from DataFrame
- `@Gvar(::symbol)`: group variable from DataFrame
!!! not specify cᵢ, w_g, uᵢₜ since there is no exogenous determinant in this model, i.e. @cᵢ(_cons), @w_g(_cons), @uᵢₜ(_cons), which is trivial 
    to specify
"""
function spec(arg::Vararg) 

    printstyled("\n###------------------------------------###\n"; color=:yellow)
    printstyled("###  Estimating hierachy panel data stochastic frontier models using Julia  ###\n"; color=:yellow)
    printstyled("###------------------------------------###\n\n"; color=:yellow)

    global _dicM 
           _dicM = Dict{Symbol, Any}()  # nullify and initiate new dictionaries when a new model is specified

     #* Create a new ini vector for each new run; otherwise, result from the
     #*   previous run will be brought here. Does not do this for _dicOPT anf 
     #*   since we allow previous OPT settings to be used in the current run.
    global _dicINI  # nullify the vector for each new run; otherwise, result from the previous run will be brought here
           _dicINI = Dict{Symbol, Any}()



        #* -- creates default values ---

        # user need to specify:(:depvar, :frontier, :timevar, :idvar, :Gvar)
        for k in (:depvar, :frontier, :timevar, :idvar, :Gvar) 
            _dicM[k] = nothing
        end

        # user can't specify (:σ²ᵤ₀, :σ²ᵤ_star, :σ²₍₀, :σ²₍_star, :σ²w⁰, :σ²w_star), since there is no exogenous determinant 
        # reason for creating default value here just for consistency
        for i in (:σ²ᵤ₀, :σ²ᵤ_star, :σ²₍₀, :σ²₍_star, :σ²w⁰, :σ²w_star)
            _dicM[i] = :_cons
        end

        
        #* -- replace the defaults with the user's value ---          

        for d in :($(arg))
            _dicM[d[1]] = d[2]
        end 

        #* -- check input

        for k in (:depvar, :frontier, :timevar, :idvar, :Gvar)
            if _dicM[k] == nothing
                throw("For hierarchical panel data stochastic frontier model, $k is required but missing in spec() function. ")
            end
        end

        for k in (:σ²ᵤ₀, :σ²ᵤ_star, :σ²₍₀, :σ²₍_star, :σ²w⁰, :σ²w_star)
            if _dicM[k] != :_cons
                throw("In hierarchical panel data stochastic frontier model, there is no exogenous determinant for $k, hence $k would not be specified in spec() function.")
            end
        end




        return _dicM # for debugging purpose

end  # end of spec()


########################################################
####                init_vec()                  ####
########################################################

function init_vec(arg::Vararg) # create a dictionary of inital vectors
  # in addition to frontier and timevar, don't forget to create initial vector for sigma2_u_0,sigma2_u_star,....
  # two way to get initial vector
  # for spec(@frontier(x1,x2), @timevar(T), @Gvar(G), @idvar(id), @depvar(y)) :
  # e.g.1. init_vec(@frontier(0.5,0.5), @timevar(0.6), @sigma2_u_0(0.1), ....) 
  # e.g.2  init_vec(all_init(1))  

    global _dicINI
    _dicINI = Dict{Symbol, Any}()

    for d in :($(arg))
        _dicINI[d[1]] = d[2]
    end        

    
    # all possible keys: :frontier,:timevar, :σ²ᵤ₀, :σ²ᵤ_star, :σ²₍₀, :σ²₍_star, :σ²ω⁰, :σ²w_star
    
    

    return _dicINI # for debugging purpose
end    


########################################################
####                opt()                   ####
########################################################


function opt(arg::Vararg) # create a dictionary of maximization options

    global _dicOPT
        _dicOPT = Dict{Symbol, Any}()

    #* -- creates the default ---

    _dicOPT[:warmstart_solver] = :(NelderMead())
    _dicOPT[:warmstart_maxIT]  =  100
    _dicOPT[:main_solver]      = :(NewtonTrustRegion())
    _dicOPT[:main_maxIT]       =  2000
    _dicOPT[:tolerance]        =  1.0e-8
    _dicOPT[:silent]           =  false

    #* -- replace the defaults with the user's value ---

    for d in :($(arg))
        _dicOPT[d[1]] = d[2]
    end    
    
    #* ---- error checking --

    if (_dicOPT[:main_solver] == nothing) || (_dicOPT[:main_maxIT] == nothing) || (_dicOPT[:tolerance] == nothing)
         throw("Neither of the keys of @main_solver(), @main_maxIT(), or @tolerance() could be empty (`()`). If have doubts, you may delete these @ keys to use default values.")
    end

    return _dicOPT # for debugging purpose

end    # end of opt()




########################################################
###                 fit()                   ####
########################################################


function fit(sfdat::DataFrame) #, D1::Dict = _dicM, D2::Dict = _dicINI, D3::Dict = _dicOPT)

    #* ###### Check if the OPT dictionary exists #####
    #  --- No need for the INI dictionary because ---- #
    #  ---   a new dict is created in spec() - #
    #  ---   for each new run. ---------------------#

    @isdefined(_dicOPT) || opt() # if not exist, create one with default values

    #* ##### Get variables from dataset #######
    # pos: (begx, endx, begT, endT, pos_log_σ²ᵤ₀ ...); variables' positions in the parameter vector.
    # num: (nofobs, nofx, nofpara); number of variables in each equation
    # eqvec: ("frontier"=2, "timevar"=6,...); named tuple of equation names and equation position in the table
    # eqvec2: (coeff_frontier=(1,3), coef_Time=(4,5), coeff_log_σ²ᵤ₀ = 6...); named tuple of equation and parameter positions, for predict()
    # varlist: ("x1", "x2",...); variable names for making table

    (minfo1, minfo2, pos, num, eqvec, eqvec2, yvar, xvar, nofG, noffirm, noftime, varlist) = getvar(sfdat)

    
    #* ### print preliminary information ########

    if _dicOPT[:silent] == false

        printstyled("*********************************\n "; color=:cyan)
        printstyled("      Model Specification:\n"; color=:cyan); 
        printstyled("*********************************\n"; color=:cyan)

        print("Model type: "); printstyled(minfo1; color=:yellow); println();println()
        printstyled(minfo2; color=:yellow); println()
    end

    #* ########## Process initial value dictionary  #####
    #* --- Get OLS results and other auxiliary values. --- #

    β0     = xvar \ yvar;         # OLS estiamte, uses a pivoted QR factorization;

    resid  = yvar - xvar*β0
    sse    = sum((resid).^2)  
    ssd    = sqrt(sse/(size(resid)[1])) # sample standard deviation; σ² = (1/N)* Σ ϵ^2
    ll_ols = sum(normlogpdf.(0, ssd, resid)) # ols log-likelihood
    sk_ols = sum((resid).^3) / ((ssd^3)*(size(resid)[1])) # skewnewss of ols residuals

    #* --- Create the dictionary -----------

    if (:all_init in keys(_dicINI))
        init_vec = _dicINI[:all_init] * ones(num.nofpara)
    else
        #*  Create ini vectors from user's values; if none, use the default.--- #      
        b_ini  = get(_dicINI, :frontier, β0)
        # :σ²ᵤ₀, :σ²ᵤ_star, :σ²₍₀, :σ²₍_star, :σ²w⁰, :σ²w_star

        log_σ²ᵤ₀_ini     = get(_dicINI, :σ²ᵤ₀, 0.1)
        log_σ²ᵤ_star_ini = get(_dicINI, :σ²ᵤ_star, 0.1)
        log_σ²₍₀_ini     = get(_dicINI, :σ²₍₀, 0.1)
        log_σ²₍_star_ini = get(_dicINI, :σ²₍_star, 0.1)
        log_σ²w⁰_ini     = get(_dicINI, :σ²w⁰, 0.1)
        log_σ²w_star_ini = get(_dicINI, :σ²w_star, 0.1)


        init_vec = vcat(b_ini, log_σ²ᵤ₀_ini, log_σ²ᵤ_star_ini, log_σ²₍₀_ini, log_σ²₍_star_ini, log_σ²w⁰_ini, log_σ²w_star_ini)  
        init_vec = vec(init_vec)    
    end # if :all_init


    #* ############ Misc.  ################     
    # --- check if the number of initial values is correct 
    (length(init_vec) == num.nofpara) || throw("The number of initial values does not match the number of parameters to be estimated. Make sure the number of init values in init() matches the number of variabls in spec(). Number of parameters required is: $(num.nofx) (frontier) + 6, while length of init_vector is $(length(init_vec)) ($(init_vec))")

    # --- Make sure there is no numerical issue arising from int vs. Float64.
    init_vec = convert(Array{Float64,1}, init_vec) 

    #* ############# process optimization dictionary  #######

    if (_dicOPT[:warmstart_solver] == nothing) || (_dicOPT[:warmstart_maxIT] == nothing)
        do_warmstart_search = 0
    else 
        do_warmstart_search = 1
        sf_ini_algo  = eval(_dicOPT[:warmstart_solver])  # warmstart search algorithms
        sf_ini_maxit = _dicOPT[:warmstart_maxIT]         # warmstart search iter limit
    end    

    # ---- main maximization algorithm -----
    sf_algo  = eval(_dicOPT[:main_solver])    # main algorithm
    sf_maxit = _dicOPT[:main_maxIT] 
    sf_tol   = _dicOPT[:tolerance] 

    #* ########  Start the Estimation  ##########

    #* ----- Define the problem's Hessian -----#


    _Hessian = TwiceDifferentiable(rho -> LL_T( 
                   yvar, xvar, pos, rho,
                   nofG, noffirm, noftime),
             init_vec;               
             autodiff=:forward); 

    #* ---- Make placeholders for dictionary recording purposes *#

    sf_init_1st_dic  = 0
    sf_init_2nd_dic  = 0
    sf_ini_algo_dic  = nothing
    sf_ini_maxit_dic = 0
    sf_total_iter    = 0

    _run = 1  # a counter; use the -if- instead of -for- to avoid using global variables

    if (do_warmstart_search == 1) && (_run == 1)  

        if _dicOPT[:silent] == false
            printstyled("The warmstart run...\n\n"; color = :green)
        end

        sf_init_1st_dic  = copy(init_vec) # for dict recording
        sf_ini_algo_dic  = sf_ini_algo
        sf_ini_maxit_dic = copy(sf_ini_maxit)

        # 
        @time mfun = optimize(_Hessian, 
                           init_vec,         # initial values  
                           sf_ini_algo,                   
                           Optim.Options(g_tol = sf_tol,
                                         iterations  = sf_ini_maxit, 
                                         store_trace = true,
                                         show_trace  = false))


        sf_total_iter += Optim.iterations(mfun) # for later use

        init_vec = Optim.minimizer(mfun)  # save as initials for the next run and update init_vec by nongradient optimize result
        _run    = 2                      # modify the flag

        if _dicOPT[:silent] == false
            println()
            # print("The warmstart run is:\n$mfun \n")  
            print("$mfun \n")
            print("The warmstart results are:\n"); printstyled(Optim.minimizer(mfun); color=:yellow); println("\n")
        end

    end  # if  (do_warmstart_search == 1) && (_run == 1)  

    if (do_warmstart_search == 0 ) || (_run == 2) # either no warmstart run and go straight here, or the 2nd run

        sf_init_2nd_dic = copy(init_vec) # for dict recording 

        if _dicOPT[:silent] == false
            println()
            printstyled("Starting the optimization run...\n\n" ; color = :green)
        end 

        @time mfun = optimize(_Hessian, 
                          init_vec,       # initial values  
                          sf_algo,       # different from search run
                          Optim.Options(g_tol = 1e-8, #sf_tol,
                                        iterations  = sf_maxit, # different from search run
                                        store_trace = true,
                                        show_trace  = false))
        sf_total_iter += Optim.iterations(mfun)

        if _dicOPT[:silent] == false
            println()
            # print("The optimization problem is:\n$mfun \n")  
            print("$mfun \n")  
            print("The resulting coefficient vector is:\n"); printstyled(Optim.minimizer(mfun); color=:yellow); println("\n")
        end 
        if Optim.iteration_limit_reached(mfun) == true
         printstyled("Caution: The number of iterations reached the limit.\n\n"; color= :red)  
        end  

    end     # if (do_warmstart_search == 0 )....


    #* ###### Post-estimation process ############### 

    _coevec            = Optim.minimizer(mfun)  # coef. vec.
    numerical_hessian  = hessian!(_Hessian, _coevec)  # Hessain

    #* ------ Check if the matrix is invertible. ----
    # default main algorithm is NewtonTrustRegion, which is Hessian required
    var_cov_matrix = try
                 inv(numerical_hessian)
              catch err 
                 throw("The Hessian matrix is not invertible, indicating the model does not converge properly. The estimation is abort.")
              end 
    #* In some cases the matrix is invertible but the resulting diagonal
    #*    elements are negative. Check.
    all( diag(var_cov_matrix) .> 0 ) || throw("Some of the diagonal elements of the var-cov matrix are non-positive, indicating problems in the convergence. The estimation is abort.")


    #* ------- Make Table ------------------

    stddev  = sqrt.(diag(var_cov_matrix)) # standard error
    t_stats = _coevec ./ stddev          # t statistics
    p_value = ones(num.nofpara)   # p values
    ci_low  = ones(num.nofpara) # confidence interval
    ci_upp  = ones(num.nofpara) 
    tt      = cquantile(Normal(0,1), 0.025)

    for i = 1:num.nofpara 
        p_value[i,1] = pvalue(TDist(num.nofobs - num.nofpara), t_stats[i,1]; tail=:both)
        ci_low[i,1] = _coevec[i,1] - tt*stddev[i,1]
        ci_upp[i,1] = _coevec[i,1] + tt*stddev[i,1]
    end  

    #* Build the table columns *#

    table = ones(num.nofpara, 7)  # 7 columns in the table
    table[:,2] = _coevec   # estiamted coefficients
    table[:,3] = stddev    # std deviation
    table[:,4] = t_stats   # t statistic
    table[:,5] = p_value   # p value
    table[:,6] = ci_low
    table[:,7] = ci_upp
    table      = [" " "Coef." "Std. Err." "z" "P>|z|" "95%CI_l" "95%CI_u"; table]  # add to top of the table

    #*  creating a column of function names 

    table[:, 1] .= ""
    for i in (1:length(eqvec))
        j = eqvec[i]
        table[j,1] = keys(eqvec)[i]
    end

    #*  Add the column of variable names

    table = hcat(varlist, table)                      # combine the variable names column (hcat, horizontal concatenate; see also vcat)
    table[:,1], table[:,2] = table[:,2], table[:,1]   # swap the first name column and the function column
    table[1,2] = "Var."

    # * ------ Print Results ----------- *#

    if _dicOPT[:silent] == false

        printstyled("*********************************\n "; color=:cyan)
        printstyled("      Estimation Results:\n"; color=:cyan); 
        printstyled("*********************************\n"; color=:cyan)

        print("Model type: "); printstyled(minfo1; color=:yellow); println()
        print("Number of observations: "); printstyled(num.nofobs; color=:yellow); println()
        print("Number of total iterations: "); printstyled(sf_total_iter; color=:yellow); println()
        if Optim.converged(mfun) == true
            print("Converged successfully: "); printstyled(Optim.converged(mfun); color=:yellow); println()
        elseif Optim.converged(mfun) == false
            print("Converged successfully: "); printstyled(Optim.converged(mfun); color=:red); println()
        end         
        print("Log-likelihood value: "); printstyled(round(-1*Optim.minimum(mfun); digits=5); color=:yellow); println()
        println()

        pretty_table(table[2:size(table)[1],:], 
                  ["", "Var.", "Coef.", "Std.Err.", "z", "P>|z|", 
                   "95%CI_l", "95%CI_u"],
                  formatters = ft_printf("%5.4f", 3:8),
                  compact_printing = true)
        println()


        # *----- Auxiliary Table, log parameters to original scales --------



        auxtable = Array{Any}(undef,6,3)
        rn = 1 # row index

        pos_key = [:pos_log_σ²ᵤ₀, :pos_log_σ²ᵤ_star, :pos_log_σ²₍₀, :pos_log_σ²₍_star, :pos_log_σ²w⁰, :pos_log_σ²w_star]
        parameters_vec = [:σ²ᵤ₀, :σ²ᵤ_star, :σ²₍₀, :σ²₍_star, :σ²w⁰, :σ²w_star]


        for key in pos_key
            auxtable[rn, 1] = parameters_vec[rn]
            auxtable[rn, 2] = exp(_coevec[pos[key]])
            auxtable[rn, 3] = exp(_coevec[pos[key]])*stddev[pos[key]]   
            rn += 1
        end


        println("Convert the constant log-parameter to its original scale, e.g., σ²ᵤ₀ = exp(log_σ²ᵤ₀):")   
        pretty_table(auxtable[1:6,:],
                  ["", "Coef.", "Std.Err."],
                  formatters = ft_printf("%5.4f", 2:3),
                  compact_printing = true)
        println()


        printstyled("***** Additional Information *********\n"; color=:cyan)
        print("OLS (frontier only) log-likelihood: "); printstyled(round(ll_ols; digits=5); color=:yellow); println("")
        print("Skewness of OLS residuals: "); printstyled(round(sk_ols; digits=5); color=:yellow); println("\n")
        println("Use `name.list` to see saved results (keys and values) where `name` is the return specified as in `name = fit(..)`. Values may be retrieved using the keys; e.g., `name.loglikelihood` returns the log-likelihood value of the model, and `name.frontier` returns a vector of frontier names specified in spec(). Use `keys(name)` to see available keys.")

        printstyled("**************************************\n\n\n"; color=:cyan)

    end  # if_silent

    #* ----- Test the convergence. Will show error if the test fails. ----

    isnan(Optim.g_residual(mfun)) || (Optim.g_residual(mfun) < 1e-4) || 
    (printstyled("Note that the estimation may not have converged properly. Either the gradient is too large or the iteration hits the limit.\n\n", color = :red))

    #* ########### create a dictionary and make a tuple for return ########### *#

    _dicRES = OrderedDict{Symbol, Any}()     
    _dicRES[:converged]          = Optim.converged(mfun)
    _dicRES[:iter_limit_reached] = Optim.iteration_limit_reached(mfun)
    _dicRES[:_______________] = "___________________"  #33
    _dicRES[:n_observations]  = num.nofobs
    _dicRES[:loglikelihood]   = -Optim.minimum(mfun)
    _dicRES[:table]           = [table][1]
    _dicRES[:coeff]           = _coevec
    _dicRES[:std_err]         = stddev
    _dicRES[:var_cov_mat]     = [var_cov_matrix][1]
    _dicRES[:OLS_loglikelihood] = ll_ols
    _dicRES[:OLS_resid_skew]    = sk_ols
    _dicRES[:_____________] = "___________________"  #31      
    _dicRES[:model]         = minfo1      
    _dicRES[:depvar]        = _dicM[:depvar]
    _dicRES[:frontier]      = _dicM[:frontier]
    _dicRES[:timevar]       = _dicM[:timevar]
    _dicRES[:idvar]         = _dicM[:idvar]
    _dicRES[:Gvar]          = _dicM[:Gvar]

    # _dicM of σ²ᵤ₀, σ²ᵤ_star, σ²₍₀, σ²₍_star, σ²w⁰, σ²w_star are actually _cons, store it into _dicRES just for consistency
    _dicRES[:σ²ᵤ₀]                  = _dicM[:σ²ᵤ₀]
    _dicRES[:log_σ²ᵤ₀]              = _dicM[:σ²ᵤ₀]
    _dicRES[:σ²ᵤ_star]              = _dicM[:σ²ᵤ_star]
    _dicRES[:log_σ²ᵤ_star]          = _dicM[:σ²ᵤ_star]
    _dicRES[:σ²₍₀]                  = _dicM[:σ²₍₀]
    _dicRES[:log_σ²₍₀]              = _dicM[:σ²₍₀]
    _dicRES[:σ²₍_star]              = _dicM[:σ²₍_star]
    _dicRES[:log_σ²₍_star]          = _dicM[:σ²₍_star]
    _dicRES[:σ²w⁰]                  = _dicM[:σ²w⁰]
    _dicRES[:log_σ²w⁰]              = _dicM[:σ²w⁰]
    _dicRES[:σ²w_star]              = _dicM[:σ²w_star]
    _dicRES[:log_σ²w_star]          = _dicM[:σ²w_star]




    for i in 1:length(eqvec2)
        _dicRES[keys(eqvec2)[i]] = _coevec[eqvec2[i]]
    end

    _dicRES[:________________]  = "___________________" #34
    _dicRES[:Hessian]           = [numerical_hessian][1]
    _dicRES[:gradient_norm]     = Optim.g_residual(mfun)
    # _dicRES[:trace]             = Optim.trace(mfun)     # comment out because not very informative and size could be large
    _dicRES[:actual_iterations] = Optim.iterations(mfun)
    _dicRES[:______________] = "______________________" #32
    _dicRES[:warmstart_solver] = sf_ini_algo_dic
    _dicRES[:warmstart_ini]    = sf_init_1st_dic
    _dicRES[:warmstart_maxIT]  = sf_ini_maxit_dic
    _dicRES[:main_solver]      = sf_algo
    _dicRES[:main_ini]         = sf_init_2nd_dic
    _dicRES[:main_maxIT]       = sf_maxit
    _dicRES[:tolerance]        = sf_tol
    _dicRES[:eqpo]             = eqvec2


    #* ----- Create a NamedTuple from the dic as the final output; 
    #* -----     put the dic in the tuple.

    _ntRES = NamedTuple{Tuple(keys(_dicRES))}(values(_dicRES))
    _ntRES = (; _ntRES..., list    = _dicRES)

    #* ---- Create a gloal dictionary for predict() ---- 

    global _eqncoe 
    _eqncoe = Dict{Symbol, Any}()  # nullify and initiate new dictionaries when a new model is specified
    
    for i in 1:length(eqvec2)
        _eqncoe[keys(eqvec)[i]]  = _coevec[eqvec2[i]]  # for predict()
    end

#* ############  make returns  ############ *#

    return _ntRES

end # fit()


########################################################
####            catching type error                 ####
########################################################


function init_vec(arg::Vector) 
    throw(ArgumentError("The initial values in init_vec() are specified incorrectly. They must be supplied in macros such as @frontier(0.1, 0.1)."))
end    

function init_vec(args::Number...) 
    throw(ArgumentError("The initial values in init_vec() are specified incorrectly. They must be supplied in macros such as @frontier(0.1, 0.1)."))
end    

function fit(sfdat::Any) 
    throw(ArgumentError("The dataset specified in fit() must be a DataFrame."))
end