# Assume noftime and noffirm is constant, i.e. nofobsinG(noftime*noffirm) is constant
# Assume no time fixed effect in this model

function getvar(dat::DataFrame)
    # intercept_exist parameter is true if dat have intercept column (i.e. _cons column)
    
    #* --- model info printout --------- 

    modelinfo1 = "A hierachy panel data stochastic frontier model for the estimation of stochastic metafrontiers"

    modelinfo2 = begin
    """
    $(_dicM[:depvar][1]) = frontier($(_dicM[:frontier])) + cᵢ + w_g(i) + uᵢₜ,
     
    where uᵢₜ = u⁰ᵢₜ - u*ᵢₜ
          uᵒᵢₜ ∼ N(0, σ²ᵤ₀)
            σ²ᵤ₀ = exp(log_σ²ᵤ₀)
                 = exp($(_dicM[:σ²ᵤ₀]));
          u*ᵢₜ ∼ N⁺(0, σ²ᵤ*)
          σ²ᵤ₀ = exp(log_σ²ᵤ*)
               = exp($(_dicM[:σ²ᵤ_star]));

          cᵢ = c⁰ᵢ - cᵢ*
          c⁰ᵢ ∼ N(0, σ²₍₀)
          σ²₍₀ = exp(log_σ²₍₀) 
               = exp($(_dicM[:σ²₍₀]));
          cᵢ* ∼ N⁺(0, σ²₍*)
          σ²₍* = exp(log_σ²₍*) 
               = exp($(_dicM[:σ²₍_star]));

          w_g = w⁰_g - w*_g
          wᵒ_g ∼ N(0, σ²w⁰)
          σ²w⁰ = exp(log_σ²w⁰)
               = exp($(_dicM[:σ²w⁰]));
          w*_g ∼ N⁺(0, σ²w*)
          σ²w* = exp(log_σ²w*)
               = exp($(_dicM[:σ²w_star]));
     """
    end
    yvar = dat[:, _dicM[:depvar]]   # still a DataFrame
    xvar = dat[:, _dicM[:frontier]]
    Gvar = dat[:, _dicM[:Gvar][1]]
    ivar = dat[:, _dicM[:idvar][1]] # need the [1] to make _dicM[:idvar][1] element => ivar is vector not DataFrame
    tvar = dat[:, _dicM[:timevar][1]]
    
    data = deepcopy(dat)  # using deepcopy to make data dataframe independent to original dat (i.e. changing data don't change dat )

    # sort data by group, id, time
    sort!(data,[_dicM[:Gvar][1], _dicM[:idvar][1], _dicM[:timevar][1]])  # the [1] make _dicM[:timevar] element, not vector

    first_group = data[1, _dicM[:Gvar][1]]
    nofobs = sum(data[:, _dicM[:Gvar][1]] .== first_group)  # number of row in first group


    unique_T = unique(tvar)
    noftime = length(unique_T)
    noffirm = Int(nofobs/noftime)
    nofG = length(unique(Gvar))


    #* --- retrieve and generate important parameters -----

    #*   number of obs and number of variables

    nofobs  = nrow(data)  
    nofx = length(_dicM[:frontier])  # nofx: number of x vars

    nofpara = nofx + 6 # coefficient β and log_σ²₍₀, log_σ²₍*, log_σ²ᵤ₀, log_σ²ᵤ*, log_σ²w⁰, log_σ²w*


    nofvar = (nofobs=nofobs, nofx=nofx, nofpara=nofpara)

    #* positions of the variables/parameters
    begx=endx=pos_log_σ²ᵤ₀=pos_log_σ²ᵤ_star=pos_log_σ²₍₀=pos_log_σ²₍_star=pos_log_σ²w⁰=pos_log_σ²w_star=0 #intitialize

    begx = 1
    endx = nofx
    
    pos_log_σ²ᵤ₀ = endx + 1
    pos_log_σ²ᵤ_star = endx + 2
    pos_log_σ²₍₀ = endx + 3
    pos_log_σ²₍_star = endx + 4
    pos_log_σ²w⁰ = endx + 5
    pos_log_σ²w_star = nofpara


  

    posvec = (begx=begx, endx=endx,
              pos_log_σ²ᵤ₀=pos_log_σ²ᵤ₀, pos_log_σ²ᵤ_star=pos_log_σ²ᵤ_star, pos_log_σ²₍₀=pos_log_σ²₍₀, pos_log_σ²₍_star=pos_log_σ²₍_star,
              pos_log_σ²w⁰=pos_log_σ²w⁰, pos_log_σ²w_star=pos_log_σ²w_star)

    #* create equation names and mark positions for making tables
    eqvec = (frontier = begx + 1, 
            log_σ²ᵤ₀  = pos_log_σ²ᵤ₀ + 1,
         log_σ²ᵤ_star = pos_log_σ²ᵤ_star + 1,
             log_σ²₍₀ = pos_log_σ²₍₀ + 1,
         log_σ²₍_star = pos_log_σ²₍_star + 1,
             log_σ²w⁰ = pos_log_σ²w⁰ + 1,
         log_σ²w_star = pos_log_σ²w_star + 1 )

    #* create equation names and mark positions 
    eqvec2 = (coeff_frontier      = (begx:endx), 
              coeff_log_σ²ᵤ₀      = pos_log_σ²ᵤ₀,
              coeff_log_σ²ᵤ_star  = pos_log_σ²ᵤ_star,
              coeff_log_σ²₍₀      = pos_log_σ²₍₀,
              coeff_log_σ²₍_star  = pos_log_σ²₍_star,
              coeff_log_σ²w⁰      = pos_log_σ²w⁰,
              coeff_log_σ²w_star  = pos_log_σ²w_star)

    #* retrieve variable names for making tables
    xnames  = names(data[!, _dicM[:frontier]])
    
    varlist = vcat(" ", xnames, 
                        ["log_σ²ᵤ₀","log_σ²ᵤ*", "log_σ²₍₀", "log_σ²₍*", "log_σ²w⁰", "log_σ²w*" ])

   
    
    #* Converting the dataframe to matrix in order to do computation
    yvar  = convert(Array{Float64}, data[:, _dicM[:depvar][1]])  # vector{Float64}
    xvar = convert(Array{Float64}, data[:, _dicM[:frontier]])

    return modelinfo1, modelinfo2, posvec, nofvar, eqvec, eqvec2, yvar, xvar, nofG, noffirm, noftime , varlist
end  # getvar