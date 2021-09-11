function getvar(dat::DataFrame)
    # intercept_exist parameter is true if dat have intercept column (i.e. _cons column)
    
    #* --- model info printout --------- 

    modelinfo1 = "A hierachy panel data stochastic frontier model for the estimation of stochastic metafrontiers"

    modelinfo2 = begin
    """
    $(_dicM[:depvar][1]) = frontier($(_dicM[:frontier])) + dₜ + cᵢ + w_g(i) + uᵢₜ,
     
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
    # yvar = dat[:, _dicM[:depvar]]   # still a DataFrame
    # xvar = dat[:, _dicM[:frontier]]
    # Gvar = dat[:, _dicM[:Gvar]]
    # ivar = dat[:, _dicM[:idvar][1]] # need the [1] otherwise causes problems down the road

    tvar = dat[:, _dicM[:timevar][1]] # the [1] make _dicM[:timevar] element, not vector ==> tvar is vector not dataframe
    
    data = deepcopy(dat)  # using deepcopy to make data dataframe independent to original dat (i.e. changing data don't change dat )

    # sort data by group, id, time
    data = sort!(data,[_dicM[:Gvar][1], _dicM[:idvar][1], _dicM[:timevar][1]])  # the [1] make _dicM[:timevar] element, not vector

    # create group_span and id_span to collect number of row in each group and id
    # e.g. 
    # data =   group     id    timevar
    #            1       1        1
    #            1       1        2
    #            1       2        1
    #            1       2        2
    #            2       3        1
    #            2       3        2
    # then group_span = {1=>4, 2=>2}, 4 is number of row in group 1
    # id_span = {1=>2, 2=>2, 3=>2}      

    # nofobs_in_sorted_G, noffirm, noftime
    nofobsinG_list = combine(groupby(data, _dicM[:Gvar][1]), _dicM[:timevar][1] => length => :nofobs_in_sortG).nofobs_in_sortG  # vector{Int64}
    noffirm_list = combine(groupby(data, _dicM[:Gvar][1]), _dicM[:idvar][1] => length ∘ unique => :n_distinct_firm_id).n_distinct_firm_id
    noftime_list = combine(groupby(data, [_dicM[:Gvar][1], _dicM[:idvar][1]]), _dicM[:timevar][1] => length => :n_distinct_firm_id).n_distinct_firm_id

    # group_span = sort!(OrderedDict(countmap(data[:,_dicM[:Gvar][1]])))
    # id_span = sort!(OrderedDict(countmap(data[:,_dicM[:idvar][1]])))

    # group_nofid = combine(groupby(data, :group_id), :firm_id => length ∘ unique => :n_distinct_firm_id)
    # group_nofid_dict = OrderedDict(zip(group_id[!, _dicM[:Gvar][1]], group_nofid.n_distinct_firm_id))

    # group_id =  combine(groupby(data, _dicM[:Gvar][1]), :firm_id => unique => :distinct_firm_id)
    # group_id_dict = OrderedDict(zip(group_id[!, _dicM[:Gvar][1]], group_id.distinct_firm_id))
    

    #* --- get dummy variable of timevar --- *#

    
    unique_T = sort(unique(tvar))  # sort to make sure time dummy variables (Tvar) are in ascending order

    Tvar = DataFrame()  # intitialize dummy time dataframe, not matrix since need to reserve time information in column label
    # e.g.  
    #       t = [            Tvar = t: 1   2   3
    #            1
    #            2                   [ 1   0   0
    #            3          ==>        0   1   0   
    #            1                     0   0   1  
    #             ]                    1   0   0 ]
    for Time in unique_T
        Tvar[!,"$Time"] = Float64.(tvar .== Time)
    end

    intercept_exist = Bool(_dicM[:is_intercept_exist][1])
    intercept_exist && (Tvar = Tvar[!,1:(end-1)])  # avoid dummy variable trap
    
    

    data = hcat(data[!,Not(_dicM[:timevar][1])], Tvar)  


    

    #* --- retrieve and generate important parameters -----

    #*   number of obs and number of variables
    nofx  = nofT = 0  # to make a complete list

    nofobs  = nrow(data)    
    nofx = length(_dicM[:frontier])  # nofx: number of x vars
    
    #* if intercept exist, then drop one of time period to avoid dummy trap, therefore nofT minus one
    intercept_exist ? (nofT = length(unique_T) - 1) :  (nofT = length(unique_T))

             
  
    nofpara = nofx + nofT + 6 # coefficient β, D and log_σ²₍₀, log_σ²₍*, log_σ²ᵤ₀, log_σ²ᵤ*, log_σ²w⁰, log_σ²w*


    nofvar = (nofobs=nofobs, nofx=nofx, nofT=nofT, nofpara=nofpara)

    #* positions of the variables/parameters
    begx=endx=begT=endT=pos_log_σ²ᵤ₀=pos_log_σ²ᵤ_star=pos_log_σ²₍₀=pos_log_σ²₍_star=pos_log_σ²w⁰=pos_log_σ²w_star=0 #intitialize

    begx = 1
    endx = nofx
    begT = endx + 1
    endT = begT + nofT -1
    

    pos_log_σ²ᵤ₀ = endT + 1
    pos_log_σ²ᵤ_star = endT + 2
    pos_log_σ²₍₀ = endT + 3
    pos_log_σ²₍_star = endT + 4
    pos_log_σ²w⁰ = endT + 5
    pos_log_σ²w_star = nofpara


  

    posvec = (begx=begx, endx=endx, begT=begT, endT=endT,
              pos_log_σ²ᵤ₀=pos_log_σ²ᵤ₀, pos_log_σ²ᵤ_star=pos_log_σ²ᵤ_star, pos_log_σ²₍₀=pos_log_σ²₍₀, pos_log_σ²₍_star=pos_log_σ²₍_star,
              pos_log_σ²w⁰=pos_log_σ²w⁰, pos_log_σ²w_star=pos_log_σ²w_star)

    #* create equation names and mark positions for making tables
    eqvec = (frontier = begx + 1, 
                timevar  = begT + 1,
            log_σ²ᵤ₀  = pos_log_σ²ᵤ₀ + 1,
         log_σ²ᵤ_star = pos_log_σ²ᵤ_star + 1,
             log_σ²₍₀ = pos_log_σ²₍₀ + 1,
         log_σ²₍_star = pos_log_σ²₍_star + 1,
             log_σ²w⁰ = pos_log_σ²w⁰ + 1,
         log_σ²w_star = pos_log_σ²w_star + 1 )

    #* create equation names and mark positions 
    eqvec2 = (coeff_frontier      = (begx:endx), 
              coeff_Time          = (begT:endT),
              coeff_log_σ²ᵤ₀      = pos_log_σ²ᵤ₀,
              coeff_log_σ²ᵤ_star  = pos_log_σ²ᵤ_star,
              coeff_log_σ²₍₀      = pos_log_σ²₍₀,
              coeff_log_σ²₍_star  = pos_log_σ²₍_star,
              coeff_log_σ²w⁰      = pos_log_σ²w⁰,
              coeff_log_σ²w_star  = pos_log_σ²w_star)

    #* retrieve variable names for making tables
    xnames  = names(data[!, _dicM[:frontier]])
    Tnames  = names(Tvar)
    
    varlist = vcat(" ", xnames, Tnames, 
                        ["log_σ²ᵤ₀","log_σ²ᵤ*", "log_σ²₍₀", "log_σ²₍*", "log_σ²w⁰", "log_σ²w*" ])

   
    
    #* Converting the dataframe to matrix in order to do computation
    yvar  = convert(Array{Float64}, data[:, _dicM[:depvar][1]])  # vector{Float64}
    xvar = convert(Array{Float64}, data[:, _dicM[:frontier]])
    Tvar = convert(Array{Float64}, Tvar)
    Gvar  = convert(Array{Float64}, data[:, _dicM[:Gvar][1]])    # vector{Float64}
    ivar  = convert(Array{Float64}, data[:, _dicM[:idvar][1]])   # vector{Float64}

    return modelinfo1, modelinfo2, posvec, nofvar, eqvec, eqvec2, yvar, xvar, Tvar, nofobsinG_list , noffirm_list, noftime_list, varlist
end  # getvar