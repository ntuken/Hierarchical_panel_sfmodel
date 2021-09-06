
function GID_Dict(Gvar, ivar)
  # output G, ID paired dictionary
    GID_dict = Dict()  # intitialize dictionary whose key is Group value(Int64) and value is all ID in given group
    GID_pair = unique(zip(Gvar, ivar))

    for (g, ID) in GID_pair
        if haskey(GID_dict, g)
            GID_dict[g] = push!(GID_dict[g], ID)
        else
            GID_dict[g] = [ID]
        end
    end
    return GID_dict
end


function IDRow_Dict(ivar)
  # output ID, Row index paired dictionary
    IDRow_dict = Dict()  # intitialize dictionary whose key is Group value(Float64) and value is all ID in given group
    unique_ID = unique(ivar)               

    for id in unique_ID
        
        IDRow_dict[id] = findall(ivar .== id)
    end

    return IDRow_dict

end


function getvar(dat::DataFrame)

    yvar = dat[:, _dicM[:depvar]]   # still a DataFrame
    xvar = dat[:, _dicM[:frontier]]

    tvar = dat[:, _dicM[:timevar][1]] # the [1] make _dicM[:timevar] element, not vector ==> tvar is vector not dataframe
    Gvar = dat[:, _dicM[:Gvar]]
    ivar = dat[:, _dicM[:idvar][1]] # need the [1] otherwise causes problems down the road

    #* --- model info printout --------- 

    modelinfo1 = "A hierachy panel data stochastic frontier model for the estimation of stochastic metafrontiers"

    modelinfo2 = begin
    """
    $(_dicM[:depvar][1]) = frontier(dₜ + $(_dicM[:frontier])) + cᵢ + w_g(i) + uᵢₜ,
     
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

    #* --- get dummy variable of timevar --- *#

    #              !!!!!!!!!!!!!   waiting for coding   !!!!!!!!!!!!!
    #                              output: Tvar matrix          #
    
    unique_T = unique(tvar)  

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

    #* --- retrieve and generate important parameters -----

    #*   number of obs and number of variables
    nofx  = nofT = 0  # to make a complete list

    nofobs  = nrow(dat)    
    nofx = size(xvar)[2]  # nofx: number of x vars
    nofT = size(Tvar)[2]  # nofT: number of unique time in data set


  
    nofpara = nofx + nofT + 6 # coefficient β, D and log_σ²₍₀, log_σ²₍*, log_σ²ᵤ₀, log_σ²ᵤ*, log_σ²w⁰, log_σ²w*


    nofvar = (nofobs=nofobs, nofx=nofx, nofT=nofT, nofpara=nofpara)

    #* positions of the variables/parameters
    begx=endx=begT=endT=pos_log_σ²ᵤ₀=pos_log_σ²ᵤ_star=pos_log_σ²₍₀=pos_log_σ²₍_star=pos_log_σ²w⁰=pos_log_σ²w_star=0 #intitialize

    begx = 1
    endx = nofx
    # begz = endx + 1
  # endz = begz + nofz -1
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
    xnames  = names(xvar)
  # znames  = names(zvar)
    Tnames  = names(Tvar)
    
    varlist = vcat(" ", xnames, Tnames, 
                        ["log_σ²ᵤ₀","log_σ²ᵤ*", "log_σ²₍₀", "log_σ²₍*", "log_σ²w⁰", "log_σ²w*" ])

   
    #* Converting the dataframe to matrix in order to do computation
    yvar  = convert(Array{Float64}, yvar)
    xvar = convert(Array{Float64}, xvar)
    Tvar = convert(Array{Float64}, Tvar)
    Gvar  = convert(Array{Float64}, Gvar)
    ivar  = convert(Array{Float64}, ivar)


    
    #* panel info and within transformation


    IDRow_dict = IDRow_Dict(ivar)   # (Nx2): col_1 is panel's row info; col_2 is panel's number of periods
    GID_dict = GID_Dict(Gvar, ivar)       # (Nx2): col_1 is group's row info; col_2 is ∑ (id_g * time_g) for each group

    

    return modelinfo1, modelinfo2, posvec, nofvar, eqvec, eqvec2, yvar, xvar, Tvar, IDRow_dict, GID_dict, varlist
end  # getvar



#=function get_rowG(Gvar)
    # output: vector of matrix
   

    nofG  = length(unique(Gvar)) # nofG = number of unique group
    G = Array{Int64}(undef,nofG,2) # first column store unique group label, second store #row each group
    G[:,1] = unique(Gvar)    # list of group label without repeatation
    for i = 1:nofG
        G[i,2]= sum(Gvar.== G[i,1])    # number of rows each group

    end
    T_g= G[:,2] 

    rowG = Vector{Vector{Int64}}()    # store row index of same group
    for i=1:nofG
        cc = findall(x-> x == G[i,1], Gvar) # row index of i'th group
        push!(rowG, cc) # put the group info in the vector
    end    
  
    rowGidt = hcat(rowG, T_g) 

    return rowGidt # (Nx2): col_1 is group's row info; col_2 is ∑ (id_g * time_g) for each group
end=#

#=function get_rowIDT(ivar)
    # output: vector and vector of matrix
    

    N  = length(unique(ivar)) # N=number of panels
    id = Array{Int64}(undef,N,2)
    id[:,1] = unique(ivar)    # list of id with no repetition
    for i = 1:N
        id[i,2]= sum(ivar.== id[i,1])    # number of periods for the panel
    end
    Tᵢ = id[:,2]

    rowID = Vector{Vector{Int64}}()    # store individual's row number
    for i=1:N
        cc = findall(x-> x == id[i,1], ivar) # row index of i'th firm
        push!(rowID, cc) # put the id info in the vector
    end    
  
    rowIDT = hcat(rowID, Tᵢ) 

    return rowIDT # (Nx2): col_1 is panel's row info; col_2 is panel's number of periods
end=#
