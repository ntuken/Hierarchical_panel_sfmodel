
########################################################
####             equation prediction                ####
########################################################

#= """
    predict(@eq(eq_name), data::DataFrame)

Return predicted values of an equation (`eq_name`) of a stochastic frontier
model after the latter is estimated on `data`.

 
The predicted value is computed based on the equation's variable list and the
estimated coefficient vector. For instance, if the `frontier` function is a
linear function of variables ``X`` and coefficient vector β,
`predict(@eq(frontier), df)` returns `β̂ₓ`. In this model, variance of errors such as σ²ᵤ₀, σ²ᵤ_star, ...are parameterized by constant
. Hence, `predict(@eq(σ²₍₀), df)` will return estimation of σ²₍₀, which is equal to exp(log_σ²₍₀)
`predict(@eq(log_σ²₍₀), df)` will return log_σ²₍₀ estimation.

# _eqncoe is global dictionary constructed in sfmodel_fit(), all possible keys are: :frontier, :timevar, :log_σ²ᵤ₀, :log_σ²ᵤ_star, :log_σ²₍₀, log_σ²₍_star, log_σ²ω⁰, log_σ²w_star  

# Examples
```julia
frontier_hat = predict(@eq(frontier), df);
sigma_c_0_hat = predict(@eq(log_σ²₍₀), df);
sigma2_c_0_hat = predict(@eq(σ²₍₀), df);
```
""" =#

function predict(eq::Vector{Symbol},  sfdata::DataFrame) 

    eqname = eq[1]

    takeexp = false

    if eqname == :frontier
        eq_var = :frontier
        eq_coe = :frontier
    elseif eqname == :timevar || eqname == :time
        eq_var = :timevar
        eq_coe = :timevar

    elseif eqname == :sigma2_u_0 || eqname == :σ²ᵤ₀
        eq_var  = :σ²ᵤ₀
        eq_coe  = :log_σ²ᵤ₀
        takeexp = true

    elseif eqname == :log_sigma2_u_0 || eqname == :log_σ²ᵤ₀
        eq_var  = :σ²ᵤ₀
        eq_coe  = :log_σ²ᵤ₀

    elseif eqname == :sigma2_u_star || eqname == :σ²ᵤ_star
        eq_var  = :σ²ᵤ_star
        eq_coe  = :log_σ²ᵤ_star
        takeexp = true

    elseif eqname == :log_sigma2_u_star || eqname == :log_σ²ᵤ_star
        eq_var  = :σ²ᵤ_star
        eq_coe  = :log_σ²ᵤ_star

    elseif eqname == :sigma2_c_0 || eqname == :σ²₍₀
        eq_var  = :σ²₍₀
        eq_coe  = :log_σ²₍₀
        takeexp = true
    
    elseif eqname == :log_sigma2_c_0 || eqname == :log_σ²₍₀
        eq_var  = :σ²₍₀
        eq_coe  = :log_σ²₍₀

    elseif eqname == :sigma2_c_star || eqname == :σ²₍_star
        eq_var  = :σ²₍_star
        eq_coe  = :log_σ²₍_star
        takeexp = true
    
    elseif eqname == :log_sigma2_c_star || eqname == :log_σ²₍_star
        eq_var  = :σ²₍_star
        eq_coe  = :log_σ²₍_star

    elseif eqname == :sigma2_w_0 || eqname == :σ²ω⁰
        eq_var  = :σ²ω⁰
        eq_coe  = :log_σ²ω⁰
        takeexp = true
    
    elseif eqname == :log_sigma2_w_0 || eqname == :log_σ²ω⁰
        eq_var  = :σ²ω⁰
        eq_coe  = :log_σ²ω⁰

    elseif eqname == :sigma2_w_star || eqname == :σ²w_star
        eq_var  = :σ²w_star
        eq_coe  = :log_σ²w_star
        takeexp = true
    
    elseif eqname == :log_sigma2_w_star || eqname == :log_σ²w_star
        eq_var  = :σ²w_star
        eq_coe  = :log_σ²w_star
       
    else 
        throw("@eq() is not specified correctly")   
    end

    #* -- begin processing the data --------
    if (eq_var in [:frontier, :timevar])
        pvar = sfdata[:, _dicM[eq_var]]
        pvar = convert(Array{Float64}, pvar)
        value = pvar*_eqncoe[eq_coe]
        return value
    else  
        # if equation is error standard deviation since _cons is only exogenous determinant of them,
        # just return estimation result directly.

        # "Since there is no other exogenous determinant for $eq_var except _cons, prediction for every observation will be same, therefore, just return one prediction")
        takeexp == true || return(_eqncoe[eq_coe])
        
        return exp(_eqncoe[eq_coe])
    end

end