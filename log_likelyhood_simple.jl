# using Random, Distributions, Statistics package

function LL_T(Y, X,
   po, rho, nofG, noffirm, noftime, nofdraw=nofdraw, halton=randseq)
    # X may include const  
    
    β  = rho[1:po.endx]
    log_σ²ᵤ₀ = rho[po.pos_log_σ²ᵤ₀]  # σ²_u⁰ = exp(δ2)
    log_σ²ᵤ_star = rho[po.pos_log_σ²ᵤ_star]

    log_σ²₍₀ = rho[po.pos_log_σ²₍₀]
    log_σ²₍_star = rho[po.pos_log_σ²₍_star]

    log_σ²w₀ = rho[po.pos_log_σ²w⁰]
    log_σ²w_star = rho[po.pos_log_σ²w_star]
    
     
    sigma2_u_0 = exp(log_σ²ᵤ₀)     
    sigma2_u_star = exp(log_σ²ᵤ_star)
    sigma2_c_0 = exp(log_σ²₍₀)
    sigma2_c_star = exp(log_σ²₍_star)
    sigma2_w_0 = exp(log_σ²w₀)
    sigma2_w_star = exp(log_σ²w_star)

    sigma_u_0 = exp(0.5*log_σ²ᵤ₀)     
    sigma_u_star = exp(0.5*log_σ²ᵤ_star)
    sigma_c_0 = exp(0.5*log_σ²₍₀)
    sigma_c_star = exp(0.5*log_σ²₍_star)
    sigma_w_0 = exp(0.5*log_σ²w₀)
    sigma_w_star = exp(0.5*log_σ²w_star)

    sigma_u = sqrt(sigma2_u_0 + sigma2_u_star)   
    lambda_u     = -sigma_u_star / sigma_u_0
    
    my_w_g = (my_w_g_1 .*sigma_w_0 ) .- (my_w_g_2 .* sigma_w_star)
    my_c_i = (my_c_i_1 .*sigma_c_0 ) .- (my_c_i_2 .* sigma_c_star)
    cexpand = ones(noftime, nofdraw) ⊗ (my_c_i)'


    lnLike = 0.0  # initilize lnlike
    Like_my2 = Array{Real}(undef, noffirm, nofdraw)  # initilize here since noffirm and nofdraw are constant

    nofobsinG = noftime*noffirm
    g_row_beg = 1
    
    # Below assume noffirm, noftime are constant, therefore, nofobsinG=noftime*noffirm are constant
    @inbounds for g in 1:nofG
        g_row_end = g_row_beg + nofobsinG - 1
        # given group g, retrieve all rows in group g
        block_g = ((Y[g_row_beg:g_row_end] .- (X[g_row_beg:g_row_end, :] * β ) ) .*ones(nofobsinG, nofdraw)) .- (my_w_g)'
        # (Y[row_g] .- (X[row_g, :] * β ) .- (T[row_g,:] * ϴ)) .*ones(nofobsinG, nofdraw) is NT X S_w matrix
        # (my_w_g)' is transpose of my_w_g, and is 1 X S_w matrix, 
        # block_g is NT X S_w matrix, each column is residual minus different draw from my_w_g

        g_row_beg = g_row_end + 1

        id_row_beg = 1
        Like_my2_index = 1

        @inbounds for id in 1:noffirm
            id_row_end = id_row_beg + noftime - 1
            
            u_it = (block_g[id_row_beg:id_row_end,:]⊗ones(1, nofdraw)) - cexpand
            #  kronecker(block_g[id_row, :], ones(1, nofdraw)) is T X (S_w X S_c) matrix
            #  (cexpand) is T X (S_w X S_c) matrix, given same colummn, the value would be same

            id_row_beg = id_row_end + 1
            Like_i_mat = (2/sigma_u) .* pdf(Normal(0,1), u_it ./ sigma_u) .* cdf(Normal(0,1), lambda_u .* u_it ./ sigma_u)

            #  Likelihood values of observations within a firm, T x (S_w * S_c)  matrix   
            
            
            Like2_ii =  exp.(sum(log.(Like_i_mat), dims=1))  # 1 X (S_w*S_c) Matrix  
            Like_ii = reshape(Like2_ii, nofdraw, nofdraw)
            #  e.g. Like2_ii = [ 1         ,2        ,3        ,4]      
            #                  -w_1-c_1  -w_1-c_2    -w_2-c_1   -w_2-c_2
            # then: Like_ii =  [1 3;
            #                   2 4]
            Like_my1 = mean(Like_ii, dims=1)  # mean over S_c simulation draw

            Like_my2[Like_my2_index,:] = Like_my1        # Like_my2 is N X s_w matrix
            Like_my2_index += 1
            
        end  # loop over id


        Like_gg = exp.(sum(log.(Like_my2), dims=1))  # Like_gg is vector(len = S_w), colummn product of Like_my2

        ans = log(mean((Like_gg)))  # given a group, ans is vector with length 1
        lnLike += ans 
        
    end     # loop over group




    return -lnLike
end