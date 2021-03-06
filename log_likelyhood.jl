# using Random, Distributions, Statistics package
using HaltonSequences

""" set random seed and nofdraw here !!  """
Random.seed!(val_random_seed)

""" create  HaltonSequences  """

nofdraw = 100
randseq_tmp = HaltonPoint(4, start=20, length=nofdraw)  # dimension is 4(w_0, w_star, c_0, c_star)

randseq = zeros(nofdraw, 4)
for i in 1:nofdraw
    randseq[i,:] = randseq_tmp[i]
end


function LL_T(Y, X, T,
   po, rho,   nofobsinG_list , noffirm_list, noftime_list, nofdraw=nofdraw, halton=randseq)
    # X may include const  
    
    β  = rho[1:po.endx]
    ϴ = rho[po.begT:po.endT]
    δ1 = rho[po.pos_log_σ²ᵤ₀]  # σ²_u⁰ = exp(δ2)
    δ2 = rho[po.pos_log_σ²ᵤ_star]

    γ1 = rho[po.pos_log_σ²₍₀]
    γ2 = rho[po.pos_log_σ²₍_star]

    ϝ1 = rho[po.pos_log_σ²w⁰]
    ϝ2 = rho[po.pos_log_σ²w_star]
    
     
    sigma2_u_0 = exp(δ1)     
    sigma2_u_star = exp(δ2)
    sigma2_c_0 = exp(γ1)
    sigma2_c_star = exp(γ2)
    sigma2_w_0 = exp(ϝ1)
    sigma2_w_star = exp(ϝ2)

    sigma_u_0 = exp(0.5*δ1)     
    sigma_u_star = exp(0.5*δ2)
    sigma_c_0 = exp(0.5*γ1)
    sigma_c_star = exp(0.5*γ2)
    sigma_w_0 = exp(0.5*ϝ1)
    sigma_w_star = exp(0.5*ϝ2)

    sigma_u = sqrt(sigma2_u_0 + sigma2_u_star)   
    lambda_u     = -sigma_u_star / sigma_u_0
    
    my_w_g = (quantile.(Normal(0,1), randseq[:,1]) .*sigma_w_0 ) .- (quantile.(Normal(0,1), (0.5 .* randseq[:,2]) .+ 0.5) .* sigma_w_star)
    my_c_i = (quantile.(Normal(0,1), randseq[:,3]) .*sigma_c_0 ) .- (quantile.(Normal(0,1), (0.5 .* randseq[:,4]) .+ 0.5) .* sigma_c_star)


    lnLike = 0.0  # initilize lnlike 

    g_row_beg = 1
    
    # nofobsinG, noffirm, noftime
    firm_index = 1
    for g in 1:length(nofobsinG_list)
        
        nofobsinG = nofobsinG_list[g]
        noffirm = noffirm_list[g]
        g_row_end = g_row_beg + nofobsinG - 1
        like_gg = Array{Real}(undef, nofdraw, nofdraw)  # initilize like_gg vector
        # given group g, retrieve all rows in group g
        block_g = (Y[g_row_beg:g_row_end] .- (X[g_row_beg:g_row_end, :] * β ) .- (T[g_row_beg:g_row_end,:] * ϴ)) .*ones(nofobsinG, nofdraw) .- (my_w_g)'
        # (Y[row_g] .- (X[row_g, :] * β ) .- (T[row_g,:] * ϴ)) .*ones(nofobsinG, nofdraw) is NT X S_w matrix
        # (my_w_g)' is transpose of my_w_g, and is 1 X S_w matrix, 
        # block_g is NT X S_w matrix, each column is residual minus different draw from my_w_g

        g_row_beg = g_row_end + 1
        
        Like_my2 = Array{Real}(undef, noffirm, nofdraw)  # initilize with Real type since ForwardDiff require

        id_row_beg = 1
        Like_my2_index = 1
        for id in 1:noffirm
            noftime = noftime_list[firm_index]
            id_row_end = id_row_beg + noftime - 1
            cexpand = kronecker(ones(noftime, nofdraw), (my_c_i)')
            # cexpand is T X(S X S) matrix, given same colummn, the value would be same
            
            u_it = kronecker(block_g[id_row_beg:id_row_end,:], ones(1, nofdraw)) - cexpand
            #  kronecker(block_g[id_row, :], ones(1, nofdraw)) is T X (S_w X S_c) matrix
            #  (cexpand) is T X (S_w X S_c) matrix

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