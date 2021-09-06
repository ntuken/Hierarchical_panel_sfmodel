function colprod(mat)
	# input R x C matrix, output 1 X C matrix (each column is column product)
    # mat = convert(Matrix{Float64}, mat)
    prodv = ((-1) .^ (sum(mat .< 0, dims=1)) .* exp.(sum(log.(abs.(mat)), dims=1)))
    return prodv
end

