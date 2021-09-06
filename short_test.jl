using CSV, DataFrames
data = CSV.read("sim_data1.csv", DataFrame, missingstring="NA", copycols=true)

first(data,6)
tvar = data[:,:time_id]
println(typeof(tvar))
unique_T = unique(tvar)
println(unique_T)
Tvar = DataFrame()  # intitialize dummy time dataframe, not matrix since need to reserve time information in column label
# e.g.  
#       t = [            Tvar = t: 1   2   3
#            1
#            2                   [ 1   0   0
#            3          ==>        0   1   0   
#            1                     0   0   1  
#             ]                    1   0   0 ]

for Time in unique_T
    println(Time)
    Tvar[!,"$Time"] = Float64.(tvar .== Time)
end