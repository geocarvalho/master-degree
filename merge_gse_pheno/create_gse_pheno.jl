using CSV
using DataFrames

# reading the file
df = DataFrame(CSV.File("/home/watson/george/master-degree/GSE109381/GSE109381_bvalues_gbm.csv"))
rename!(df, Dict(:Column1 => :samples))