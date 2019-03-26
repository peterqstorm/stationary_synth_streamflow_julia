
#load libraries
#note that if you need to install a library from Github, use the following two lines
#using Pkg
#Pkg.add("name")
using Base
using CSV
using Tables
using LinearAlgebra
#import functions from local library and adding to PATH
push!( LOAD_PATH, "./" )
using stationary_generator


datadir = pwd() * "/data/"

Qdaily = CSV.read(datadir*"Qdaily.txt", delim=" ")
Qdaily = convert(Matrix, Qdaily)[:, 1:4]

Qdaily[:,4] = log.(Qdaily[:,4])
sites = ["qMarietta", "qMuddyRun", "qLateral", "evapConowingo"]

Nyears = size(Qdaily, 1) / 365
Nsites = size(Qdaily, 2)

num_realizations = [100, 1000]
num_years = [100, 1]

if isdir(datadir*"\\validation\\")
    x=1
else
    mkdir(datadir*"\\validation\\")
end

for k = 1:size(num_realizations, 1)
  Qd_cg = stationary_generator.combined_generator(Qdaily, num_realizations[k], num_years[k])

  #back-transform data

  Qd_cg[:][:, 4] = log.(Qd_cg[:][:, 4])

  #write simulations to file
  # for i = 1:Nsites
  # end


end
