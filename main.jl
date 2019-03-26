
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

@time for k = 1:length(num_realizations)
    @time Qd_cg = stationary_generator.combined_generator(Qdaily, num_realizations[k], num_years[k])

    #back-transform data
    for i = 1:size(Qd_cg, 1)
        Qd_cg[i][:, 4] = log.(Qd_cg[i][:, 4])
    end

    Qd2 = zeros(365 * num_years[k] * num_realizations[k])
    q_ = zeros(num_realizations[k], 365 * num_years[k])
    Q_monthly = zeros(num_realizations[k] * num_years[k], 12)

    for i = 1: Nsites
        #put into array of [realizations, 365*num_yrs]
        for j = 1: num_realizations[k]

            q_[i] = Qd_cg[j][:, i]'[1]
        end

        #write to csv for daily
        file_name = datadir * "\\validation\\" * sites[i] *string(num_realizations[k]) * "x" * string(num_years[k]) * "_daily.csv"
        CSV.write(file_name, Tables.table(q_), writeheader=false)

        #convert to monthly, then write to csv
        Qd2 = reshape(q_, :, 1)
        Q_monthly[:, :] = stationary_generator.convert_data_to_monthly(Qd2)[1]
        file_name = datadir * "\\validation\\" * sites[i] *string(num_realizations[k]) * "x" * string(num_years[k]) * "_monthly.csv"
        CSV.write(file_name, Tables.table(Q_monthly), writeheader=false)
    end

end
