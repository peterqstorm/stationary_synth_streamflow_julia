
#load libraries
#note that if you need to install a library from Github, use the following two lines
#using Pkg
#Pkg.add("name")
using Base
using CSV
using Tables
using LinearAlgebra
using Statistics
using ElasticArrays
#import functions from local library and adding to PATH
push!( LOAD_PATH, "./" )
using stationary_generator


#######################################
#### USER INPUT FOR OUTPUT FILES#######
#### Iterates through both arrays######
#######################################
num_realizations = [10, 10, 10 , 10  , 100, 100, 100, 1000, 1000, 1000]
num_years =        [1 , 10, 100, 1000, 1  , 10 , 100, 1   , 10  , 100 ]
#######################################
#######################################

# Importing data
datadir     = pwd() * "/data/"
Qdaily      = CSV.read(datadir*"Qdaily.txt", delim=" ")
Qdaily      = convert(Matrix, Qdaily)[:, 1:4]
Qdaily[:,4] = log.(Qdaily[:,4])

sites  = ["qMarietta", "qMuddyRun", "qLateral", "evapConowingo"]
Nyears = size(Qdaily, 1) / 365
Nsites = size(Qdaily, 2)


#check if directory for output is available. Create if it isn't.
if isdir(datadir*"\\output\\")
    x=1
else
    mkdir(datadir*"\\output\\")
end


#Note that @time is used to calculate the time for each function
#the final time will be displayed at end of all outputs
@time for k = 1:length(num_realizations)

    #put data and parameters into generator
    @time Qd_cg = combined_generator(Qdaily, num_realizations[k], num_years[k])

    #back-transform data
    for i = 1:size(Qd_cg, 1)
        Qd_cg[i][:, 4] = log.(Qd_cg[i][:, 4])
    end

    #initialize blank arrays
    Qd2       = zeros(365 * num_years[k] * num_realizations[k])
    q_        = Array{Float64, 2}(undef, num_realizations[k], 365 * num_years[k])
    Q_monthly = zeros(num_realizations[k] * num_years[k], 12)

    for i = 1: Nsites
        # put into array of [realizations, 365*num_yrs]
        for j = 1: num_realizations[k]
            q_[j, :] = Qd_cg[j][:, i]'
        end


        # write to csv for daily
        file_name = datadir * "\\output\\" * sites[i] *string(num_realizations[k]) * "x" * string(num_years[k]) * "_daily.csv"

        CSV.write(file_name, Tables.table(q_), writeheader=false)


        #convert to monthly, then write to csv
        Qd2             = reshape(q_, :, 1)
        Q_monthly[:, :] = convert_data_to_monthly(Qd2)[1]
        file_name       = datadir * "\\output\\" * sites[i] *string(num_realizations[k]) * "x" * string(num_years[k]) * "_monthly.csv"

        CSV.write(file_name, Tables.table(Q_monthly), writeheader=false)

        println("Finished $(num_realizations[k]) realizations over $(num_years[k]) years for $(sites[i]).")

    end

    #output entire daily record to csv
    file_name   = datadir * "\\output\\" * "Qdaily" * string(num_realizations[k]) * "x" * string(num_years[k]) * ".csv"
    output_file = deepcopy(Qd_cg[1])

    for i = 2:size(Qd_cg, 1)
        output_file = [output_file; Qd_cg[i]]
    end

    CSV.write(file_name, Tables.table(output_file), writeheader=false)
end
