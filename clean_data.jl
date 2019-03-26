
# This loads the mat files and combines the Q [cfs] time series into one
# txt file. Removes leap years by averaging with the prior date().

using Tables
using Base


datadir = pwd() * "\\data\\"
files = ["qMarietta_1932-2001.csv", "qMuddyRun_1932-2001.csv",
    "qLateral_1932-2001.csv", "evapConowingo_1932-2001.csv",
    "evapMuddyRun_1932-2001.csv"]

hist_data = []

# find indices of leap years
# this is specific to the Susquehanna; not general
leaps = Int.(60:365*3+366:365*(2001-1932+1)+ceil(ceil(2001-1932)/4))
all = Int.(1:1:365*(2001-1932+1)+ceil(ceil(2001-1932)/4) - 1)
non_leaps = Int.(setdiff(all, leaps))

Qfinal = zeros(Int(length(non_leaps) + 1), length(files))


for i = 1:length(files)
    temp_file = CSV.read(datadir*files[1], header=false)[1]
    println(length(temp_file))
    temp = Any[]
    println(datadir*files[1])
    x = 0
    for index in 1:Int(length(all) + 1)
        if index in leaps
            x += 1
        else
            push!(temp, temp_file[index])
        end
    end
    Qfinal[:, i] = temp

end


CSV.write(datadir*"QDaily.txt", Tables.table(Qfinal), delim=" ")


# reshape into nyears x 365 and nyears x 12 for daily and monthly
# statistical validation figures
Qfinal_monthly = stationary_generator.convert_data_to_monthly(Qfinal)


# divide evaporation by 86400 (s/day) to get total monthly evap in in/month
Qfinal_monthly[4] = Qfinal_monthly[4] ./ 86400
Qfinal_monthly[5] = Qfinal_monthly[5] ./ 86400

# create directories to write files to
if isdir(datadir*"cleaned_data")
    x=1
else
    mkdir(datadir*"cleaned_data")
end

new_dir = datadir*"\\cleaned_data\\"

for i = 1:length(files)
    ff = files[i]
    file_name = new_dir*ff[1:length(ff) - 10]
    println(file_name)
    q_nx365 = reshape(Qfinal[:, i], 365, :)
    CSV.write(file_name*"-daily.csv", Tables.table(q_nx365), delim=" ")
    CSV.write(file_name*"-monthly.csv", Tables.table(Qfinal_monthly[i]), delim=" ")
end
