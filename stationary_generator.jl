# all code adapted from Matteo Giuliani, Jon Herman and Julianne Quinn
# at https://github.com/julianneq/Kirsch-Nowak_Streamflow_Generator
# specifically within this directory: https://github.com/julianneq/Kirsch-Nowak_Streamflow_Generator/tree/master/stationary_generator
# most comments are copied from the code on the github


module stationary_generator
    using LinearAlgebra
    using Statistics
    function cholesky_corr(Z)
        # Computes the cholesky decomp of correlation matrix of columns of Z
        # Then attempts to repair non-positive-definite matrics
        # Code adapted from https://github.com/julianneq/Kirsch-Nowak_Streamflow_Generator/blob/master/stationary_generator/chol_corr.m
        # http://www.mathworks.com/matlabcentral/answers/6057-repair-non-positive-definite-correlation-matrix
        # rank-1 update followed by rescaling to get unit diagonal entries


        R = Statistics.cor(Z)
        U = cholesky(R, check=true)

        #check if positive definite, otherwise modify slightly until true
        while issuccess(U) == false
            k = min([real(eigh(R)) - 1 * eps()])
            R = R - k * Matrix{Float64}(I, size(R), size(R))
            R = R / R[1, 1]
            U = chol(R)
        end

        return U
    end


    function convert_data_to_monthly(Qt)

    
        num_years = size(Qt, 1) / 365   #first dimension of input array
        num_sites = size(Qt, 2)

        num_months = 12
        days_in_each_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        #potential area to improve speed
        Qmonthly = Any[]
        for i = 1:num_sites
           push!(Qmonthly, zeros(Float64, convert(Int, round(num_years)), convert(Int, round(num_months))))
        end

        for yr = 1:num_years
            for mo = 1:num_months

                start = convert(Int, (yr - 1) * 365 + sum(days_in_each_month[1:(mo - 1)]) + 1)

                for i = 1:num_sites
                    Qmonthly[i][convert(Int, yr), mo] = 86400 * sum(Qt[start:start + days_in_each_month[mo] - 1, i])
                end

            end
        end

        return Qmonthly
    end


    function KNN_identification(Z, Qtotals, month)
        # [KNN_id, W] = KNN_identification[ Z, Qtotals, month, k ]
        #
        # Identification of K-nearest neighbors of Z in the historical annual data
        # z and computation of the associated weights W.
        #
        # Input:    Z = synthetic datum [scalar]
        #           Qtotals = total monthly flows at all sites for all historical months 
        #             within +/- 7 days of the month being disaggregated
        #           month = month being disaggregated
        #           k = number of nearest neighbors (by default k=n_year^0.5
        #             according to Lall and Sharma [1996])
        # Output:   KNN_id = indices of the first K-nearest neighbors of Z in the
        #             the historical annual data z
        #           W = nearest neighbors weights, according to Lall and Sharma
        #             (1996): W[i] = (1/i) / (sum(1/i)) 
        #
        #

        # Ntotals is the number of historical monthly patterns used for disaggregation.
        # A pattern is a sequence of ndays of daily flows, where ndays is the
        # number of days in the month being disaggregated. Patterns are all()
        # historical sequences of length ndays beginning within 7 days before or
        # after the 1st day of the month being disaggregated.    
    
    
        n_totals = size(Qtotals[month], 1)
        k = round(sqrt(n_totals))

    
        # nearest neighbors identification
        # only look at neighbors from the same month +/- 7 days
        n_sites = size(Qtotals[month], 1)
        delta = zeroes([n_totals, 1])     # first and last month have 7 less possible shifts

        for i = 1:n_totals
            for j = 1:n_sites
                delta[i] = delta[i] + (Qtotals[month][i, j] - Z[1, 1, j]) ^ 2
            end
        end

        Y = [collect(1:size(delta, 1))', delta]
        sort!(Y, by = x -> x[1])

        KNN_id = Y[1:k, 1]

    
        # computation of the weights
        f = [1:k]
        f = 1 ./ f
        weights = v ./ sum(f)

        return KNN_id, weights
    end
    

    function KNN_sampling(KNN_id, indices, weights_cumm, Qdaily, month)
        # py = KNN_sampling[ KKN_id, indices, Wcum, Qdaily, month ]
        #
        # Selection of one KNN according to the probability distribution defined by
        # the weights W.
        #
        # Input:    KNN_id = indices of the first K-nearest neighbors
        #           indices = n x 2 matrix where n is the number of monthly totals
        #             and the 2 columns store the historical year in which each
        #             monthly total begins, and the number of shift index
        #             where 1 is 7 days earlier and 15 is 7 days later
        #           Wcum = cumulated probability for each nearest neighbor
        #           Qdaily = historical data
        #           month = month being disaggregated
        # Output:   py = selected proportion vector corresponding to the sampled
        #             shifted historical month
        #           yearID = randomly selected monthly total [row to select from indices]
        #
        # 

        #Randomly select one of the k-NN using the Lall and Sharma density
        #estimator
    
    
        r = rand()
        prepend!(0, weights_cumm)

        for i = 1:length(weights_cumm)-1
            if(r > weights_cumm[i]) && (r <= weights_cumm[i + 1])
                KNNs = i
            end      
        end
        yearID = KNN_id[KNNs]

        # concatenate last 7 days of last year before first 7 days of first year
        # and first 7 days of first year after last 7 days of last year
        nrows = size(Qdaily, 1)
        QDaily = [Qdaily[nrows-7:nrows,:]; Qdaily; Qdaily[1:8,:]]

    
        #shift historical data to get nearest neighbor corresponding to yearID
        year = indices(yearID, 1)
        k = indices(yearID, 2)
        shifted_Qdaily = Qdaily[k:k + nrows - 1, :]


        days_in_each_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        start = (year - 1) * 365 + sum(days_in_each_month[1:(month - 1)]) + 1
        daily_flows = shifted_Qdaily[start:start + days_in_each_month[month] - 1, :]

        py = zeros(size(daily_flows))
        for i = 1:size(QDaily, 2)
            py[:, i] = daily_flows[:, i] / sum(daily_flows[:, i])
        end

        return py, yearID
    end


    function monthly_gen(Q_historical, num_years)

        num_points = length(Q_historical)
        num_Q_hist = length(Q_historical[1][:,1])
        
        #error checking
        for i = 2:num_points
            if length(Q_historical[i][:,1]) != num_Q_hist
                error("All matrices in Q_historical must be the same size.")
            end
        end
            
        num_years = num_years + 1      #adjusts for their new corr technique
    
        nQ = num_Q_hist
        random_matrix = rand(1:nQ, num_years, 12)

        Qs = Any[]

        for k = 1: num_points
            Q_matrix = Q_historical[k]

            logQ = log.(Q_matrix)

            monthly_mean  = zeros(1, 12)
            monthly_stdev = zeros(1, 12)
            Z             = zeros(nQ, 12)


            for i = 1:12
                monthly_mean[i]  = mean(logQ[:, i])
                monthly_stdev[i] = Statistics.std(logQ[:, i])
                Z[:,i] = (logQ[:, i] .- monthly_mean[i]) ./ monthly_stdev[i]
            end

            Z_vector = reshape(Z', 1, :)
            Z_shifted = reshape(Z_vector[7:(nQ * 12 - 6)], 12, :)
           
        
            # The correlation matrices should use the historical Z's
            # (the "appended years" do not preserve correlation)
            U = cholesky_corr(Z[1:num_Q_hist, :])
            U_shifted = cholesky_corr(Z_shifted[1: num_Q_hist - 1, :])

            Qs_uncorr = Any[]
            for i = 1:12
                push!(Qs_uncorr, Z[random_matrix[:, 3][3]])
            end


            Qs_uncorr_vector        = reshape(Qs_uncorr[:, :]', 1, :)
            Qs_uncorr_shifted[:, :] = reshape(Qs_uncorr_vector[7:(num_years * 12 - 6)], 12, :)'
            Qs_corr[:, :]           = Qs_uncorr[:, :] * U
            Qs_corr_shifted         = Qs_uncorr_shifted[:, :] * U_shifted

            Qs_log[:, 1:6]  = Qs_corr_shifted[:, 7:12]
            Qs_log[:, 7:12] = Qs_corr[2:num_years, 7:12]


            Qsk = Any[]
            for i = 1:12
                push!(Qsk, exp.(Qs_log[:, i] .* monthly_stdev[i] .+ monthly_mean[i]))
            end

            push!(Qs, Qsk)
        end


       return Qs 
    end


    function monthly_main( hist_data, nR, nY )

        num_years  = size(hist_data, 1) / 365
        num_sites = size(hist_data, 2)
        
        # from daily to monthly
        Qh = convert_data_to_monthly(hist_data)
        
        # initialize output
        qq = Any[]
        for i = 1:num_sites
           push!(qq, zeros(nR, nY * 12))
        end

    
        # generate data
        for r = 1:nR
            Qs = monthly_gen(Qh, nY)
            for k = 1:num_sites
                qq[k][r, :] = reshape(Qs[k]', 1, :)
            end
        end

        #output matrix
        Qgen = Any[]

        for k = 1:num_sites
           push!(Qgen, qq[k]) 
        end

        return Qgen
    end


    function combined_generator(hist_data, nR, nY)

        num_sites = size(hist_data, 2)

    
        # generation of monthly data via Kirsch et al. (2013):
        # Kirsch, B. R., G. W. Characklis, and H. B. Zeff [2013], 
        # Evaluating the impact of alternative hydro-climate scenarios on transfer 
        # agreements: Practical improvement for generating synthetic streamflows, 
        # Journal of Water Resources Planning and Management, 139[4], 396–406.
        QQg       = monthly_main(hist_data, nR, nY)
        Qh        = convert_data_to_monthly(hist_data)
        num_years = size(Qh[1], 1)

    
    
        # disaggregation from monthly to daily time step as in Nowak et al. (2010):
        # Nowak, K., Prairie, J., Rajagopalan, B., & Lall, U. (2010). 
        # A nonparametric stochastic approach for multisite disaggregation of 
        # annual to daily streamflow. Water Resources Research, 46[8].

        # Find K-nearest neighbors [KNN] in terms of total monthly flow and 
        # randomly select one for disaggregation. Proportionally scale the flows in
        # the selected neighbor to match the synthetic monthly total. To
        # disaggregate Jan Flows, consider all historical January totals +/- 7
        # days, etc.
        Dt = 3600 * 24
        days_in_each_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        D = zeros(nR, 365 * nY, num_sites)

    
        # concatenate last 7 days of last year before first 7 days of first year
        # and first 7 days of first year after last 7 days of last year
        nrows = size(hist_data, 1)
        extra_hist_data = [hist_data[nrows-7:nrows,:]; hist_data; hist_data[1:8,:]]

        Qtotals  = Any[]
        Qindices = Any[]

    
        # find monthly totals for all months +/- 7 days
        for i = 1:12
            count = 1

            if i == 1 || i == 12
                nTotals = num_years * 15 - 7     # 7 fewer shifts in first and last month
            else
                nTotals = num_years * 15
            end

            Qmonthly_shifted = zeroes(nTotals, num_sites)
            indices          = zeros(nTotals, 2)

            for k = 1:15
                shifted_hist_data = extra_hist_data[k: k + nrows - 1, :]
                Qh = convert_data_to_monthly(shifted_hist_data)

                for j = 1:num_sites
                    if i == 1 && k < 8
                        Qh[j] = Qh[j][2:size(Qh[j], 1), i]     # remove first year
                    elseif i == 12 && k > 8
                        Qh[j] = Qh[j][1:size(Qh[j], 1) - 1, i] # remove last year
                    end
                    Qmonthly_shifted[count:(count + size(Qh[j], 1) - 1), 1] = Qh[j][:, 1]
                end

                if i == 1 && k < 8
                    indices[count: (count + size(Qh[j], 1) - 1), 1] = 2:(size(Qh[j], 1) + 1)
                else
                    indices[count: (count + size(Qh[j], 1) - 1), 1] = 1:(size(Qh[j], 1))
                end
                indices[count:(count + size(Qh[j], 1) - 1), 2]

                count = count + size(Qh[j], 1)
            end

            !push(Qtotals, Qmonthly_shifted)
            !push(Qindices, indices)

        end



        for r = 1:nR
            dd = Any[]
            for i = 1:Ny * 12
                # monthly value for all sites
                Z = QQg(r, i, :)
                    
                #KNN and weights
                month = mod(i, 12)
                if month == 0
                    month = 12
                end
                KNN_id, W = KNN_identification(Z, Qtotals, month)
                Wcum = cumsum(W)

                   #sampling of one KNN
                py = KNN_sampling(KNN_id, Qindices[month], Wcum, hist_data, month)
                d = zeros(num_sites, days_in_each_month[month])

                for j = 1:num_sites
                    d[j, :] = py[:, j] .* Z[1, 1, j]
                end

                !push(dd, d)
            end



            D[r, :, :] = dd' ./ Dt
        end



        return D

    end

end