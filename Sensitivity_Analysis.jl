
# Import relevant packages
using Sobol, Plots

# Include the right-hand-side and objective functions
include(joinpath(@__DIR__, "Functions.jl"))

# Do a sensitivity analysis for each parameter starting from the optimized results

# Optimized blackbox results:
e1_res = [0.019325275, 0.001042151, 0.001449269, 4.65E-05, 0.001082172, 1.72E-16, 0, 0.003622651, 0.020911826, 6.57E-05, 0.029629188]
e2_res = [0.030722039, 0.011610798, 3.35E-06, 0.007322919, 0.000480144, 0.001548954, 3.21E-23, 6.96E-20, 1.63E-25, 2.75E-22, 0.791162583, 0.163629988, 0.038254347, 0.291499494, 0.033444937, 6.79E-20, 0.111915253]
e3_res = [0.098836516, 0.006872884, 0.034535793, 0.002792464, 2.58E-05, 3.14E-18, 9.04E-86, 0.096750287, 0.163481954, 0.198059637, 0.012763079]
e4_res = [0.073396763, 0.00416242, 0.001255579, 3.86E-06, 1.57E-19, 0.002260102, 0.20353776, 0.726192307]
e5_res = [0.141863434, 0.091993842, 0.003077546, 0.033697918, 0.019062568, 0.014813437, 8.52E-24, 9.14E-22, 3.85E-22, 1.41E-24, 0.096787124, 0.454861078, 0.049914551, 0.574747917, 0.011027763, 1.70E-05, 0.650334741]
e6_res = [0.213966914, 0.01870633, 0.105161052, 0.003919695, 0.046129702, 1.98E-22, 0.006357671, 0.258936789, 5.05E-15, 0.366243266, 0.022720123]
e7_res = [0.099067857, 0.100899746, 0.007525859, 0.033043242, 0.007350668, 0.000204659, 4.40E-23, 6.11E-21, 7.74E-05, 1.14E-20, 0.019877077, 0.540830275, 2.60E-16, 0.83745931, 1.13E-12, 4.68E-22, 0.999838314]
e8_res = [0.095129397, 0.007828674, 0.071871613, 5.74E-18, 0.007933755, 4.55E-20, 1.57E-17, 1.82E-01, 1.73E-16, 3.98E-01, 2.24E-03]

sensitivities = zeros(8,17)
targets = [e1_res, e2_res, e3_res, e4_res, e5_res, e6_res, e7_res, e8_res]
Sobol_points = 10000
search_dist = 1e-5
local_search_dist = 1e-7

# Now we make the big matrix of sensitivities.
for i = 1:8
    seq = SobolSeq(length(targets[i]))

    lo = max.(0.0, targets[i].-search_dist)
    hi = min.(1.0, targets[i].+search_dist)

    # Search through an arbitrary number of Sobol points
    for j = 1:Sobol_points
        x = next!(seq)
        x_new = x.*(hi-lo)+lo
        # @show x_new

        # Check the local sensitivity of each parameter and add its absolute value
        # to the total
        for k = 1:length(targets[i])
            local_lo = copy(x_new)
            local_hi = copy(x_new)
            local_lo[k] -= local_search_dist
            local_hi[k] += local_search_dist

            if i==1
                sensitivities[i,k] += abs((objective_func_e1(local_hi) - objective_func_e1(local_lo))/(2*local_search_dist))
            elseif i==2
                sensitivities[i,k] += abs((objective_func_e2(local_hi) - objective_func_e2(local_lo))/(2*local_search_dist))
            elseif i==3
                sensitivities[i,k] += abs((objective_func_e3(local_hi) - objective_func_e3(local_lo))/(2*local_search_dist))
            elseif i==4
                sensitivities[i,k] += abs((objective_func_e4(local_hi) - objective_func_e4(local_lo))/(2*local_search_dist))
            elseif i==5
                sensitivities[i,k] += abs((objective_func_e5(local_hi) - objective_func_e5(local_lo))/(2*local_search_dist))
            elseif i==6
                sensitivities[i,k] += abs((objective_func_e6(local_hi) - objective_func_e6(local_lo))/(2*local_search_dist))
            elseif i==7
                sensitivities[i,k] += abs((objective_func_e7(local_hi) - objective_func_e7(local_lo))/(2*local_search_dist))
            elseif i==8
                sensitivities[i,k] += abs((objective_func_e8(local_hi) - objective_func_e8(local_lo))/(2*local_search_dist))
            end
        end
    end

    # Now each one is the sum of the number of Sobol points. Divide by that value to
    # get an average sensitivity
    sensitivities[i,:] = sensitivities[i,:] ./ Sobol_points
end

@show sensitivities
