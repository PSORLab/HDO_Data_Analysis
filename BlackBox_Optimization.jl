
# Import the black-box optimization package
using BlackBoxOptim

# Include the right-hand-side and objective functions
include(joinpath(@__DIR__, "Functions.jl"))

# Loop through each reaction, performing the resampling inhereitance memetic search for 1 hour on each
for i=1:8
    println(i)
    if i==1
        bbsolution = bboptimize(objective_func_e1; SearchRange=(0.0, 1.0), NumDimensions=11, Method=:resampling_inheritance_memetic_search, MaxTime=3600.0)
    elseif i==2
        bbsolution = bboptimize(objective_func_e2; SearchRange=(0.0, 1.0), NumDimensions=17, Method=:resampling_inheritance_memetic_search, MaxTime=3600.0)
    elseif i==3
        bbsolution = bboptimize(objective_func_e3; SearchRange=(0.0, 1.0), NumDimensions=11, Method=:resampling_inheritance_memetic_search, MaxTime=3600.0)
    elseif i==4
        bbsolution = bboptimize(objective_func_e4; SearchRange=(0.0, 1.0), NumDimensions=8, Method=:resampling_inheritance_memetic_search, MaxTime=3600.0)
    elseif i==5
        bbsolution = bboptimize(objective_func_e5; SearchRange=(0.0, 1.0), NumDimensions=17, Method=:resampling_inheritance_memetic_search, MaxTime=3600.0)
    elseif i==6
        bbsolution = bboptimize(objective_func_e6; SearchRange=(0.0, 1.0), NumDimensions=11, Method=:resampling_inheritance_memetic_search, MaxTime=3600.0)
    elseif i==7
        bbsolution = bboptimize(objective_func_e7; SearchRange=(0.0, 1.0), NumDimensions=17, Method=:resampling_inheritance_memetic_search, MaxTime=3600.0)
    elseif i==8
        bbsolution = bboptimize(objective_func_e8; SearchRange=(0.0, 1.0), NumDimensions=11, Method=:resampling_inheritance_memetic_search, MaxTime=3600.0)
    end
    println("result $i is:", best_candidate(bbsolution))
end
   

