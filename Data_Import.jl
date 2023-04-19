
# This file loads data from the CSV data file into the form of a RxnData struct,
# for data analysis in later files.
using CSV, Statistics, DataFrames

# Create a quick data structure to store relevant information
mutable struct RxnData
    T::Float64 #Temperature [degC]
    P::Float64 #Pressure [psi] (PSIA or PSIG?)
    solvent::String #Either water or decane
    solv_vol::Float64 #Solvent volume [mL]
    cat::String #Catalyst type (full identifier)
    cat_metal::String #Metal type in catalyst
    cat_supp::String #Catalyst support
    cat_wt::Float64 #Catalyst weight [g]
    time::Vector{Float64} #Time [minutes]
    c_guaiacol::Vector{Float64} #guiacol concentration [mol/L]
    c_methoxycyclohexanone::Vector{Float64} #methoxycyclohexanone concentration [mol/L]
    c_methoxycyclohexane::Vector{Float64} #methoxycyclohexane concentration [mol/L]
    c_phenol::Vector{Float64} #phenol concentration [mol/L]
    c_cyclohexanone::Vector{Float64} #cyclohexanone concentration [mol/L]
    c_cyclohexanol::Vector{Float64} #cyclohexanol concentration [mol/L]
    c_cyclohexane::Vector{Float64} #cyclohexane concentration [mol/L]
    c_coke::Vector{Float64} #Coked concentration, equal to starting_conc-sum(all other concs)
    # Standard deviations for all species:
    s_guaiacol::Vector{Float64} #STDEV of guiacol concentration [mol/L]
    s_methoxycyclohexanone::Vector{Float64} #STDEV of methoxycyclohexanone concentration [mol/L]
    s_methoxycyclohexane::Vector{Float64} #STDEV of methoxycyclohexane concentration [mol/L]
    s_phenol::Vector{Float64} #STDEV of phenol concentration [mol/L]
    s_cyclohexanone::Vector{Float64} #STDEV of cyclohexanone concentration [mol/L]
    s_cyclohexanol::Vector{Float64} #STDEV of cyclohexanol concentration [mol/L]
    s_cyclohexane::Vector{Float64} #STDEV of cyclohexane concentration [mol/L]
    s_coke::Vector{Float64} #STDEV of Coked concentration, equal to starting_conc-sum(all other concs)
end

# Import data from the data CSV file
data = Matrix{Any}(CSV.read(joinpath(@__DIR__, "HDO_Data.csv"), DataFrame))

# Set up lists of catalysts and reaction species
catalyst_list = ["10%Ni-CAC in water",
                 "10%Ni-FWAC-meso in decane",
                 "5%Ru-CAC in decane",
                 "5%Ru-CAC in water",
                 "5%Ru-FWAC-meso in decane",
                 "5%Ru-FWAC-meso in water",
                 "5%Ru-FWAC-micro in decane",
                 "5%Ru-FWAC-micro in water"]

species_list = ["guaiacol",
                "methoxycyclohexanone",
                "methoxycyclohexane",
                "phenol",
                "cyclohexanone",
                "cyclohexanol",
                "cyclohexane"]

# Loop through the data using lists of catalysts and reaction species
for i = eachindex(catalyst_list)
    # Extract only data for the current catalyst
    temp_data = data[(data[:,1].==catalyst_list[i]),:]

    # Trim out the missing's
    temp_data = temp_data[:, 1:(size(temp_data)[2]-sum(ismissing.(temp_data)[1,:]))]

    # Get the times, and create a concentration list that we will update with species averages
    time_vector = convert(Vector{Float64}, temp_data[1,3:end])
    temp_concs = zeros(length(species_list)+1, size(temp_data)[2]-2)
    temp_stds = zeros(length(species_list)+1, size(temp_data)[2]-2)

    # Loop through species and adjust temp_concs based on the species average
    for j = 1:length(species_list)
        species_data = temp_data[(temp_data[:,2].==species_list[j]),3:end]
        if size(species_data)[1] != 0
            temp_concs[j,:] = sum(species_data, dims=1)./3
        end
        if ~isempty(species_data)
            temp_stds[j,:] = Statistics.std(species_data, dims=1)
        end
    end

    # Create the "coking" amount
    start_point = temp_concs[1,1]
    temp_concs[end,:] = start_point .- sum(temp_concs[1:end-1, :], dims=1)
    temp_stds[end,:] = sqrt.(sum(temp_stds[1:end-1,:].^2, dims=1))

    # Fill in the RxnData data structure..?
    if catalyst_list[i]=="10%Ni-CAC in water"
        global e1 = RxnData(-1.0, -1.0, "water", -1.0, "10%Ni-CAC in water", "Ni", "CAC", -1.0, time_vector, temp_concs[1,:], temp_concs[2,:], temp_concs[3,:],
                            temp_concs[4,:], temp_concs[5,:], temp_concs[6,:], temp_concs[7,:], temp_concs[8,:], temp_stds[1,:], temp_stds[2,:], temp_stds[3,:],
                            temp_stds[4,:], temp_stds[5,:], temp_stds[6,:], temp_stds[7,:], temp_stds[8,:])
    elseif catalyst_list[i]=="5%Ru-CAC in decane"
        global e2 = RxnData(-1.0, -1.0, "decane", -1.0, "5%Ru-CAC in decane", "Ru", "CAC", -1.0, time_vector, temp_concs[1,:], temp_concs[2,:], temp_concs[3,:],
                            temp_concs[4,:], temp_concs[5,:], temp_concs[6,:], temp_concs[7,:], temp_concs[8,:], temp_stds[1,:], temp_stds[2,:], temp_stds[3,:],
                            temp_stds[4,:], temp_stds[5,:], temp_stds[6,:], temp_stds[7,:], temp_stds[8,:])
    elseif catalyst_list[i]=="5%Ru-CAC in water"
        global e3 = RxnData(-1.0, -1.0, "water", -1.0, "5%Ru-CAC in water", "Ru", "CAC", -1.0, time_vector, temp_concs[1,:], temp_concs[2,:], temp_concs[3,:],
                            temp_concs[4,:], temp_concs[5,:], temp_concs[6,:], temp_concs[7,:], temp_concs[8,:], temp_stds[1,:], temp_stds[2,:], temp_stds[3,:],
                            temp_stds[4,:], temp_stds[5,:], temp_stds[6,:], temp_stds[7,:], temp_stds[8,:])
    elseif catalyst_list[i]=="10%Ni-FWAC-meso in decane"
        global e4 = RxnData(-1.0, -1.0, "decane", -1.0, "10%Ni-FWAC-meso in decane", "Ni", "FWAC-meso", -1.0, time_vector, temp_concs[1,:], temp_concs[2,:], temp_concs[3,:],
                            temp_concs[4,:], temp_concs[5,:], temp_concs[6,:], temp_concs[7,:], temp_concs[8,:], temp_stds[1,:], temp_stds[2,:], temp_stds[3,:],
                            temp_stds[4,:], temp_stds[5,:], temp_stds[6,:], temp_stds[7,:], temp_stds[8,:])
    elseif catalyst_list[i]=="5%Ru-FWAC-meso in decane"
        global e5 = RxnData(-1.0, -1.0, "decane", -1.0, "5%Ru-FWAC-meso in decane", "Ru", "FWAC-meso", -1.0, time_vector, temp_concs[1,:], temp_concs[2,:], temp_concs[3,:],
                            temp_concs[4,:], temp_concs[5,:], temp_concs[6,:], temp_concs[7,:], temp_concs[8,:], temp_stds[1,:], temp_stds[2,:], temp_stds[3,:],
                            temp_stds[4,:], temp_stds[5,:], temp_stds[6,:], temp_stds[7,:], temp_stds[8,:])
    elseif catalyst_list[i]=="5%Ru-FWAC-meso in water"
        global e6 = RxnData(-1.0, -1.0, "water", -1.0, "5%Ru-FWAC-meso in water", "Ru", "FWAC-meso", -1.0, time_vector, temp_concs[1,:], temp_concs[2,:], temp_concs[3,:],
                            temp_concs[4,:], temp_concs[5,:], temp_concs[6,:], temp_concs[7,:], temp_concs[8,:], temp_stds[1,:], temp_stds[2,:], temp_stds[3,:],
                            temp_stds[4,:], temp_stds[5,:], temp_stds[6,:], temp_stds[7,:], temp_stds[8,:])
    elseif catalyst_list[i]=="5%Ru-FWAC-micro in decane"
        global e7 = RxnData(-1.0, -1.0, "decane", -1.0, "5%Ru-FWAC-micro in decane", "Ru", "FWAC-micro", -1.0, time_vector, temp_concs[1,:], temp_concs[2,:], temp_concs[3,:],
                            temp_concs[4,:], temp_concs[5,:], temp_concs[6,:], temp_concs[7,:], temp_concs[8,:], temp_stds[1,:], temp_stds[2,:], temp_stds[3,:],
                            temp_stds[4,:], temp_stds[5,:], temp_stds[6,:], temp_stds[7,:], temp_stds[8,:])
    elseif catalyst_list[i]=="5%Ru-FWAC-micro in water"
        global e8 = RxnData(-1.0, -1.0, "water", -1.0, "5%Ru-FWAC-micro in water", "Ru", "FWAC-micro", -1.0, time_vector, temp_concs[1,:], temp_concs[2,:], temp_concs[3,:],
                            temp_concs[4,:], temp_concs[5,:], temp_concs[6,:], temp_concs[7,:], temp_concs[8,:], temp_stds[1,:], temp_stds[2,:], temp_stds[3,:],
                            temp_stds[4,:], temp_stds[5,:], temp_stds[6,:], temp_stds[7,:], temp_stds[8,:])
    end
end