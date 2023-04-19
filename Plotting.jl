
# Note: Need PyPlot installed
using XLSX, Plots

# Import all the data
data = XLSX.readxlsx(joinpath(@__DIR__, "Plot_data.xlsx"))
Ru_CAC_d = data["Ru_CAC_d"][:]
Ru_FWAC_meso_d = data["Ru_FWAC_meso_d"][:]
Ru_FWAC_micro_d = data["Ru_FWAC_micro_d"][:]
Ni_FWAC_meso_d = data["Ni_FWAC_meso_d"][:]
Ru_CAC_w = data["Ru_CAC_w"][:]
Ru_FWAC_meso_w = data["Ru_FWAC_meso_w"][:]
Ru_FWAC_micro_w = data["Ru_FWAC_micro_w"][:]
Ni_CAC_w = data["Ni_CAC_w"][:]

# Now we have matrices of the data. All of the data is formatted the same in terms =
# of columns, so we can just save the column names and cut them all down by one row
cols = Ru_CAC_d[1,:]

Ru_CAC_d = Ru_CAC_d[2:end,:]
Ru_FWAC_meso_d = Ru_FWAC_meso_d[2:end,:]
Ru_FWAC_micro_d = Ru_FWAC_micro_d[2:end,:]
Ni_FWAC_meso_d = Ni_FWAC_meso_d[2:end,:]
Ru_CAC_w = Ru_CAC_w[2:end,:]
Ru_FWAC_meso_w = Ru_FWAC_meso_w[2:end,:]
Ru_FWAC_micro_w = Ru_FWAC_micro_w[2:end,:]
Ni_CAC_w = Ni_CAC_w[2:end,:]

# Now we can plot, based on the columns in "cols" as our known points.
Plots.pyplot()

# Create the decane plot
p1 = plot()
p2 = plot()
p3 = plot()
p4 = plot()
for i = 1:4
    if i==1
        plt = p1
        x = Ru_CAC_d
        stop = 9
    elseif i==2
        plt = p2
        x = Ru_FWAC_meso_d
        stop = 11
    elseif i==3
        plt = p3
        x = Ru_FWAC_micro_d
        stop = 10
    elseif i==4
        plt = p4
        x = Ni_FWAC_meso_d
        stop = 9
    end

    # Define the colors
    c = [:cornflowerblue, # Guaiacol
         :crimson,        # Methoxycyclohexanone
         :coral,          # Methoxycyclohexane [not in Ni]
         :goldenrod1,     # Cyclohexanone
         :palegreen3,     # Cyclohexanol
         :mediumorchid1,  # Cyclohexane [Not in Ni]
         :black]          # Fouling (Placeholder)

    # Create plots of all the experimental data as solid lines, with error bars
    plot!(plt, x[1:stop,10], x[1:stop,11], linestyle=:solid, color=c[1], yerr=x[1:stop,20], markerstrokecolor=c[1], legend=false)
    plot!(plt, x[1:stop,10], x[1:stop,12], linestyle=:solid, color=c[2], yerr=x[1:stop,21], markerstrokecolor=c[2], legend=false)
    if i != 4
        plot!(plt, x[1:stop,10], x[1:stop,13], linestyle=:solid, color=c[3], yerr=x[1:stop,22], markerstrokecolor=c[3], legend=false)
    end
    plot!(plt, x[1:stop,10], x[1:stop,15], linestyle=:solid, color=c[4], yerr=x[1:stop,24], markerstrokecolor=c[4], legend=false)
    plot!(plt, x[1:stop,10], x[1:stop,16], linestyle=:solid, color=c[5], yerr=x[1:stop,25], markerstrokecolor=c[5], legend=false)
    if i != 4
        plot!(plt, x[1:stop,10], x[1:stop,17], linestyle=:solid, color=c[6], yerr=x[1:stop,26], markerstrokecolor=c[6], legend=false)
    end
    plot!(plt, x[1:stop,10], x[1:stop,19], linestyle=:solid, color=c[7], yerr=x[1:stop,27], markerstrokecolor=c[7], legend=false)

    # Create plots of all the model data as dashed lines
    plot!(plt, x[:,1], x[:,2], linestyle=:dash, color=c[1], legend=false)
    plot!(plt, x[:,1], x[:,3], linestyle=:dash, color=c[2], legend=false)
    if i != 4
        plot!(plt, x[:,1], x[:,4], linestyle=:dash, color=c[3], legend=false)
    end
    plot!(plt, x[:,1], x[:,6], linestyle=:dash, color=c[4], legend=false)
    plot!(plt, x[:,1], x[:,7], linestyle=:dash, color=c[5], legend=false)
    if i != 4
        plot!(plt, x[:,1], x[:,8], linestyle=:dash, color=c[6], legend=false)
    end
    plot!(plt, x[:,1], x[:,9], linestyle=:dash, color=c[7], legend=false)

    # Create labels
    plot!(plt, xlims=(0, 125), ylims=(-0.04, 0.24), xticks=0:20:120, xtickfontsize=10, ytickfontsize=10)

    if i==1
        title!("Ru-CAC")
    elseif i==2
        title!("Ru-FWAC-meso")
    elseif i==3
        title!("Ru-FWAC-micro")
    elseif i==4
        title!("Ni-FWAC-meso")
    end
end

# Define the layout to make a nicer-looking 4-in-1 plot
my_layout = @layout [ a{0.0001w} [grid(2,2); b{0.0001h}] c{0.3w} ]

# Create a legend "plot" for the right-most place
legend_plot = plot((1:7)', labels=["Guaiacol" "Methoxycyclohexanone" "Methoxycyclohexane" "Cyclohexanone" "Cyclohexanol" "Cyclohexane" "Placeholder"], 
                    color=[:cornflowerblue :crimson :coral :goldenrod1 :palegreen3 :mediumorchid1 :black], legend=:right, framestyle=:none, legendfontsize=15, background_color_legend=nothing)

# Create separate "plots" for the xlabel and ylabel
ylabel_plot = plot(framestyle=:none, ylabel="Concentration (mol/L)", ylabelfontsize=15)
xlabel_plot = plot(framestyle=:none, xlabel="Time (min)", xlabelfontsize=15)

# Combine all of the above into a single, big plot, and save it
decane_plot = plot(ylabel_plot, p1, p2, p3, p4, xlabel_plot, legend_plot, layout=my_layout, size=(1200,750))

# Uncomment to save:
# savefig(joinpath(@__DIR__, "decane_plot.svg"))


# Repeat the above for water
p1 = plot()
p2 = plot()
p3 = plot()
p4 = plot()
for i = 1:4
    if i==1
        plt = p1
        x = Ru_CAC_w
        stop = 8
    elseif i==2
        plt = p2
        x = Ru_FWAC_meso_w
        stop = 10
    elseif i==3
        plt = p3
        x = Ru_FWAC_micro_w
        stop = 10
    elseif i==4
        plt = p4
        x = Ni_CAC_w
        stop = 9
    end

    # Define the colors
    c = [:cornflowerblue, # Guaiacol
         :crimson,        # Methoxycyclohexanone
         :teal,           # Phenol
         :goldenrod1,     # Cyclohexanone
         :palegreen3,     # Cyclohexanol
         :black]          # Fouling (Placeholder)

    # Create plots of all the experimental data as solid lines, with error bars
    plot!(plt, x[1:stop,10], x[1:stop,11], linestyle=:solid, color=c[1], yerr=x[1:stop,20], markerstrokecolor=c[1], legend=false)
    plot!(plt, x[1:stop,10], x[1:stop,12], linestyle=:solid, color=c[2], yerr=x[1:stop,21], markerstrokecolor=c[2], legend=false)
    plot!(plt, x[1:stop,10], x[1:stop,14], linestyle=:solid, color=c[3], yerr=x[1:stop,23], markerstrokecolor=c[3], legend=false)
    plot!(plt, x[1:stop,10], x[1:stop,15], linestyle=:solid, color=c[4], yerr=x[1:stop,24], markerstrokecolor=c[4], legend=false)
    plot!(plt, x[1:stop,10], x[1:stop,16], linestyle=:solid, color=c[5], yerr=x[1:stop,25], markerstrokecolor=c[5], legend=false)
    plot!(plt, x[1:stop,10], x[1:stop,19], linestyle=:solid, color=c[6], yerr=x[1:stop,27], markerstrokecolor=c[6], legend=false)

    # Create plots of all the model data as dashed lines
    plot!(plt, x[:,1], x[:,2], linestyle=:dash, color=c[1], legend=false)
    plot!(plt, x[:,1], x[:,3], linestyle=:dash, color=c[2], legend=false)
    plot!(plt, x[:,1], x[:,5], linestyle=:dash, color=c[3], legend=false)
    plot!(plt, x[:,1], x[:,6], linestyle=:dash, color=c[4], legend=false)
    plot!(plt, x[:,1], x[:,7], linestyle=:dash, color=c[5], legend=false)
    plot!(plt, x[:,1], x[:,9], linestyle=:dash, color=c[6], legend=false)

    # Create labels and such
    if i==1
        plot!(plt, xlims=(0,61), ylims=(-0.01, 0.15), xticks=0:10:60, xtickfontsize=10, ytickfontsize=10, title="Ru-CAC")
    elseif i==2
        plot!(plt, xlims=(0,31), ylims=(-0.01, 0.15), xticks=0:5:30, xtickfontsize=10, ytickfontsize=10, title="Ru-FWAC-meso")
    elseif i==3
        plot!(plt, xlims=(0,41), ylims=(-0.01, 0.15), xticks=0:5:40, xtickfontsize=10, ytickfontsize=10, title="Ru-FWAC-micro")
    elseif i==4
        plot!(plt, xlims=(0,125), ylims=(-0.01, 0.15), xticks=0:20:120, xtickfontsize=10, ytickfontsize=10, title="Ni-CAC")
    end
end

# Define the layout to make a nicer-looking 4-in-1 plot
my_layout = @layout [ a{0.0001w} [grid(2,2); b{0.0001h}] c{0.3w} ]

# Create a legend "plot" for the right-most place
legend_plot = plot((1:6)', labels=["Guaiacol" "Methoxycyclohexanone" "Phenol" "Cyclohexanone" "Cyclohexanol" "Placeholder"], 
                    color=[:cornflowerblue :crimson :teal :goldenrod1 :palegreen3 :black], legend=:right, framestyle=:none, legendfontsize=15, background_color_legend=nothing)

#Old color: mediumorchid1, wheat4
# Create separate "plots" for the xlabel and ylabel
ylabel_plot = plot(framestyle=:none, ylabel="Concentration (mol/L)", ylabelfontsize=15)
xlabel_plot = plot(framestyle=:none, xlabel="Time (min)", xlabelfontsize=15)

# Combine all of the above into a single, big plot, and save it
water_plot = plot(ylabel_plot, p1, p2, p3, p4, xlabel_plot, legend_plot, layout=my_layout, size=(1200,750))
# savefig(joinpath(@__DIR__, "water_plot.svg"))