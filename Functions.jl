
using DifferentialEquations

# Define right-hand-side functions and objective functions for each reaction. Begin
# by importing the data
include(joinpath(@__DIR__, "Data_Import.jl"))


###############################################################################
############ Reaction 1: 10%Ni-CAC in water
###############################################################################
function rhs_e1!(dx, x, p ,t) 
    T = 250+273.15 #Kelvin
    P = 750 #psi
    dx[1] = -(p[1]+p[2]+p[3]+p[4]+p[5])*x[1] + p[11]*x[6] #Guaiacol
    dx[2] = p[2]*x[1] + -(p[6]+p[7])*x[2] #Methoxycyclohexanone
    dx[3] = p[3]*x[1] + -(p[8]+p[9])*x[3] #Phenol
    dx[4] = p[4]*x[1] + p[6]*x[2] + p[8]*x[3] -(p[10])*x[4] #Cyclohexanone
    dx[5] = p[5]*x[1] + p[7]*x[2] + p[9]*x[3] + p[10]*x[4] #Cyclohexanol
    dx[6] = p[1]*x[1] - p[11]*x[6] #Coking
    nothing
end
function objective_func_e1(new_p)
    yData = e1
    x0_1 = [0.143846808902391; 0.0; 0.0; 0.0; 0.0; 0.0]
    tspan = (0.0, 120.0)
    new_prob = ODEProblem(rhs_e1!, x0_1, tspan, new_p) #Solution is in terms of the parameters p
    sol = DifferentialEquations.solve(new_prob, tstops=[10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 120.0])

    error = zero(typeof(new_p[1]))
    counter = 2

    for t in [10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 120.0]
        sol_val = sol(t)
        error += (sol_val[1] - yData.c_guaiacol[counter])^2
        error += (sol_val[2] - yData.c_methoxycyclohexanone[counter])^2
        error += (sol_val[3] - yData.c_phenol[counter])^2
        error += (sol_val[4] - yData.c_cyclohexanone[counter])^2
        error += (sol_val[5] - yData.c_cyclohexanol[counter])^2
        error += (sol_val[6] - yData.c_coke[counter])^2
        counter += 1
    end
    return error
end


###############################################################################
############ Reaction 2: 5%Ru-CAC in decane
###############################################################################
function rhs_e2!(dx, x, p ,t)
    T = 250+273.15 #Kelvin
    P = 750 #psi

    dx[1] = -(p[1]+p[2]+p[3]+p[4]+p[5]+p[6])*x[1] + p[17]*x[7] #Guaiacol
    dx[2] = p[2]*x[1] - (p[7]+p[8]+p[9]+p[10])*x[2] #Methoxycyclohexanone
    dx[3] = p[3]*x[1] + p[7]*x[2] - (p[11]+p[12]+p[13])*x[3] #Methoxycyclohexane
    dx[4] = p[4]*x[1] + p[8]*x[2] + p[11]*x[3] - (p[14]+p[15])*x[4] #Cyclohexanone
    dx[5] = p[5]*x[1] + p[9]*x[2] + p[12]*x[3] + p[14]*x[4] - (p[16])*x[5] #Cyclohexanol
    dx[6] = p[6]*x[1] + p[10]*x[2] + p[13]*x[3] + p[15]*x[4] + p[16]*x[5] #Cyclohexane
    dx[7] = p[1]*x[1] - p[17]*x[7] #Coking
    nothing
end
function objective_func_e2(new_p)
    yData = e2
    x0_2 = [0.232874906630344; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
    tspan = (0.0, 120.0)
    new_prob = ODEProblem(rhs_e2!, x0_2, tspan, new_p) #Solution is in terms of the parameters p
    sol = DifferentialEquations.solve(new_prob, tstops=[10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 120.0])

    error = zero(typeof(new_p[1]))
    counter = 2

    for t in [10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 120.0]
        sol_val = sol(t)
        error += (sol_val[1] - yData.c_guaiacol[counter])^2
        error += (sol_val[2] - yData.c_methoxycyclohexanone[counter])^2
        error += (sol_val[3] - yData.c_methoxycyclohexane[counter])^2
        error += (sol_val[4] - yData.c_cyclohexanone[counter])^2
        error += (sol_val[5] - yData.c_cyclohexanol[counter])^2
        error += (sol_val[6] - yData.c_cyclohexane[counter])^2
        error += (sol_val[7] - yData.c_coke[counter])^2
        counter += 1
    end
    return error
end


###############################################################################
############ Reaction 3: 5%Ru-CAC in water
###############################################################################
function rhs_e3!(dx, x, p ,t)
    T = 250+273.15 #Kelvin
    P = 750 #psi

    dx[1] = -(p[1]+p[2]+p[3]+p[4]+p[5])*x[1] + p[11]*x[6] #Guaiacol
    dx[2] = p[2]*x[1] - (p[6]+p[7])*x[2] #Methoxycyclohexanone
    dx[3] = p[3]*x[1] - (p[8]+p[9])*x[3] #Phenol
    dx[4] = p[4]*x[1] + p[6]*x[2] + p[8]*x[3] - (p[9])*x[4] #Cyclohexanone
    dx[5] = p[5]*x[1] + p[7]*x[2] + p[9]*x[3] + p[10]*x[4] #Cyclohexanol
    dx[6] = p[1]*x[1] - p[11]*x[6] #Coking
    nothing
end
function objective_func_e3(new_p) 
    yData = e3
    x0_3 = [0.143846808902391; 0.0; 0.0; 0.0; 0.0; 0.0]
    tspan = (0.0, 60.0)
    new_prob = ODEProblem(rhs_e3!, x0_3, tspan, new_p) #Solution is in terms of the parameters p
    sol = DifferentialEquations.solve(new_prob, tstops=[5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0])

    error = zero(typeof(new_p[1]))
    counter = 2

    for t in [5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0]
        sol_val = sol(t)
        error += (sol_val[1] - yData.c_guaiacol[counter])^2
        error += (sol_val[2] - yData.c_methoxycyclohexanone[counter])^2
        error += (sol_val[3] - yData.c_phenol[counter])^2
        error += (sol_val[4] - yData.c_cyclohexanone[counter])^2
        error += (sol_val[5] - yData.c_cyclohexanol[counter])^2
        error += (sol_val[6] - yData.c_coke[counter])^2
        counter += 1
    end
    return error
end


###############################################################################
############ Reaction 4: 10%Ni-FWAC-meso in decane
###############################################################################
function rhs_e4!(dx, x, p ,t)
    T = 250+273.15 #Kelvin
    P = 750 #psi

    dx[1] = -(p[1]+p[2]+p[3]+p[4])*x[1] + p[8]*x[5] #Guaiacol
    dx[2] = (p[2]*x[1]) + -(p[5]+p[6])*x[2] #Methoxycyclohexanone
    dx[3] = (p[3]*x[1] + p[5]*x[2]) + -(p[7])*x[3] #Cyclohexanone
    dx[4] = (p[4]*x[1] + p[6]*x[2] + p[7]*x[3]) #Cyclohexanol
    dx[5] = p[1]*x[1] - p[8]*x[5] #Coking
    nothing
end
function objective_func_e4(new_p)
    yData = e4
    x0_4 = [0.232874906630344; 0.0; 0.0; 0.0; 0.0]
    tspan = (0.0, 120.0)
    new_prob = ODEProblem(rhs_e4!, x0_4, tspan, new_p) #Solution is in terms of the parameters p
    sol = DifferentialEquations.solve(new_prob, tstops=[10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 120.0])

    error = zero(typeof(new_p[1]))
    counter = 2

    for t in [10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 120.0]
        sol_val = sol(t)
        error += (sol_val[1] - yData.c_guaiacol[counter])^2
        error += (sol_val[2] - yData.c_methoxycyclohexanone[counter])^2
        error += (sol_val[3] - yData.c_cyclohexanone[counter])^2
        error += (sol_val[4] - yData.c_cyclohexanol[counter])^2
        error += (sol_val[5] - yData.c_coke[counter])^2
        counter += 1
    end
    return error
end


###############################################################################
############ Reaction 5: 5%Ru-FWAC-meso in decane
###############################################################################
function rhs_e5!(dx, x, p ,t)
    T = 250+273.15 #Kelvin
    P = 750 #psi

    dx[1] = -(p[1]+p[2]+p[3]+p[4]+p[5]+p[6])*x[1] + p[17]*x[7] #Guaiacol
    dx[2] = (p[2]*x[1]) + -(p[7]+p[8]+p[9]+p[10])*x[2] #Methoxycyclohexanone
    dx[3] = (p[3]*x[1] + p[7]*x[2]) + -(p[11]+p[12]+p[13])*x[3] #Methoxycyclohexane
    dx[4] = (p[4]*x[1] + p[8]*x[2] + p[11]*x[3]) + -(p[14]+p[15])*x[4] #Cyclohexanone
    dx[5] = (p[5]*x[1] + p[9]*x[2] + p[12]*x[3] + p[14]*x[4]) + -(p[16])*x[5] #Cyclohexanol
    dx[6] = (p[6]*x[1] + p[10]*x[2] + p[13]*x[3] + p[15]*x[4] + p[16]*x[5]) #Cyclohexane
    dx[7] = p[1]*x[1] - p[17]*x[7] #Coking
    nothing
end
function objective_func_e5(new_p)
    yData = e5
    x0_5 = [0.232874906630344; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
    tspan = (0.0, 120.0)
    new_prob = ODEProblem(rhs_e5!, x0_5, tspan, new_p) #Solution is in terms of the parameters p
    sol = DifferentialEquations.solve(new_prob, tstops=[2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 80.0, 120.0])

    error = zero(typeof(new_p[1]))
    counter = 2

    for t in [2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 80.0, 120.0]
        sol_val = sol(t)
        error += (sol_val[1] - yData.c_guaiacol[counter])^2
        error += (sol_val[2] - yData.c_methoxycyclohexanone[counter])^2
        error += (sol_val[3] - yData.c_methoxycyclohexane[counter])^2
        error += (sol_val[4] - yData.c_cyclohexanone[counter])^2
        error += (sol_val[5] - yData.c_cyclohexanol[counter])^2
        error += (sol_val[6] - yData.c_cyclohexane[counter])^2
        error += (sol_val[7] - yData.c_coke[counter])^2
        counter += 1
    end
    return error
end


###############################################################################
############ Reaction 6: 5%Ru-FWAC-meso in water
###############################################################################
function rhs_e6!(dx, x, p ,t)
    T = 250+273.15 #Kelvin
    P = 750 #psi

    dx[1] = -(p[1]+p[2]+p[3]+p[4]+p[5])*x[1] + p[11]*x[6] #Guaiacol
    dx[2] = (p[2]*x[1]) + -(p[6]+p[7])*x[2] #Methoxycyclohexanone
    dx[3] = (p[3]*x[1]) + -(p[8]+p[9])*x[3] #Phenol
    dx[4] = (p[4]*x[1] + p[6]*x[2] + p[8]*x[3]) + -(p[10])*x[4] #Cyclohexanone
    dx[5] = (p[5]*x[1] + p[7]*x[2] + p[9]*x[3] + p[10]*x[4]) #Cyclohexanol
    dx[6] = p[1]*x[1] - p[11]*x[6] #Coking
    nothing
end
function objective_func_e6(new_p)
    yData = e6
    x0_6 = [0.143846809; 0.0; 0.0; 0.0; 0.0; 0.0]
    tspan = (0.0, 30.0)
    new_prob = ODEProblem(rhs_e6!, x0_6, tspan, new_p) #Solution is in terms of the parameters p
    sol = DifferentialEquations.solve(new_prob, tstops=[1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0])

    error = zero(typeof(new_p[1]))
    counter = 2

    for t in [1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0]
        sol_val = sol(t)
        error += (sol_val[1] - yData.c_guaiacol[counter])^2
        error += (sol_val[2] - yData.c_methoxycyclohexanone[counter])^2
        error += (sol_val[3] - yData.c_phenol[counter])^2
        error += (sol_val[4] - yData.c_cyclohexanone[counter])^2
        error += (sol_val[5] - yData.c_cyclohexanol[counter])^2
        error += (sol_val[6] - yData.c_coke[counter])^2
        counter += 1
    end
    return error
end


###############################################################################
############ Reaction 7: 5%Ru-FWAC-micro in decane
###############################################################################
function rhs_e7!(dx, x, p ,t)
    T = 250+273.15 #Kelvin
    P = 750 #psi

    dx[1] = -(p[1]+p[2]+p[3]+p[4]+p[5]+p[6])*x[1] + p[17]*x[7] #Guaiacol
    dx[2] = (p[2]*x[1]) + -(p[7]+p[8]+p[9]+p[10])*x[2] #Methoxycyclohexanone
    dx[3] = (p[3]*x[1] + p[7]*x[2]) + -(p[11]+p[12]+p[13])*x[3] #Methoxycyclohexane
    dx[4] = (p[4]*x[1] + p[8]*x[2] + p[11]*x[3]) + -(p[14]+p[15])*x[4] #Cyclohexanone
    dx[5] = (p[5]*x[1] + p[9]*x[2] + p[12]*x[3] + p[14]*x[4]) + -(p[16])*x[5] #Cyclohexanol
    dx[6] = (p[6]*x[1] + p[10]*x[2] + p[13]*x[3] + p[15]*x[4] + p[16]*x[5]) #Cyclohexane
    dx[7] = p[1]*x[1] - p[17]*x[7] #Coking
    nothing
end
function objective_func_e7(new_p)
    yData = e7
    x0_7 = [0.232874906630344; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]
    tspan = (0.0, 120.0)
    new_prob = ODEProblem(rhs_e7!, x0_7, tspan, new_p) #Solution is in terms of the parameters p
    sol = DifferentialEquations.solve(new_prob, tstops=[5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 80.0, 120.0])

    error = zero(typeof(new_p[1]))
    counter = 2

    for t in [5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 80.0, 120.0]
        sol_val = sol(t)
        error += (sol_val[1] - yData.c_guaiacol[counter])^2
        error += (sol_val[2] - yData.c_methoxycyclohexanone[counter])^2
        error += (sol_val[3] - yData.c_methoxycyclohexane[counter])^2
        error += (sol_val[4] - yData.c_cyclohexanone[counter])^2
        error += (sol_val[5] - yData.c_cyclohexanol[counter])^2
        error += (sol_val[6] - yData.c_cyclohexane[counter])^2
        error += (sol_val[7] - yData.c_coke[counter])^2
        counter += 1
    end
    return error
end


###############################################################################
############ Reaction 8: 5%Ru-FWAC-micro in water
###############################################################################
function rhs_e8!(dx, x, p ,t)
    T = 250+273.15 #Kelvin
    P = 750 #psi

    dx[1] = -(p[1]+p[2]+p[3]+p[4]+p[5])*x[1] + p[11]*x[6] #Guaiacol
    dx[2] = p[2]*x[1] + -(p[6]+p[7])*x[2] #Methoxycyclohexanone
    dx[3] = p[3]*x[1] + -(p[8]+p[9])*x[3] #Phenol
    dx[4] = p[4]*x[1] + p[6]*x[2] + p[8]*x[3] + -(p[10])*x[4] #Cyclohexanone
    dx[5] = p[5]*x[1] + p[7]*x[2] + p[9]*x[3] + p[10]*x[4] #Cyclohexanol
    dx[6] = p[1]*x[1] - p[11]*x[6] #Coking
    nothing
end
function objective_func_e8(new_p)
    yData = e8
    x0_8 = [0.143846809; 0.0; 0.0; 0.0; 0.0; 0.0]
    tspan = (0.0, 40.0)
    new_prob = ODEProblem(rhs_e8!, x0_8, tspan, new_p) #Solution is in terms of the parameters p
    sol = DifferentialEquations.solve(new_prob, tstops=[1.0, 2.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0])

    error = zero(typeof(new_p[1]))
    counter = 2 

    for t in [1.0, 2.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 40.0]
        sol_val = sol(t)
        error += (sol_val[1] - yData.c_guaiacol[counter])^2
        error += (sol_val[2] - yData.c_methoxycyclohexanone[counter])^2
        error += (sol_val[3] - yData.c_phenol[counter])^2
        error += (sol_val[4] - yData.c_cyclohexanone[counter])^2
        error += (sol_val[5] - yData.c_cyclohexanol[counter])^2
        error += (sol_val[6] - yData.c_coke[counter])^2
        counter += 1
    end
    return error
end