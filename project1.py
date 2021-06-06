#  Python Code for a piping problem which can capable to solve three different types of cases.
#  All formulas are taken from pages 144 and 145 of the Fluid Mechanics for Chemical Engineers 2nd Edition.

"""
The FPS unit system part of the program has run properly based on the values in Example 3.2-3.3 from pages 146-150
of the Fluid Mechanics for Chemical Engineers 2nd Edition.
"""

import math  # "math" library imported in order to use the pi constant and the natural log(ln) function.
import sys  # "sys" imported since the program exits if the flow is in transition region.

in2_to_ft2 = 1 / 144  # conversion factor for the use of turning in^2 into ft^2

# gravitational accelerations
g_SI = 9.8066  # m/s^2
g_FPS = 32.174  # ft/s^2

# conversion factor from lbm/(ft.s^2) to psi while handling with the FPS unit system
psi_conversion = (g_FPS / in2_to_ft2)

"""
All defined functions are 

1- velocity_mean, 
2- reynolds, 
3- friction_factor, 
4- error_func 
"""


def velocity_mean(q, d):
    # This function calculates the mean velocity for given values.

    u = float((4 * q) / (math.pi * d ** 2))

    return u


def reynolds(d, rho, u_mean, mu):
    # This function calculates the Reynolds number for given values.

    reynolds1 = float((d * rho * u_mean) / mu)

    return reynolds1


def friction_factor(eps_d, rey):
    # This function calculates the friction factor for given values.

    if rey <= 2000:  # Laminar flow

        ffactor = float(16 / rey)

    elif 2000 < rey <= 4000:  # Transition region

        sys.exit('The flow is not laminar nor turbulent. Based on the Reynolds number calculated, this is the transition region. Please reboot the program and set your values properly.')

    elif rey > 4000:  # Turbulent flow

        ffactor = float(((-1.737) * math.log((0.269 * eps_d) - (2.185 / rey) * math.log((0.269 * eps_d) + (14.5 / rey)))) ** (-2))

    return ffactor


def error_func(old, new):
    # This function calculates the error in Reynolds number which is obtained by iteration.

    error = float(abs((old - new) / old) * 100)

    return error


# Firstly, inputs about case types and unit systems are requested from user
"""
The first two if loops check whether user inputs are in correct form for case number and unit system or not.
If they are not in the correct form, the program shows a warning and exits.
"""

case_type = float(input("PLease enter the case type(1,2, or 3): "))
if case_type.is_integer() and 3 >= case_type >= 1:
    # .is_integer() function checks whether number which is requested from user as case number is an integer or not.

    unit_system = input("Please enter the unit system you want to use (SI for m/kg/s or FPS for ft/lbm/s): ")
    if unit_system == "SI" or unit_system == "FPS":

        if unit_system == "SI":
            # The below parameters are identical inputs for all the case types while using SI unit system.

            rho = float(input("\nPlease specify the density of the fluid ρ (kg/m^3): "))
            mu = float(input("Please specify the dynamic viscosity of the fluid μ (kg/(m·s): "))
            eps = float(input("Please specify the pipe surface roughness ε (m): "))
            deltaz = float(input("Please specify the elevation difference Δz (m): "))
            L = float(input("Please specify the pipe length L (m): "))

            if case_type == 1:

                D = float(input("Please specify the pipe diameter D (m): "))
                Q = float(input("Please specify the flow rate Q (m^3/s): "))

                u_mean = velocity_mean(Q, D)
                Re = reynolds(D, rho, u_mean, mu)
                eps_D = eps / D  # It is the roughness ratio ε/D
                fF = friction_factor(eps_D, Re)

                deltaP = ((2 * fF * rho * (u_mean ** 2) * (L / D)) + (rho * g_SI * deltaz)) * (-1)

                print("\nValue of the pressure drop Δp(p1 - p2) in SI is:", deltaP, "Pa")


            elif case_type == 2:

                deltaP = float(input("Please specify the pressure drop(p1 - p2)(negative) Δp (Pa): "))
                D = float(input("Please specify the pipe diameter D (m): "))

                tol = 0.001  # The convergence tolerance of the result gotten from the iteration method.
                Re_ig = 100000.0  # It is initial guess for Reynolds number
                error = 1.0  # For the initial computation in the while loop, 1 is merely a random value that satisfies error>tol.

                while error > tol:
                    eps_D = eps / D  # It is the roughness ratio ε/D
                    fF = friction_factor(eps_D, Re_ig)
                    u_mean = ((-deltaP - (rho * g_SI * deltaz)) * (D / (2 * fF * rho * L))) ** (1 / 2)
                    Re = reynolds(D, rho, u_mean, mu)
                    error = error_func(Re_ig, Re)
                    Re_ig = Re

                Q = (u_mean * math.pi * D ** 2) / 4

                print("\nValue of the flow rate Q in SI is:", Q, "m^3/s")


            elif case_type == 3:

                Q = float(input("Please specify the flow rate Q (m^3/s): "))
                deltaP = float(input("Please specify the pressure drop(p1 - p2)(negative) Δp (Pa): "))

                tol = 0.001
                D = 10  # It is initial guess for pipe diameter D (m)
                error = 1  # For the initial computation in the while loop, 1 is merely a random value that satisfies error>tol.

                while error > tol:
                    u_mean = velocity_mean(Q, D)
                    Re = reynolds(D, rho, u_mean, mu)
                    eps_D = eps / D  # It is the roughness ratio ε/D
                    fF = friction_factor(eps_D, Re)

                    D = ((32 * fF * rho * Q ** 2 * L) / (math.pi ** 2 * (-deltaP - (rho * g_SI * deltaz)))) ** (1 / 5)

                    u_mean = velocity_mean(Q, D)
                    Re_new = reynolds(D, rho, u_mean, mu)
                    error = error_func(Re, Re_new)

                    Re = Re_new

                print("\nValue of the pipe diameter D in SI is: ", D, "(m)")


        elif unit_system == "FPS":
            # The below parameters are identical inputs for all the case types while using FPS unit system.

            rho = float(input("\nPlease specify the density of the fluid ρ (lbm/ft^3): "))
            mu = float(input("Please specify the dynamic viscosity of the fluid μ (lbm/(ft·s): "))
            eps = float(input("Please specify the pipe surface roughness ε (ft): "))
            deltaz = float(input("Please specify the elevation difference Δz (ft): "))
            L = float(input("Please specify the pipe length L (ft): "))

            if case_type == 1:
                D = float(input("Please specify the pipe diameter D (ft): "))
                Q = float(input("Please specify the flow rate Q (ft^3/s): "))

                u_mean = velocity_mean(Q, D)
                Re = reynolds(D, rho, u_mean, mu)
                eps_D = eps / D  # It is the roughness ratio ε/D
                fF = friction_factor(eps_D, Re)

                deltaP = (((2 * fF * rho * (u_mean ** 2) * (L / D)) + (rho * g_FPS * deltaz)) * (-1)) / psi_conversion

                print("\nValue of the pressure drop(p1 - p2) Δp in FPS is:", deltaP, "psi")

            elif case_type == 2:

                deltaP = float(input("Please specify the pressure drop(p1 - p2)(negative) Δp (psi): "))
                D = float(input("Please specify the pipe diameter D (ft): "))

                tol = 0.001
                Re_ig = 100000.0  # It is initial guess for Reynolds number
                error = 1.0  # For the initial computation in the while loop, 1 is merely a random value that satisfies error>tol.

                while error > tol:
                    eps_D = eps / D  # It is the roughness ratio ε/D
                    fF = friction_factor(eps_D, Re_ig)
                    u_mean = ((-deltaP * psi_conversion - (rho * g_FPS * deltaz)) * (D / (2 * fF * rho * L))) ** (1 / 2)
                    Re = reynolds(D, rho, u_mean, mu)
                    error = error_func(Re_ig, Re)
                    Re_ig = Re

                Q = (u_mean * math.pi * D ** 2) / 4

                print("\nValue of the flow rate Q in FPS is:", Q, "ft^3/s")


            elif case_type == 3:

                Q = float(input("Please specify the flow rate Q (ft^3/s): "))
                deltaP = float(input("Please specify the pressure drop(p1 - p2)(negative) Δp (psi): "))

                tol = 0.001
                D = 10  # It is initial guess for pipe diameter D (m)
                error = 1  # For the initial computation in the while loop, 1 is merely a random value that satisfies error>tol.

                while error > tol:
                    u_mean = velocity_mean(Q, D)
                    Re = reynolds(D, rho, u_mean, mu)
                    eps_D = eps / D  # It is the roughness ratio ε/D
                    fF = friction_factor(eps_D, Re)

                    D = ((32 * fF * rho * Q ** 2 * L) / (math.pi ** 2 * (-deltaP * psi_conversion - (rho * g_FPS * deltaz)))) ** (
                            1 / 5)

                    u_mean = velocity_mean(Q, D)
                    Re_new = reynolds(D, rho, u_mean, mu)
                    error = error_func(Re, Re_new)

                    Re = Re_new

                print("\nValue of the pipe diameter D in FPS is: ", D, "(ft)")

    else:
        sys.exit('The unit system must be "SI" or "FPS". Please reboot the program and try again.')
else:
    sys.exit('The case type must be "1", "2" or "3". Please reboot the program and try again.')
