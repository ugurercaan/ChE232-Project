#  Python Code for heat effects problems. It can calculate the adiabatic flame temperature, the heat requirement, and the entropy change.
#  All formulas are taken from pages 139,183,184 of the INTRODUCTION TO CHEMICAL ENGINEERING THERMODYNAMICS, EIGHTH EDITION.
#  Parameters for Cp are from Table C.1 of the INTRODUCTION TO CHEMICAL ENGINEERING THERMODYNAMICS, EIGHTH EDITION.

import math  # "math" library imported in order to use the pi constant and the natural log(ln) function.
import sys  # "sys" imported since the program exits if the flow is in transition region.
import scipy.constants as sc  # "math" library imported in order to use the R constant (universal gas constant) with units of J/mol.K.

import pandas as pd
data=pd.read_excel('pythontable.xlsx',header=None,index_col=0,skiprows=1)
properties=data.values

T_C = 25
T = T_C + 273.15

def MCPH(A, B, C, D, T, T0):
    
    CpH = float(sc.R * ( A + B*(T+T0)/2 + C*(T**2+T0**2+T*T0)/3 + D/(T*T0) ))
    
    return CpH

def standard_enthalpy_of_formation(n, CO2, H2O, a):
    # the 'a' parameter is reffering the alkane in this chemical equation.
    # n is the mole number of the alkane burned.
    # all CO2, H2O, and a parameters mean its stoichiometric coefficient times the standard ethalpy of formation of that substance in the gas form.
    # we read the standard ethalpy of formation values from Table C.4 in the refferring textbook above.
    
    DH298 = float(n * ( CO2 + H2O - a))
    
    return DH298

def enthalpy_change_of_the_products_adiabatic(DHair, DHa, DH298):
    
    # the parameters are the enthalpy change of the air, the alkane burned, and the change of standard enthalpy of formation, respectively.
    
    DHp = float(-DHair - DHa - DH298)  # this equation applies due to the total enthalpy change equals zero(0) for an adiabatic system.
    
    return DHp

def Tadiabatic(DHp, Cp, T0):
    
    # in order to calculate the adiabatic flame temperature, we iterate the calculations for Cp, until convergence with an initial guess (T=2*T0).

    Tad = float(DHp / Cp + T0)
    
    return Tad

def total_entropy_change( DSair, DSa, DS298, DSp):
    
    DS = float(DSair + DSa + DS298 + DSp)
    
    return DS

def mean_heat_capacity_S(A, B, C, D, T, T0):
    
    CpS = float(sc.R * ( A + ( B + ((T + T0) / 2) * (C + D / (T**2 * T0**2)) )  * ((T-T0) / math.log((T/T0)) )))

    return CpS

def standard_entropy(n, CO2, H2O, a, O2):
    
    DS298 = float(n * (CO2 + H2O - a - O2))
    
    return DS298

def total_enthalpy_change(DHair, DHa, DH298, DHp):
    # this function calculates the heat requirement.
    
    DH = DHair + DHa + DH298 + DHp
    
    return DH

def error_func(old, new):
    # This function calculates the error which is obtained by iteration.

    error = float(abs((old - new) / old) * 100)

    return error



calculation_type = input("\nPlease enter the calculation type ( the letter 'T' for adiabatic flame temperature / the letter 'H' for heat requirement ) : ")


if calculation_type == "T" or calculation_type == "H":
    
    burned_alkane = input("\nPlease enter the alkane molecule you want to burn ( the letter 'E' for ethane / the letter 'M' for methane ) : ")
    
    if burned_alkane == "E" or burned_alkane == "M":
        
        mole_of_alkane = float(input("\nPlease enter the mole of the burned alkane: "))
        
        if mole_of_alkane > 0:
            
            excess_air = float(input("\nPlease enter the percentage amount of excess air(%): "))
            
            if 100 >= excess_air >= 0:
                
                if calculation_type == "T":
        
                    T0_a_C = float(input("\nPlease enter the initial temperature of the burned alkane(°C): "))
                    T0_a_K = T0_a_C + 273.15
                    
                    if T0_a_K >= 0:
                        
                        T0_air_C = float(input("\nPlease enter the initial temperature of the air(°C): "))
                        T0_air_K = T0_air_C + 273.15
                        
                        if T0_air_K >= 0:
                            
                            
                            if burned_alkane == "E":
                                
                                T_ethane = T0_a_K
                                
                                
                                                                
                                coef_ethane = 1
                                coef_O2 = 7/2
                                coef_CO2 = 2
                                coef_H2O = 3                                
                                
                                feed_mol_ethane = mole_of_alkane * coef_ethane
                                feed_mol_O2 = mole_of_alkane * coef_O2 * (excess_air/100 + 1)
                                feed_mol_N2 = feed_mol_O2 * (79/21)
                                
                                product_mol_O2 = feed_mol_O2 - (mole_of_alkane * coef_O2)
                                product_mol_N2 = feed_mol_N2
                                product_mol_CO2 = mole_of_alkane * coef_CO2
                                product_mol_H2O = mole_of_alkane * coef_H2O
                                
                                A_air = feed_mol_O2 * properties[1][0] + feed_mol_N2 * properties[2][0]
                                B_air = (feed_mol_O2 * properties[1][1] + feed_mol_N2 * properties[2][1]) * 10**-3
                                C_air = (feed_mol_O2 * properties[1][2] + feed_mol_N2 * properties[2][2]) * 10**-6
                                D_air = (feed_mol_O2 * properties[1][3] + feed_mol_N2 * properties[2][3]) * 10**5
                            
                                if T0_air_K == T:
                                    
                                    DH_air = 0
                                
                                else:
                                    DH_air = MCPH(A_air, B_air, C_air, D_air, T, T0_air_K) * (T-T0_air_K)
                                
                                                               
                                A_ethane = feed_mol_ethane * properties[0][0]
                                B_ethane = (feed_mol_ethane * properties[0][1]) * 10**-3
                                C_ethane = (feed_mol_ethane * properties[0][2]) * 10**-6
                                D_ethane = (feed_mol_ethane * properties[0][3]) * 10**5

                                if T_ethane == T:
                                    
                                    DH_ethane = 0
                                
                                else:
                                    DH_ethane = MCPH(A_ethane, B_ethane, C_ethane, D_ethane, T, T_ethane) * (T-T_ethane)
                                
                                DH298 = standard_enthalpy_of_formation(mole_of_alkane, coef_CO2*properties[3][4], coef_H2O*properties[4][4], coef_ethane*properties[0][4])
                                
                                DHp = enthalpy_change_of_the_products_adiabatic(DH_air, DH_ethane, DH298)
                                                                
                                A_p = product_mol_O2 * properties[1][0] + product_mol_N2 * properties[2][0] + product_mol_CO2 * properties[3][0] + product_mol_H2O * properties[4][0]
                                B_p = (product_mol_O2 * properties[1][1] + product_mol_N2 * properties[2][1] + product_mol_CO2 * properties[3][1] + product_mol_H2O * properties[4][1]) * 10**-3
                                C_p = (product_mol_O2 * properties[1][2] + product_mol_N2 * properties[2][2] + product_mol_CO2 * properties[3][2] + product_mol_H2O * properties[4][2]) * 10**-6
                                D_p = (product_mol_O2 * properties[1][3] + product_mol_N2 * properties[2][3] + product_mol_CO2 * properties[3][3] + product_mol_H2O * properties[4][3]) * 10**5
                                

                                tolerance = 10**-4
                                T_guess = 2*T
                                error = 100
                                

                                while error >= tolerance:
  
                                    CpH_product = MCPH(A_p, B_p, C_p, D_p, T_guess, T)
                                    T_ad = Tadiabatic(DHp, CpH_product, T)
                                    error = error_func(T_guess, T_ad)
                                    
                                    T_guess = T_ad
                                
                                print("\nThe adiabatic flame temperature that can be reached for", mole_of_alkane, "moles of ethane at", T0_a_C, "°C burned with", excess_air, "% excess air at", T0_air_C, "°C is", T_ad, "K.")

                                if T0_air_K == T:
                                    
                                    DS_air = 0
                                
                                else:
                                    DS_air = mean_heat_capacity_S(A_air, B_air, C_air, D_air, T, T0_air_K) * math.log(T/T0_air_K)
                                    
                                print("\nThe entropy change of air from ", T0_air_C," °C to 25°C is ", DS_air,"J/K.")
                                
                                if T_ethane == T:
                                    
                                    DS_ethane = 0
                                
                                else:
                                    DS_ethane = mean_heat_capacity_S(A_ethane, B_ethane, C_ethane, D_ethane, T, T0_a_K) * math.log(T/T0_a_K)
                                
                                print("\nThe entropy change of ethane from ", T0_a_C," °C to 25°C is ", DS_ethane,"J/K.")
                                   
                                DS298 = standard_entropy(mole_of_alkane, coef_CO2*properties[3][5], coef_H2O*properties[4][5], coef_ethane*properties[0][5], coef_O2*properties[1][5])
                                
                                print("\nThe standard entropy change while the ethane is burning is ", DS298," J/K.")
                                
                                DS_product = mean_heat_capacity_S(A_p, B_p, C_p, D_p, T_ad, T) * math.log(T_ad/T)
                                
                                print("\nThe entropy change of the products from ", T_C," °C to ", T_ad," is ", DS_ethane,"J/K.")
                                
                                DS = total_entropy_change(DS_air, DS_ethane, DS298, DS_product)
                                
                                print("\nThe entropy change that can be reached for", mole_of_alkane, "moles of ethane at", T0_a_C, "°C burned with", excess_air, "% excess air at", T0_air_C, "°C is", DS, "J/K.")
                                
                            elif burned_alkane == "M":
                    
                                T_methane = T0_a_K
                                
                                coef_methane = 1
                                coef_O2 = 2
                                coef_CO2 = 1
                                coef_H2O = 2
                                
                                feed_mol_methane = mole_of_alkane * coef_methane
                                feed_mol_O2 = mole_of_alkane * coef_O2 * (excess_air/100 + 1)
                                feed_mol_N2 = feed_mol_O2 * (79/21)
                                
                                product_mol_O2 = feed_mol_O2 - (mole_of_alkane * coef_O2)
                                product_mol_N2 = feed_mol_N2
                                product_mol_CO2 = mole_of_alkane * coef_CO2
                                product_mol_H2O = mole_of_alkane * coef_H2O
                                
                                A_air = feed_mol_O2 * properties[1][0] + feed_mol_N2 * properties[2][0]
                                B_air = (feed_mol_O2 * properties[1][1] + feed_mol_N2 * properties[2][1]) * 10**-3
                                C_air = (feed_mol_O2 * properties[1][2] + feed_mol_N2 * properties[2][2]) * 10**-6
                                D_air = (feed_mol_O2 * properties[1][3] + feed_mol_N2 * properties[2][3]) * 10**5
                                
                                if T0_air_K == T:
                                    
                                    DH_air = 0
                                
                                else:
                                    DH_air = MCPH(A_air, B_air, C_air, D_air, T, T0_air_K) * (T-T0_air_K)
                                                               
                                A_methane = feed_mol_methane * properties[5][0]
                                B_methane = (feed_mol_methane * properties[5][1]) * 10**-3
                                C_methane = (feed_mol_methane * properties[5][2]) * 10**-6
                                D_methane = (feed_mol_methane * properties[5][3]) * 10**5
                                
                                if T_methane == T:
                                    
                                    DH_methane = 0
                                
                                else:
                                    DH_methane = MCPH(A_methane, B_methane, C_methane, D_methane, T, T_methane) * (T-T_methane)

                                
                                
                                DH298 = standard_enthalpy_of_formation(mole_of_alkane, coef_CO2*properties[3][4], coef_H2O*properties[4][4], coef_methane*properties[5][4])
                                
                                DHp = enthalpy_change_of_the_products_adiabatic(DH_air, DH_methane, DH298)
                                                                
                                A_p = product_mol_O2 * properties[1][0] + product_mol_N2 * properties[2][0] + product_mol_CO2 * properties[3][0] + product_mol_H2O * properties[4][0]
                                B_p = (product_mol_O2 * properties[1][1] + product_mol_N2 * properties[2][1] + product_mol_CO2 * properties[3][1] + product_mol_H2O * properties[4][1]) * 10**-3
                                C_p = (product_mol_O2 * properties[1][2] + product_mol_N2 * properties[2][2] + product_mol_CO2 * properties[3][2] + product_mol_H2O * properties[4][2]) * 10**-6
                                D_p = (product_mol_O2 * properties[1][3] + product_mol_N2 * properties[2][3] + product_mol_CO2 * properties[3][3] + product_mol_H2O * properties[4][3]) * 10**5
                                

                                tolerance = 10**-4
                                T_guess = 2*T
                                error = 100
                                

                                while error >= tolerance:
  
                                    CpH_product = MCPH(A_p, B_p, C_p, D_p, T_guess, T)
                                    T_ad = Tadiabatic(DHp, CpH_product, T)
                                    error = error_func(T_guess, T_ad)
                                    
                                    T_guess = T_ad
                                
                                print("\nThe adiabatic flame temperature that can be reached for", mole_of_alkane, "moles of methane at", T0_a_C, "°C burned with", excess_air, "% excess air at", T0_air_C, "°C is", T_ad, "K.")

                                if T0_air_K == T:
                                    
                                    DS_air = 0
                                
                                else:
                                    DS_air = mean_heat_capacity_S(A_air, B_air, C_air, D_air, T, T0_air_K) * math.log(T/T0_air_K)
                                    
                                print("\nThe entropy change of air from ", T0_air_C," °C to 25°C is ", DS_air,"J/K.")
                                    
                                if T_methane == T:
                                    
                                    DS_methane = 0
                                
                                else:
                                    DS_methane = mean_heat_capacity_S(A_methane, B_methane, C_methane, D_methane, T, T0_a_K) * math.log(T/T0_a_K)
                                    
                                print("\nThe entropy change of methane from ", T0_a_C," °C to 25°C is ", DS_methane,"J/K.")

                                DS298 = standard_entropy(mole_of_alkane, coef_CO2*properties[3][5], coef_H2O*properties[4][5], coef_methane*properties[5][5], coef_O2*properties[1][5])
                                
                                print("\nThe standard entropy change while the ethane is burning is ", DS298," J/K.")
                                
                                DS_product = mean_heat_capacity_S(A_p, B_p, C_p, D_p, T_ad, T) * math.log(T_ad/T)
                                
                                print("\nThe entropy change of the products from ", T_C," °C to ", T_ad," is ", DS_methane,"J/K.")
                                
                                DS = total_entropy_change(DS_air, DS_methane, DS298, DS_product)
                                
                                print("\nThe entropy change that can be reached for", mole_of_alkane, "moles of methane at", T0_a_C, "°C burned with", excess_air, "% excess air at", T0_air_C, "°C is", DS, "J/K.")  
                            
                            
                            
                        else:
                            sys.exit('\nYou entered a negative value for the initial temperature of the air. Absolute zero (zero(0) Kelvin) is the lowest limit of the thermodynamic temperature scale. Please reboot the program and try again.')
          
                    else:
                        sys.exit('\nYou entered a negative value for the initial temperature of the substance. Absolute zero (zero(0) Kelvin) is the lowest limit of the thermodynamic temperature scale. Please reboot the program and try again.')
                                                          
                                                                             
                elif calculation_type == "H":
                    
                    T0_a_K = float(input("\nPlease enter the initial temperature of the burned alkane(below 2000 K)(K): "))
                    
                    if 2000 > T0_a_K >= 0:
                        
                        T0_air_K = float(input("\nPlease enter the initial temperature of the air(below 2000 K)(K): "))
                        
                        if 2000 > T0_air_K >= 0:
                    
                            if burned_alkane == "E":
                                
                                T_ethane = T0_a_K
                                
                                
                                                                
                                coef_ethane = 1
                                coef_O2 = 7/2
                                coef_CO2 = 2
                                coef_H2O = 3                                
                                
                                feed_mol_ethane = mole_of_alkane * coef_ethane
                                feed_mol_O2 = mole_of_alkane * coef_O2 * (excess_air/100 + 1)
                                feed_mol_N2 = feed_mol_O2 * (79/21)
                                
                                product_mol_O2 = feed_mol_O2 - (mole_of_alkane * coef_O2)
                                product_mol_N2 = feed_mol_N2
                                product_mol_CO2 = mole_of_alkane * coef_CO2
                                product_mol_H2O = mole_of_alkane * coef_H2O
                                
                                A_air = feed_mol_O2 * properties[1][0] + feed_mol_N2 * properties[2][0]
                                B_air = (feed_mol_O2 * properties[1][1] + feed_mol_N2 * properties[2][1]) * 10**-3
                                C_air = (feed_mol_O2 * properties[1][2] + feed_mol_N2 * properties[2][2]) * 10**-6
                                D_air = (feed_mol_O2 * properties[1][3] + feed_mol_N2 * properties[2][3]) * 10**5
                            
                                if T0_air_K == T:
                                    
                                    DH_air = 0
                                
                                else:
                                    DH_air = MCPH(A_air, B_air, C_air, D_air, T, T0_air_K) * (T-T0_air_K)
                                
                                                               
                                A_ethane = feed_mol_ethane * properties[0][0]
                                B_ethane = (feed_mol_ethane * properties[0][1]) * 10**-3
                                C_ethane = (feed_mol_ethane * properties[0][2]) * 10**-6
                                D_ethane = (feed_mol_ethane * properties[0][3]) * 10**5

                                if T_ethane == T:
                                    
                                    DH_ethane = 0
                                
                                else:
                                    DH_ethane = MCPH(A_ethane, B_ethane, C_ethane, D_ethane, T, T_ethane) * (T-T_ethane)
                                
                                DH298 = standard_enthalpy_of_formation(mole_of_alkane, coef_CO2*properties[3][4], coef_H2O*properties[4][4], coef_ethane*properties[0][4])
                                
                                
                                                                
                                A_p = product_mol_O2 * properties[1][0] + product_mol_N2 * properties[2][0] + product_mol_CO2 * properties[3][0] + product_mol_H2O * properties[4][0]
                                B_p = (product_mol_O2 * properties[1][1] + product_mol_N2 * properties[2][1] + product_mol_CO2 * properties[3][1] + product_mol_H2O * properties[4][1]) * 10**-3
                                C_p = (product_mol_O2 * properties[1][2] + product_mol_N2 * properties[2][2] + product_mol_CO2 * properties[3][2] + product_mol_H2O * properties[4][2]) * 10**-6
                                D_p = (product_mol_O2 * properties[1][3] + product_mol_N2 * properties[2][3] + product_mol_CO2 * properties[3][3] + product_mol_H2O * properties[4][3]) * 10**5
                                
                                DH_p = MCPH(A_p, B_p, C_p, D_p, 2000, T) * (2000 - T)
                                
                                DH = total_enthalpy_change(DH_air, DH_ethane, DH298, DH_p)
                                
                                print("\nThe heat requirement for the products to reach a temperature of 2000 K for", mole_of_alkane,"moles of ethane at ", T0_a_K,"K burned with ", excess_air,"% excess air at ", T0_air_K," K is ", DH,"J.")
                                
                                if T0_air_K == T:
                                    
                                    DS_air = 0
                                
                                else:
                                    DS_air = mean_heat_capacity_S(A_air, B_air, C_air, D_air, T, T0_air_K) * math.log(T/T0_air_K)
                                    
                                print("\nThe entropy change of air from ", T0_air_K," K to 298.15 K is ", DS_air,"J/K.")
                                
                                if T_ethane == T:
                                    
                                    DS_ethane = 0
                                
                                else:
                                    DS_ethane = mean_heat_capacity_S(A_ethane, B_ethane, C_ethane, D_ethane, T, T0_a_K) * math.log(T/T0_a_K)
                                
                                print("\nThe entropy change of ethane from ", T0_a_K," K to 298.15 K is ", DS_ethane,"J/K.")
                                   
                                DS298 = standard_entropy(mole_of_alkane, coef_CO2*properties[3][5], coef_H2O*properties[4][5], coef_ethane*properties[0][5], coef_O2*properties[1][5])
                                
                                print("\nThe standard entropy change while the ethane is burning is ", DS298," J/K.")
                                
                                DS_product = mean_heat_capacity_S(A_p, B_p, C_p, D_p, 2000, T) * math.log(2000/T)
                                
                                print("\nThe entropy change of the products from ", T," K to ", 2000," is ", DS_ethane,"J/K.")
                                
                                DS = total_entropy_change(DS_air, DS_ethane, DS298, DS_product)
                                
                                print("\nThe entropy change that can be reached for", mole_of_alkane, "moles of ethane at", T0_a_K, "K burned with", excess_air, "% excess air at", T0_air_K, "K is", DS, "J/K.")

                                
                            elif burned_alkane == "M":
                                
                                T_methane = T0_a_K
                                
                                coef_methane = 1
                                coef_O2 = 2
                                coef_CO2 = 1
                                coef_H2O = 2
                                
                                feed_mol_methane = mole_of_alkane * coef_methane
                                feed_mol_O2 = mole_of_alkane * coef_O2 * (excess_air/100 + 1)
                                feed_mol_N2 = feed_mol_O2 * (79/21)
                                
                                product_mol_O2 = feed_mol_O2 - (mole_of_alkane * coef_O2)
                                product_mol_N2 = feed_mol_N2
                                product_mol_CO2 = mole_of_alkane * coef_CO2
                                product_mol_H2O = mole_of_alkane * coef_H2O
                                
                                A_air = feed_mol_O2 * properties[1][0] + feed_mol_N2 * properties[2][0]
                                B_air = (feed_mol_O2 * properties[1][1] + feed_mol_N2 * properties[2][1]) * 10**-3
                                C_air = (feed_mol_O2 * properties[1][2] + feed_mol_N2 * properties[2][2]) * 10**-6
                                D_air = (feed_mol_O2 * properties[1][3] + feed_mol_N2 * properties[2][3]) * 10**5
                                
                                if T0_air_K == T:
                                    
                                    DH_air = 0
                                
                                else:
                                    DH_air = MCPH(A_air, B_air, C_air, D_air, T, T0_air_K) * (T-T0_air_K)
                                                               
                                A_methane = feed_mol_methane * properties[5][0]
                                B_methane = (feed_mol_methane * properties[5][1]) * 10**-3
                                C_methane = (feed_mol_methane * properties[5][2]) * 10**-6
                                D_methane = (feed_mol_methane * properties[5][3]) * 10**5
                                
                                if T_methane == T:
                                    
                                    DH_methane = 0
                                
                                else:
                                    DH_methane = MCPH(A_methane, B_methane, C_methane, D_methane, T, T_methane) * (T-T_methane)

                                
                                
                                DH298 = standard_enthalpy_of_formation(mole_of_alkane, coef_CO2*properties[3][4], coef_H2O*properties[4][4], coef_methane*properties[5][4])
                                
                                                                
                                A_p = product_mol_O2 * properties[1][0] + product_mol_N2 * properties[2][0] + product_mol_CO2 * properties[3][0] + product_mol_H2O * properties[4][0]
                                B_p = (product_mol_O2 * properties[1][1] + product_mol_N2 * properties[2][1] + product_mol_CO2 * properties[3][1] + product_mol_H2O * properties[4][1]) * 10**-3
                                C_p = (product_mol_O2 * properties[1][2] + product_mol_N2 * properties[2][2] + product_mol_CO2 * properties[3][2] + product_mol_H2O * properties[4][2]) * 10**-6
                                D_p = (product_mol_O2 * properties[1][3] + product_mol_N2 * properties[2][3] + product_mol_CO2 * properties[3][3] + product_mol_H2O * properties[4][3]) * 10**5
                                

                                DH_p = MCPH(A_p, B_p, C_p, D_p, 2000, T) * (2000 - T)
                                
                                DH = total_enthalpy_change(DH_air, DH_methane, DH298, DH_p)
                                
                                print("\nThe heat requirement for the products to reach a temperature of 2000 K for", mole_of_alkane,"moles of methane at ", T0_a_K," K burned with ", excess_air,"% excess air at ", T0_air_K," K is ", DH,"J.")

                                
                                if T0_air_K == T:
                                    
                                    DS_air = 0
                                
                                else:
                                    DS_air = mean_heat_capacity_S(A_air, B_air, C_air, D_air, T, T0_air_K) * math.log(T/T0_air_K)
                                    
                                print("\nThe entropy change of air from ", T0_air_K," K to 298.15 K is ", DS_air,"J/K.")
                                    
                                if T_methane == T:
                                    
                                    DS_methane = 0
                                
                                else:
                                    DS_methane = mean_heat_capacity_S(A_methane, B_methane, C_methane, D_methane, T, T0_a_K) * math.log(T/T0_a_K)
                                    
                                print("\nThe entropy change of methane from ", T0_a_K," K to 298.15 K is ", DS_methane,"J/K.")

                                DS298 = standard_entropy(mole_of_alkane, coef_CO2*properties[3][5], coef_H2O*properties[4][5], coef_methane*properties[5][5], coef_O2*properties[1][5])
                                
                                print("\nThe standard entropy change while the ethane is burning is ", DS298," J/K.")
                                
                                DS_product = mean_heat_capacity_S(A_p, B_p, C_p, D_p, 2000, T) * math.log(2000/T)
                                
                                print("\nThe entropy change of the products from ", T," K to ", 2000," is ", DS_methane,"J/K.")
                                
                                DS = total_entropy_change(DS_air, DS_methane, DS298, DS_product)
                                
                                print("\nThe entropy change that can be reached for", mole_of_alkane, "moles of methane at", T0_a_K, "K burned with", excess_air, "% excess air at", T0_air_K, "K is", DS, "J/K.")  

                        
                        else:
                            sys.exit('\nYou entered a value, that is greater than or equal to 2000 Kelvin, or a negative value, for the initial temperature of the substance. The initial temperature must be below 2000 K since the final temperature will be 2000 K. Also, absolute zero (zero(0) Kelvin) is the lowest limit of the thermodynamic temperature scale. Please reboot the program and try again.')          
                   
                    else:
                        sys.exit('\nYou entered a value, that is greater than or equal to 2000 Kelvin, or a negative value, for the initial temperature of the substance. The initial temperature must be below 2000 K since the final temperature will be 2000 K. Also, absolute zero (zero(0) Kelvin) is the lowest limit of the thermodynamic temperature scale. Please reboot the program and try again.')          
                
            else:
                sys.exit('\nThe percentage amount of excess air must be between zero(0) and one hundred in order to be a realistic value. Please reboot the program and try again.')
            
        else:
            sys.exit('\nThe entered mole of burned alkane must be greater than zero(0) in order to be a realistic value. Please reboot the program and try again.')
            
    else:
        sys.exit('\nThe entered burned substance must be "E" (ethane) or "M" (methane). Please reboot the program and try again.')

else:
    sys.exit('\nThe entered calculation type must be "T" (adiabatic flame temperature) or "H" (heat requirement). Please reboot the program and try again.')