#  Python Code for heat effects problems. It can calculate the adiabatic flame temperature, the heat requirement, and the entropy change.
#  All formulas are taken from pages 139,183,184 of the INTRODUCTION TO CHEMICAL ENGINEERING THERMODYNAMICS, EIGHTH EDITION.
#  Parameters for Cp are from Table C.1 of the INTRODUCTION TO CHEMICAL ENGINEERING THERMODYNAMICS, EIGHTH EDITION.

import math  # "math" library imported in order to use the pi constant and the natural log(ln) function.
import sys  # "sys" imported since the program exits if the flow is in transition region.
import scipy.constants as sc  # "math" library imported in order to use the R constant (universal gas constant) with units of J/mol.K.

import pandas as pd
data=pd.read_excel('pythontable.xlsx',header=None,index_col=0,skiprows=1)
cpvalues=data.values

DH298_ethane = -83820
DH298_methane = -74520
DH298_CO2 = -393509
DH298_H2O = -241818

def mean_heat_capacity_H(A, B, C, D, T, T0, R):
    
    CpH = float(R * ( A + B*(T+T0)/2 + C*(T**2+T0**2+T*T0)/3 + D/(T*T0) ))
    
    return CpH

def standard_enthalpy_of_f(n, CO2, H2O, a):
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

def mean_heat_capacity_S(A, B, C, D, T, T0, R):
    
    CpS = float(R * ( A + ( B + ((T + T0) / 2) * (C + D / (T**2 * T0**2)) )  * ((T-T0) / math.log((T/T0)) )))

    return CpS

def standard_entropy(n, CO2, H2O, a, O2):
    
    DS298 = float(n * (CO2 + H2O - a - O2))
    
    return DS298

def total_enthalpy_change(DHair, DHa, DH298, DHp):
    # this function calculates the heat requirement.
    
    DH = DHair + DHa + DH298 + DHp
    
    return DH



calculation_type = input("Please enter the calculation type ( the letter 'T' for adiabatic flame temperature / the letter 'H' for heat requirement ) : ")


if calculation_type == "T" or calculation_type == "H":
    
    burned_alkane = input("Please enter the alkane molecule you want to burn ( the letter 'E' for ethane / the letter 'M' for methane ) : ")
    
    if burned_alkane == "E" or burned_alkane == "M":
        
        mole_of_alkane = float(input("Please enter the mole of the burned alkane: "))
        
        if mole_of_alkane > 0:
            
            excess_air = float(input("Please enter the percentage amount of excess air: "))
            
            if 100 >= excess_air >= 0:
                
                if calculation_type == "T":
        
                    T0a = float(input("Please enter the initial temperature of the burned alkane: "))

                    if T0a >= 0:
                        
                        T0air = float(input("Please enter the initial temperature of the air: "))
                        
                        if T0air >= 0:
                    
                            if burned_alkane == "E":
                                
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
                                
                                

                                
                            elif burned_alkane == "M":
                    
                                print(burned_alkane)
                                print(T0air)
                                print(T0a)
                                print(excess_air)
                                print(mole_of_alkane)
                                print(calculation_type)
                        
                        else:
                            sys.exit('You entered a negative value for the initial temperature of the air. Absolute zero (zero(0) Kelvin) is the lowest limit of the thermodynamic temperature scale. Please reboot the program and try again.')
          
                    else:
                        sys.exit('You entered a negative value for the initial temperature of the substance. Absolute zero (zero(0) Kelvin) is the lowest limit of the thermodynamic temperature scale. Please reboot the program and try again.')
                                                          
                                                                             
                elif calculation_type == "H":
                    
                    T0a = float(input("Please enter the initial temperature of the burned alkane(below 2000 K): "))

                    if 2000 > T0a >= 0:
                        
                        T0air = float(input("Please enter the initial temperature of the air(below 2000 K): "))
                        
                        if 2000 > T0air >= 0:
                    
                            if burned_alkane == "E":
                                
                                print(burned_alkane)
                                print(T0air)
                                print(T0a)
                                print(excess_air)
                                print(mole_of_alkane)
                                print(calculation_type)
                                
                            elif burned_alkane == "M":
                                
                                print(burned_alkane)
                                print(T0air)
                                print(T0a)
                                print(excess_air)
                                print(mole_of_alkane)
                                print(calculation_type)
                        
                        else:
                            sys.exit('You entered a value, that is greater than or equal to 2000 Kelvin, or a negative value, for the initial temperature of the substance. The initial temperature must be below 2000 K since the final temperature will be 2000 K. Also, absolute zero (zero(0) Kelvin) is the lowest limit of the thermodynamic temperature scale. Please reboot the program and try again.')          
                   
                    else:
                        sys.exit('You entered a value, that is greater than or equal to 2000 Kelvin, or a negative value, for the initial temperature of the substance. The initial temperature must be below 2000 K since the final temperature will be 2000 K. Also, absolute zero (zero(0) Kelvin) is the lowest limit of the thermodynamic temperature scale. Please reboot the program and try again.')          
                
            else:
                sys.exit('The percentage amount of excess air must be between zero(0) and one hundred in order to be a realistic value. Please reboot the program and try again.')
            
        else:
            sys.exit('The entered mole of burned alkane must be greater than zero(0) in order to be a realistic value. Please reboot the program and try again.')
            
    else:
        sys.exit('The entered burned substance must be "E" (ethane) or "M" (methane). Please reboot the program and try again.')

else:
    sys.exit('The entered calculation type must be "T" (adiabatic flame temperature) or "H" (heat requirement). Please reboot the program and try again.')