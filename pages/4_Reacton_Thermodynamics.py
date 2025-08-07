import sys
sys.path.append('../')

## streamlit 
import streamlit as st
import streamlit.components.v1 as components

import os 
from functions.prepare_input import *
from equilibrator_api import ComponentContribution, Q_



# the current dir is where the main app.py is 
main_dir= "./"
# output dir
output_dir=os.path.join(main_dir,"output")
data_dir = os.path.join(main_dir,"data")


# Form title
st.title("Reaction Thermodynamics")

#------------------------------------------#
# Basic ΔG’ calculations
#------------------------------------------#


# 1. condition
# optional: changing the aqueous environment parameters
# https://pmc.ncbi.nlm.nih.gov/articles/PMC8728285/
st.markdown("### Input")
st.markdown("##### Aqueous environment parameters")

# Create two columns; adjust the ratio to your liking
col1, col2 = st.columns([2,2]) 

# Use the first column for text input
with col1:
    p_h_input = st.text_input("PH", 7.4)

# Use the second column for the submit button
with col2:
    temperature_c = st.text_input("Temperature in Celsius (°C)", 25)
    


# kelvin  (25 degree)
# The Kelvin temperature scale is an absolute scale where 0 K (Kelvin) represents absolute zero, 
# the point at which all molecular motion ceases 
# 0°C is equal to 273.15 K, and 100°C is equal to 373.15 K. 
temperature_k = float(temperature_c) + 273.15

# mg
#mg_input = st.text_input("Concentration of Mg2 + ions", 3.0)
#ionic_strength_input = st.text_input("""Ionic Strength – total concentration of ions in a solution)""", 0.25M)
#p_mg = Q_(3.0)

# Ionic strength \(\left(I\right)\) is a measure of the total concentration of ions in a solution. 
# 0.25 M K2SO4 Potassium sulfate
# ionic_strength = Q_("0.25M")



# 2. reaction 
# input

st.markdown("##### Template (Example: 1A + 1B = 1C + 1D)")

# Create two columns; adjust the ratio to your liking
col3, col4 = st.columns([2,2]) 

# Use the first column for text input
with col3:
    num1 = st.text_input("Stoichiometric coefficient A", "1")
    num2 = st.text_input("Stoichiometric coefficient B", "1")
    num3 = st.text_input("Stoichiometric coefficient C", "1")
    num4 = st.text_input("Stoichiometric coefficient D", "1")

# Use the second column for the submit button
with col4:
    meta1 = st.text_input("Metabolite A:", "atp")
    meta2 = st.text_input("Metabolite B:", "h2o")
    meta3 = st.text_input("Metabolite C:", "adp")
    meta4 = st.text_input("Metabolite D:", "pi")


if meta1 and meta2 and meta3 and st.button("Submit"):
    cc = ComponentContribution()
    
    # ph 
    cc.p_h = Q_(p_h_input)
    cc.temperature = Q_(f"{temperature_k}K")
    
    # Parsing reaction with non-trivial stoichiometric coefficients is simple. Just add the coefficients before each compound ID (if none is given, it is assumed to be 1)
    input_meta1=prepare_input(meta1)
    input_meta2=prepare_input(meta2)
    input_meta3=prepare_input(meta3)


    # reaction by name 
    if meta4 and num4:
        input_meta4=prepare_input(meta4)
        atpase_reaction_name = f"{num1} {meta1} + {num2} {meta2} = {num3} {meta3} + {num4} {meta4}"
        atpase_reaction = cc.parse_reaction_formula(
            f"{num1} {input_meta1} + {num2} {input_meta2}={num3} {input_meta3} + {num4} {input_meta4}"
        )
    else:
        atpase_reaction_name = f"{num1} {meta1} + {num2} {meta2} = {num3} {meta3}"
        atpase_reaction = cc.parse_reaction_formula(
            f"{num1} {input_meta1} + {num2} {input_meta2}={num3} {input_meta3}"
            )


    # Result
    with st.expander("Result"):
        st.markdown("#### The Reaction")
        st.markdown(f"""
                    - Conditions: PH {p_h_input}; Kelvin Temperature: {temperature_k}
                    - {atpase_reaction_name}
                    - {atpase_reaction}

                    """) 
        st.markdown("#### Is the reaction balanced?")


        # balance result
        if atpase_reaction.is_balanced():
            reaction_balance = "balanced"
        else:
            reaction_balance = "not balanced"
        st.markdown(f"""
                    - The reaction is {reaction_balance}.
                    """)
            
        st.markdown("#### Estimate the Gibbs free energy (ΔG)")
        dG0_prime = cc.standard_dg_prime(atpase_reaction)
        st.markdown(f"ΔG'° = {dG0_prime}")
        
        if dG0_prime < 0:
            st.markdown(f"The reactions is likely to be spontaneous, without the need for external energy input.")
       #else:
        #    st.markdown(f"The reactions is likely requiring external energy input.")
        
        
        st.markdown("#### How is the reversibility of the reaction?")
        # The reversibility index is a measure of the degree of the reversibility of the reaction that is normalized for stoichiometry.
        # If you are interested in assigning reversibility to reactions we recommend this measure because 1:2 reactions are much “easier” to reverse than reactions with 1:1 or 2:2 reactions. You can see our paper for more information.
        # takes the mean concentration of metabolites and the number of substrates and products into consideration
        #  The reversibility index can help modelers decide which reactions are reversible in physiological conditions.
        #It helps determine if a reaction is likely to be reversible under physiological conditions, even if the standard Gibbs free energy change (ΔG°) suggests otherwise. This is because ΔG° doesn't account for the impact of reactant and product concentrations on the overall free energy change (ΔG'). 
        # If a reaction has a high reversibility index, it suggests that the reaction is likely to be reversible under typical cellular conditions, meaning it can proceed in both the forward and reverse directions. 
        st.markdown(f"""
                    - ln(Reversibility Index) = {cc.ln_reversibility_index(atpase_reaction)}
                    """)