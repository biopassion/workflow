import sys
sys.path.append('../')

# functions
## download
from functions.get_binary_file_downloader import *
from md2docx_python.src.md2docx_python import markdown_to_word

## streamlit 
import streamlit as st
import streamlit.components.v1 as components

# the current dir is where the main app.py is 
main_dir= "./"
# output dir
output_dir=os.path.join(main_dir,"output")
data_dir = os.path.join(main_dir,"data")


# Form title
st.title("Characterization of a Metabolite")

#------------------------------------------#
# load a model 
#------------------------------------------#
from cobra.io import read_sbml_model
model_name = 'core.xml'
model_path = os.path.join(data_dir, model_name)
model=read_sbml_model(model_path)


with st.form("my_form",clear_on_submit=False):
    metabolite = st.text_input("Enter a specific metabolite using metabolite identifier:", "g3p_c")
    
    submit_button = st.form_submit_button("Submit")
    
    reaction_list = []
    for reaction in model.metabolites.get_by_id(metabolite).reactions:
        reaction_list.append([reaction.id, reaction, reaction.name])
        
    if submit_button:
        with st.expander("Result:"):
            st.markdown(f"The name of {metabolite} is:" + "\n" +model.metabolites.get_by_id(metabolite).name)
            st.markdown(f"The compartment of {metabolite} is: \n" + model.metabolites.get_by_id(metabolite).compartment)
            st.markdown(f"The formula of {metabolite} is:" + model.metabolites.get_by_id(metabolite).formula)
            st.markdown(f"The elements of {metabolite} include:" + str(model.metabolites.get_by_id(metabolite).elements))
            st.markdown(f"The molecular weight of {metabolite} is:" + str(model.metabolites.get_by_id(metabolite).formula_weight))
            st.markdown(f"The annotation of {metabolite} is:" + str(model.metabolites.get_by_id(metabolite).annotation))
            st.markdown(f"The metabolite {metabolite} is involved in these reactions:" + str(model.metabolites.get_by_id(metabolite).reactions))
            st.markdown(f"The reactions details:" + str(reaction_list))
       

        
        

        


