# protocol ep3 

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
st.title("Metabolic Model")

#------------------------------------------#
# show animated map 
#------------------------------------------#
flux_html_name = "flux_map.html"
flux_html_path = os.path.join(data_dir, flux_html_name)


with st.expander("Flux Balance Analysis"):
    HtmlFile = open(flux_html_path, 'r', encoding='utf-8')
    source_code = HtmlFile.read() 
    components.html(source_code,height=800, scrolling=True)
                