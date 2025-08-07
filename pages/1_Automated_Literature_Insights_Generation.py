import sys
sys.path.append('../')

# functions
## search
from functions.pubmed_search import *
from functions.pubmed_fetch_details import *
from functions.prepare_tag import *
from functions.download_csv import *

## llm
from functions.prepare_llm import *
from functions.llm_ask_a_question import *
from functions.llm_generate_summary import *

## network
from functions.trim_text import trim_text
from functions.network_generate import *
from functions.network_search import * 
from functions.network_write_to_html import *
from functions.network_get_node_names import network_get_node_names

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
print(output_dir)

# Form title
st.title("🔎 Automated Literature Search and Insight Generation")

# Form 
with st.form("my_form",clear_on_submit=False):
    st.markdown(body="## 1. Literature Search", unsafe_allow_html=False, help=None, width="stretch")
    #------------------------------------------#
    # Input
    #------------------------------------------#
    keyword = st.text_input("1. Enter Keywords:", "AI AND Pathways AND Retrosynthesis")
   
    # slider
    max_records=st.slider("Maximal Records", 1, 100, 2)
    
    # question to llm
    system_text = f"""You are specialized for analyzing scientific paper abstracts, \
    focusing on identifying specific entities related to biological studies, \
    such as pathways, species, genes, methods of genetic engineering, enzymes, proteins, performance, \
    and bioprocess conditions (e.g., growth conditions), and determining causal relationships between them. \
    It outputs all possible combinations of causal relationships among identified entities in structured pairs. \
    The output strictly follows the format: (Entity A ,Entity B), with no additional text."""
    
    #------------------------------------------#
    # Frontend 
    #------------------------------------------#
    st.markdown(body="## 2. LLM Insight Generation (Optional):", unsafe_allow_html=False, help=None, width="stretch")
    
    # 2.1 search literature
    st.markdown(body="### 2.1 Ask a question to Agent:", unsafe_allow_html=False, help=None, width="stretch")
    
    # question to llm text area
    question_input = st.text_area("Question to Agent:", system_text)
 
    # use llm checkbox
    checkbox_val = st.checkbox("Ask llm")

    # 2.2 network
    st.markdown(body="### 2.2 Generate a Network", unsafe_allow_html=False, help=None, width="stretch")
    keyword_network  = st.text_input("Filter network by filter keywords:", "AI")
    filter_network_by_keywords = st.checkbox("Filter network by the above filter keywords")
    
    # 2.3 network summary
    st.markdown(body="### 2.3 Generate a Network Summary", unsafe_allow_html=False, help=None, width="stretch")
    generate_network_summary = st.checkbox("Generate a network summary")
    
    # submit 
    submit_button = st.form_submit_button("Submit")
   
    #------------------------------------------#
    # Backend
    #------------------------------------------#
    # if click download
    if submit_button:
        
        # remove space from keyword
        key_word_no_space = keyword.replace(' ', '_')
        
        # tag 
        tag_string = "pubmed_"+ prepare_tag(key_word_no_space)
        

        
        #------------------------------------------#
        # 0. Search literature
        #------------------------------------------#
        # pubmed search 
        id_list = pubmed_search(keyword, retmax=max_records)
        df_articles = pubmed_fetch_details(id_list)
        
        # checkbox: use llm to generate insights or not
        # no llm
        if not checkbox_val:
            # csv name
            csv_name = f"{tag_string}.csv"
            # save a copy 
            csv_path = os.path.join(output_dir, csv_name)
            df_articles.to_csv(csv_path,index=False)
            
            # show df 
            with st.expander("List of Papers"):
                st.dataframe(df_articles)
                
                # download csv
                st.markdown(get_binary_file_downloader(csv_path, 'csv'), unsafe_allow_html=True)
        else: 
            # llm
            # Save the DataFrame with new responses back to an csv file
            csv_name =  f"{tag_string}_llm_response.csv"
    
            #------------------------------------------#
            # 1. llm answer a question
            #------------------------------------------#
            start_row = 0  # Adjusted to start from row number xxx
            end_row = len(df_articles)
            
            # llm
            llm = prepare_llm(llm_model_id ="llama-3.3-70b-versatile", 
                temperature=0, 
                max_tokens=2000)
            
            for i in range(start_row, 
                end_row):  
                row = df_articles.iloc[i]
                abstract = row['Abstract']
                response = llm_ask_a_question(llm =llm, 
                                                    input_text=abstract, 
                                                    question=question_input, 
                                                    system_text=system_text)
                
                df_articles.at[i, 'Answer to Question'] = response
                
                # Optional: Print the response and progress to monitor execution
                #print(f"Row {i+1} Response: {response}")
                #progress = ((i + 1 - (start_row - 1)) / (end_row - (start_row - 1))) * 100
                #st.write(f"Progress: {progress:.1f}% completed\n")
                
                # save csv
                csv_path = os.path.join(output_dir, csv_name)
                df_articles.to_csv(csv_path,index=False)
                
            
            with st.expander("List of Papers"):
                st.dataframe(df_articles)
                
                # download csv
                st.markdown(get_binary_file_downloader(csv_path, 
                                                       'csv'), 
                            unsafe_allow_html=True)
            
            #------------------------------------------#
            # 2. llm generate network
            #------------------------------------------#
            #if download_html_checkbox:
            Graph = network_generate(DF = df_articles, 
                                            source_col ='Title',
                                            value_col='Answer to Question')

            if filter_network_by_keywords:
                # filter
                Graph_out = network_search(Graph, keyword_network)
                # Save and show the network
                filterd_html_name =f"{tag_string}_network_filtered_by_{keyword_network}.html" 
                
            else:
                Graph_out = Graph
                
            # Save and show the network
            html_name =f"{tag_string}_network.html" 
            html_path = os.path.join(output_dir, html_name) 
            network_write_to_html(input_graph=Graph,html_path=html_path)
            
            # node names (use for network summary)
            node_names_text =network_get_node_names(input_graph =Graph)
            
            # display network and download link
            with st.expander("Network Analysis using Abstracts (Zoom in/out)"):
                HtmlFile = open(html_path, 'r', encoding='utf-8')
                source_code = HtmlFile.read() 
                components.html(source_code,height=800, scrolling=True)
                
                # download link
                st.markdown(get_binary_file_downloader(html_path, 'html'), unsafe_allow_html=True)
                    
            #------------------------------------------#
            # 3. llm generate summary of network
            #------------------------------------------#
            # The code snippet you provided is calling the function `network_get_node_names` with the following
            # parameters:
            # - `input_graph`: This parameter is set to `csv_path`, which is the path to the CSV file
            # containing the response data from the literature search and LLM processing.
            # - `cut_off_chunk_size`: This parameter is set to `30000`, which specifies the maximum size of each
            # chunk of text when trimming the node names.
            # - `trim_node_name_text_by_chunk`: This parameter is set to `"yes"`, indicating that the function
            # should trim the node names by chunk.
                        # Apply the trimming function to node_names_text
            if generate_network_summary:
                cut_off_chunk_size = 30000 
                trimmed_node_names_text = trim_text(node_names_text, 
                                                    max_length =cut_off_chunk_size)
                
                # # Define your question, we already set the system prompt so left blank here
                question_summary = " "
                
                # generate a summary related to network
                response = llm_generate_summary(llm=llm,
                                        filepath =csv_path, 
                                        keyword_network=keyword_network, 
                                        trimmed_node_names_text=trimmed_node_names_text,
                                        network_summary="yes")
                # save markdown
                md_name= f"{tag_string}_network.md" 
                md_path= os.path.join(output_dir, md_name)
                with open(md_path, 'w') as f:
                    f.write(response)
                
                # convert md to word
                word_name = f"{tag_string}_network.docx" 
                word_path = os.path.join(output_dir, word_name)
                markdown_to_word(md_path, word_path)
                    
                with st.expander("Summary of Network Components"):
                    st.markdown(response, unsafe_allow_html=True)
                    st.markdown(get_binary_file_downloader(md_path, 'md'), unsafe_allow_html=True)
                    st.markdown(get_binary_file_downloader(word_path, 'docx'), unsafe_allow_html=True)





