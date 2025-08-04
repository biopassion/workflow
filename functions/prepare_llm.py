import os
from langchain_groq import ChatGroq
import streamlit as st

try:
    os.environ['GROQ_API_KEY'] = st.secrets['GROQ_API_KEY']
except:
    import dotenv
    dotenv.load_dotenv()
    

def prepare_llm(llm_model_id ="llama-3.3-70b-versatile", 
                temperature=0, 
                max_tokens=2000):

    # Use llm model from Groq 
    # Get an API here: https://console.groq.com/keys
    llm=ChatGroq(model_name=llm_model_id,
                temperature=temperature,
                max_tokens=max_tokens   )
    return llm


