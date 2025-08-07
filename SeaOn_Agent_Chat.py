import sys
sys.path.append('./')

## llm
from functions.prepare_llm import *
from functions.llm_ask_a_question import *

import streamlit as st

#with st.sidebar:
#    openai_api_key = st.text_input("OpenAI API Key", key="chatbot_api_key", type="password")
#    "[Get an OpenAI API key](https://platform.openai.com/account/api-keys)"
#    "[View the source code](https://github.com/streamlit/llm-examples/blob/main/Chatbot.py)"
#    "[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/streamlit/llm-examples?quickstart=1)"

st.title("SeaOn Bio-Design Agent")
st.caption("🚀 Powered by LLM")
if "messages" not in st.session_state:
    st.session_state["messages"] = [{"role": "assistant", "content": "How can I help you?"}]

for msg in st.session_state.messages:
    st.chat_message(msg["role"]).write(msg["content"])

if prompt := st.chat_input():
    # llm
    llm = prepare_llm(llm_model_id ="llama-3.3-70b-versatile", 
        temperature=0, 
        max_tokens=2000)

    st.session_state.messages.append({"role": "user", "content": prompt})
    st.chat_message("user").write(prompt)
    
    
    #response = llm.chat.completions.create(model="gpt-3.5-turbo", messages=st.session_state.messages)
    
    question_input = st.session_state.messages
    # question to llm
    system_text = f"""You are specialized for answering questions, and if needed collecting information online and generate \
        sophisticated insights."""
    
    response = llm_ask_a_question(llm =llm, 
                                  input_text=" ", 
                                  question=question_input, 
                                  system_text=system_text)
    
    msg = response
    st.session_state.messages.append({"role": "assistant", "content": msg})
    st.chat_message("assistant").write(msg)
