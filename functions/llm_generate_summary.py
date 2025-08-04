# import
from langchain_core.prompts import ChatPromptTemplate

# Construct the prompt with the potentially trimmed node_names_text
def llm_generate_summary(llm, filepath, keyword_network="", trimmed_node_names_text="", network_summary="yes"):
    if network_summary=="yes":
    
        prompt_text = "These are the terms related to " + filepath + keyword_network + ", categorize them and write a summary report." + trimmed_node_names_text
    else:
        prompt_text = "These are the terms related to " + filepath + keyword_network + ", write a summary report." 
    system_text = f"""This GPT is specialized for analyzing scientific results"""
    
    prompt = ChatPromptTemplate.from_messages(
            [
                (
                    "system",
                    system_text,
                ),
                ("human", "{input}"),
            ]
        )
         
    chain = prompt | llm
    
    response = chain.invoke(
        {
            "input": prompt_text,
        }
    )

    return response.content.strip()
