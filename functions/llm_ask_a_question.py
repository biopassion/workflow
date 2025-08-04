import pandas as pd
from langchain_core.prompts import ChatPromptTemplate

# function to ask question
def llm_ask_a_question(llm, input_text, question, system_text):
    if len(input_text) > 3:
        prompt_text = question + " " + str(input_text)
    else:
         prompt_text = question

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