
from Bio import Entrez

def pubmed_search(keyword, retmax):  

    
    # ******Define your email to use with NCBI Entrez
    Entrez.email = "your_email@sample.com"
    
    # Adjust the search term to focus on abstracts
    search_term = f"{keyword}[Abstract]"
    #****** You can modify retmax to change maximum number of search results
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=retmax) #******
    record = Entrez.read(handle)
    handle.close()
    # Get the list of Ids returned by the search
    id_list = record["IdList"]
    return id_list