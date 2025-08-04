import pandas as pd
from Bio import Entrez
import os
# fetch details in pubmed
def pubmed_fetch_details(id_list):
    ids = ','.join(id_list)
    # ******Define your email to use with NCBI Entrez
    Entrez.email = "your_email@sample.com"

    handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    # Create a list to hold our article details
    articles = []

    for pubmed_article in records['PubmedArticle']:
        article = {}
        article_data = pubmed_article['MedlineCitation']['Article']
        try:
            article['Year'] = article_data.get('ArticleDate', {})[0].get('Year')
        except:
            article['Year'] = "Not_known"
        article['ELocationID'] = str(article_data.get('ELocationID', {})[0])    
        article['Title'] = article_data.get('ArticleTitle')
        
        
        # Directly output the abstract
        abstract_text = article_data.get('Abstract', {}).get('AbstractText', [])
        if isinstance(abstract_text, list):
            abstract_text = ' '.join(abstract_text)
        article['Abstract'] = abstract_text

        article['Journal'] = article_data.get('Journal', {}).get('Title')

        articles.append(article)
        
    # Convert our list of articles to a DataFrame
    DF = pd.DataFrame(articles)

    
    return DF


