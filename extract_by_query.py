from Bio import Entrez
import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Fetch environment variables
NUM_ARTICLE = int(os.getenv('NUM_ARTICLE'))
QUERY = os.getenv('QUERY')
QUERY_FILE_DIR = os.getenv('QUERY_FILE_DIR')
EMAIL = os.getenv('EMAIL')
Entrez.email = EMAIL

# Function to search PubMed and fetch article details
def fetch_articles(query, retmax):
    search_handle = Entrez.esearch(db='pubmed', sort='relevance', retmax=retmax, retmode='xml', term=query)
    id_list = Entrez.read(search_handle)['IdList']
    search_handle.close()
    
    fetch_handle = Entrez.efetch(db='pubmed', retmode='xml', id=','.join(id_list))
    papers = Entrez.read(fetch_handle)['PubmedArticle']
    fetch_handle.close()
    
    return id_list, papers

# Function to save article details to a text file
def save_article(article_id, title, abstract, journal, language, year, month):
    if not os.path.exists(QUERY_FILE_DIR):
        os.makedirs(QUERY_FILE_DIR)
    
    file_path = os.path.join(QUERY_FILE_DIR, f"{article_id}.txt")
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(f"PubMed ID: {article_id}\nTitle: {title}\nJournal: {journal}\nLanguage: {language}\n")
        f.write(f"Year: {year}\nMonth: {month}\n\nAbstract:\n{abstract}\n")

# Main function to process articles
def process_articles():
    ids, papers = fetch_articles(QUERY, NUM_ARTICLE)
    saved_count = 0  # Counter for saved articles
    
    for article_id, paper in zip(ids, papers):
        article = paper['MedlineCitation']['Article']
        title = article['ArticleTitle']
        abstract = article.get('Abstract', {}).get('AbstractText', ['No Abstract'])[0]
        journal = article['Journal']['Title']
        language = article['Language'][0]
        pubdate = article['Journal']['JournalIssue']['PubDate']
        year = pubdate.get('Year', 'No Data')
        month = pubdate.get('Month', 'No Data')
        
        save_article(article_id, title, abstract, journal, language, year, month)
        saved_count += 1  # Increment the counter
    
    print(f"Total articles saved: {saved_count}")

# Run the process
if __name__ == "__main__":
    process_articles()
