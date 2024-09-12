from Bio import Entrez
import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Fetch environment variables
LIST_FILE = os.getenv('LIST_FILE')  # The file containing PubMed IDs
SAVED_FILE_DIR = os.getenv('LIST_DIR')
EMAIL = os.getenv('EMAIL')
Entrez.email = EMAIL

# Function to fetch article details using a list of PubMed IDs
def fetch_articles(id_list):
    fetch_handle = Entrez.efetch(db='pubmed', retmode='xml', id=','.join(id_list))
    papers = Entrez.read(fetch_handle)['PubmedArticle']
    fetch_handle.close()
    return papers

# Function to save article details to a text file
def save_article(article_id, title, abstract, journal, language, year, month):
    if not os.path.exists(SAVED_FILE_DIR):
        os.makedirs(SAVED_FILE_DIR)
    
    file_path = os.path.join(SAVED_FILE_DIR, f"{article_id}.txt")
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(f"PubMed ID: {article_id}\nTitle: {title}\nJournal: {journal}\nLanguage: {language}\n")
        f.write(f"Year: {year}\nMonth: {month}\n\nAbstract:\n{abstract}\n")

# Main function to process articles based on an existing list of IDs
def process_articles_from_list():
    with open(LIST_FILE, 'r') as file:
        id_list = [line.strip() for line in file]  # Read IDs from the file

    papers = fetch_articles(id_list)
    saved_count = 0  # Counter for saved articles
    
    for article_id, paper in zip(id_list, papers):
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
    process_articles_from_list()
