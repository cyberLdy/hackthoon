import requests
from dotenv import load_dotenv
import os
import xml.etree.ElementTree as ET

# Load environment variables from .env file
load_dotenv()
NUM_ARTICLE = int(os.getenv('NUM_ARTICLE')) # 
QUERY = os.getenv('QUERY')  # 

def get_pubmed_ids(QUERY, NUM_ARTICLE):
    # API URL with parameters
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&sort=relevance&term={QUERY}&retmax={NUM_ARTICLE}"
    response = requests.get(url)
    if response.status_code == 200:
        root = ET.fromstring(response.content)
        # Extract all <Id> elements
        id_list = [id_elem.text for id_elem in root.findall(".//Id")]
        return id_list
    else:
        print(f"Error: {response.status_code}")



#EXAMPLE
-------------------------------------------------------------
        # Define Prompt Template directly in the script
        PROMPT_TEMPLATE = """
        Answer the question based only on the following context:

        {context}

        ---

        Answer the question based on the above context: {question}

        Return no answer if you cannot find an answer from the context.
        """

        #PROMPT_TEMPLATE = " Please write a complete, well-structued summary for this document"
        prompt_template = ChatPromptTemplate.from_template(PROMPT_TEMPLATE)
        prompt = prompt_template.format(context=context_text, question=query_text)

        # Load Remote LLM
        ollama = Ollama(
            base_url=os.getenv('OLLAMA_URL'),
            model=os.getenv('OLLAMA_MODEL'),
        )
        response = ollama.invoke(prompt)
--------------------------------------------------------------------------
if __name__ == "__main__":
    id_list  = get_pubmed_ids(QUERY, NUM_ARTICLE)
    for x in id_list:
        print(x)
