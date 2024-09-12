from openai import OpenAI
import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Fetch environment variables
LIST_DIR = os.getenv('LIST_DIR')  # Directory with the text files
GOOD_RESULT_DIR = os.getenv('GOOD_RESULT_DIR')  # Directory to save good results
BAD_RESULT_DIR = os.getenv('BAD_RESULT_DIR')  # Directory to save bad results
OPENAI_API_KEY = os.getenv('OPENAI_API_KEY')
client = OpenAI(api_key=OPENAI_API_KEY)

# Ensure result directories exist
os.makedirs(GOOD_RESULT_DIR, exist_ok=True)
os.makedirs(BAD_RESULT_DIR, exist_ok=True)

# Function to call OpenAI's GPT model
def evaluate_content(content):
    response = client.chat.completions.create(model="gpt-3.5-turbo",  # You can use another model if needed
    messages=[
        {"role": "system", "content": "You are an AI that helps evaluate the quality of text content."},
        {"role": "user", "content": f"Please evaluate the following content and tell me if it's good information:\n\n{content}"}
    ])
    return response.choices[0].message.content

# Function to determine whether the content is good or bad based on the response
def is_good_information(evaluation):
    # Adjust this condition based on how OpenAI's response indicates good/bad information
    return "good" in evaluation.lower()

# Function to process files and save original and evaluation results
def process_files():
    for filename in os.listdir(LIST_DIR):
        file_path = os.path.join(LIST_DIR, filename)

        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()

        # Evaluate content using the OpenAI GPT API
        evaluation = evaluate_content(content)

        # Determine if the content is good or bad
        if is_good_information(evaluation):
            result_dir = GOOD_RESULT_DIR
        else:
            result_dir = BAD_RESULT_DIR

        # Save the original content and the GPT response in the respective directory
        original_file_path = os.path.join(result_dir, f"{filename}")
        evaluation_file_path = os.path.join(result_dir, f"evaluation_{filename}")

        # Save the original file content
        with open(original_file_path, 'w', encoding='utf-8') as original_file:
            original_file.write(content)

        # Save the GPT evaluation response
        with open(evaluation_file_path, 'w', encoding='utf-8') as eval_file:
            eval_file.write(evaluation)

        print(f"Processed and saved result for: {filename} in {result_dir}")

# Run the process
if __name__ == "__main__":
    process_files()
