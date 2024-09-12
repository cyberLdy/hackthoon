setup.env: parameters
extract_by_query.py: extract pubmed article content and save to local files by taking $query and $NUM_ARTICLE from env, and save it to $QUERY_FILE_DIR
extract_by_list.py:extract pubmed article content and save to local files by provided list_id $LIST_FILE, and save it to $LIST_DIR
filter.py: use openai_key to simple evaulate each document content provided from $LIST_DIR , then filter and save the feedback to $GOOD_RESULT_DIR and $BAD_RESULT_DIR
                  [the prompt will be customized and change it later accordingly]                  [the filter step will be customized and change it later accordingly]
