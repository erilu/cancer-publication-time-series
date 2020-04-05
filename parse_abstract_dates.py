import pandas as pd
import re
import calendar

# read in a csv file of abstracts from pubmed_extractor.py
abstract_df = pd.read_csv("abstracts.csv", usecols=[0])

# parse dates from each abstract
abstract_df['year'] = [int(re.search(".\s(\d{4})\s(\w{3})",row).group(1)) for row in abstract_df['Journal']]
abstract_df['month'] = [re.search(".\s(\d{4})\s(\w{3})",row).group(2) for row in abstract_df['Journal']]
abstract_df['month_num'] = [list(calendar.month_abbr).index(x) for x in abstract_df['month']]

# count number of abstracts for each month of each year from min(year) to max(year)
papers_by_month_df = pd.DataFrame(columns=['year','month','count'])
for year in range(min(abstract_df['year']),max(abstract_df['year'])+1):
    entries = abstract_df.loc[abstract_df['year']==year,]
    for month in range(1,13):
        count = len(entries.loc[entries['month_num']==month,])
        papers_by_month_df = papers_by_month_df.append(pd.Series([year,month,count], index = papers_by_month_df.columns), ignore_index=True)

# store counts per month in a csv file
papers_by_month_df.to_csv("papers_published_per_month.csv")
