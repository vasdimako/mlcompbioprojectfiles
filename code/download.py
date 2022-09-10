import pandas as pd
from multiprocessing.pool import ThreadPool
import requests
from os.path import exists


def download_url(args):
    filename, url = args[0], args[1]
    if not exists(filename):
        try:
            r = requests.get(url)
            with open(filename, 'wb') as f:
                f.write(r.content)
            return url
        except Exception as e:
            print('Exception in download_url():', e)


def download_parallel(args):
    results = ThreadPool(4).imap_unordered(download_url, args)
    for result in results:
        print('url:', result)


data_df = pd.read_csv("tbi_data.csv")

data_df["UID"] = data_df["donor_name"].apply(lambda x: x.split(".")[2]) + "_" + data_df["structure_acronym"].apply(lambda x: str(x)) + "_" + data_df["rna_well"].apply(lambda x: str(x))
bam_urls = []
bai_urls = []

for ind, row in data_df.iterrows():
    if row["structure_acronym"] == "FWM":
        bam_urls.append(["D:/ML/data/FWM/" + row["UID"] + ".sorted.bam", "http://aging.brain-map.org" + row["anonymized_bam_file_link"]])
        bai_urls.append(["D:/ML/data/FWM/" + row["UID"] + ".bai", "http://aging.brain-map.org" + row["anonymized_bam_index_file_link"]])

download_parallel(bam_urls)