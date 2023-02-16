# import webbrowser

# start_codon = sequen.find('M')
# stop_codon = sequen[start_codon: ].find('_')
    
# for all in range(len(sequen)):
#     cds_ = ''
#     start_codon = sequen.find('M')
#     stop_codon = sequen[start_codon: ].find('_')
#     for i in range(start_codon, stop_codon):
#         aminoCodon = sequen[i:i+3]
#         cds_ += aminoCodon

# print(sequen.split('_'))
from bs4 import BeautifulSoup
import requests

# from urllib.request import urlopen
# from urllib.parse import quote

# protein_url = 'https://www.ncbi.nlm.nih.gov/protein/?term=cry'
# nucleotide_url = 'https://www.ncbi.nlm.nih.gov/nuccore/?term=cry'
# gene_url = https://www.ncbi.nlm.nih.gov/gene/?term=cry
# out = urlopen(url).read().decode('utf-8')
# webbrowser.open(url)


def get_ID(query, database_name='protein'):

    if database_name == 'protein':
        url = f'https://www.ncbi.nlm.nih.gov/protein/?term={query}'
    elif database_name == 'neucleotide':
        url = f'https://www.ncbi.nlm.nih.gov/nuccore/?term={query}'
    elif database_name == 'gene':
        url = f'https://www.ncbi.nlm.nih.gov/gene/?term={query}'

    try:
        soup = BeautifulSoup(requests.get(url).text, 'lxml')
        id_list = []
        for dd in soup.select('.rprtid dd'):
            id_list.append(dd)
        return id_list

    except:
        return 'ERR: reterive failed:\n1.check your internet connection first\n2.access wrong database or mistakes!\nyou can only access database like protein, gene, neucleotide'


print(get_ID('cry'))