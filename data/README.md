# Data 📊

**Note** some of the data files in this directory are git*ignored*.

Contents:
- `genes.list` (*ignored*)
  - purpose
    - identify the genes to scrape from pubmed
  - file properties:
    - 1 gene per row
    - keep all caps if possible (case shouldn't matter, but hasn't been tested)
    - duplicates will be ignored
- `gene_data.csv` (*ignored*)
  - purpose
    - table of the sequence information and coding regions of the genes
  - file properties:
    - header column
    - 1 gene per row
    - data has been validated & exceptions dealt with
    - minimal preprocessing has been done on the data files (check CDS starts)
