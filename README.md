# An Automatic Pipeline Approach for Semantic Exploration of Genetic Variant-Disease Associations

The code of this repository is a pipeline approach that automatically integrates diverse biomedical data-bases to enable semantic exploration of genetic variants in disease associations. The pipeline leverages multiple data sources, including NCBI Gene, MeSH, LitVar2, PubTator and SynVar databases, to retrieve comprehensive information about genes, variants, diseases and associated literature. The pipeline consists of multiple stages, encompassing querying and searching across the different databases, extracting relevant data and applying filters to refine the results. 

It aims to bridge the gap between genetic variant information and disease associations by providing a systematic framework for discovering and analyzing relevant literature. Although the pipeline has limitations, it still uncovers additional articles not referenced in expert reports that mention the genetic variants of interest. 

##Process
Given a variant, a gene and a disease:
1. The pipeline begins by querying the NCBI Gene database using the gene of interest to retrieve gene information; 
2. Using the gene aliases obtained in the previous step, the pipeline queries the SynVar database to generate a list of unique syntactic variations 
associated with the specified variant; 
3. The pipeline retrieves and stores disease synonyms from the MeSH database based on the specified disease; 
4. The pipeline searches the LitVar2 database for variant-associated publications. 
It retrieves articles that mention the gene and the variant; 
5. Finally, the pipeline filters the retrieved articles based on disease associations. 
Articles that do not mention the specified disease or its synonymous terms are excluded from further analysis.

##Inputs and outputs
Input: A csv file with the following columns: variant, gene, disease. Each row is a different variant.
Output: A folder for each variant with the following files:
  - data.json: a json file with the following information:
      - variant: the variant
      - gene: the gene
      - gene_aliases: a list with the gene aliases
      - syntactic_variations: a list with the syntactic variations
      - disease: the disease
      - disease_synonyms: a list with the disease synonyms
      - pmids: a list with the pmids of the articles retrieved from litvar2
      - pmcids: a list with the pmcids of the articles retrieved from litvar2
  - log.txt: a log file with the following information:
      - date and time of the execution
      - the variant, gene and disease of the query
      - the number of articles retrieved from litvar2
      - the number of articles retrieved from pubtator
  For each article retrieved from litvar2, creates a folder with the pmid of the article, containing the following files:
      - pubtator_annotations.json: a json file with the annotations of the articles retrieved from pubtator, for each article.
  - filter.json: a json file with the list of articles that mention the disease or its synonyms, for each article.

##Requirements
This pipeline requires an Entrez API key. To get one, go to https://www.ncbi.nlm.nih.gov/account/settings/ and create an API key.

##Instructions

