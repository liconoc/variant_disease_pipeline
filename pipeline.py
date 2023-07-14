# -*- coding: utf-8 -*-
"""
Created on 2023

@author: Lidia Contreras-Ochando

Given a variant, a gene and a disease:
(1) The pipeline begins by querying the NCBI Gene database using the gene of interest to retrieve gene information; 
(2) Using the gene aliases obtained in the previous step, the pipeline queries the SynVar database to generate a list of unique syntactic variations 
associated with the specified variant; 
(3) The pipeline retrieves and stores disease synonyms from the MeSH database based on the specified disease; 
(4) The pipeline searches the LitVar2 database for variant-associated publications. 
It retrieves articles that mention the gene and the variant; 
(5) Finally, the pipeline filters the retrieved articles based on disease associations. 
Articles that do not mention the specified disease or its synonymous terms are excluded from further analysis.

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

This pipeline requires an Entrez API key. To get one, go to https://www.ncbi.nlm.nih.gov/account/settings/ and create an API key.

"""

import datetime
import os
import requests
import json
from Bio import Entrez
from xml.etree import ElementTree


###################################################
###  (1) function to get the aliases of a gene  ###
###################################################

def get_aliases(folder_name, gene, variant, api_key):

    #Search the gene in gene
    handle = Entrez.esearch(db="gene", term=gene+"[Gene Name] AND Homo Sapiens[Organism]", api_key=api_key)
    record = Entrez.read(handle)
    handle.close()

    #list of ids in the record
    id_list=record['IdList']

    aliases=[]
    #for id in id_list, save the alias of the gene
    for gene_id in id_list:
        handle = Entrez.esummary(db="gene", id=gene_id, api_key=api_key)
        gene_data=Entrez.read(handle)
        handle.close()
        newaliases=gene_data['DocumentSummarySet']['DocumentSummary'][0]['OtherAliases']
        aliases.extend(newaliases.split(', '))    

    #Create a file including the variant and all the aliases of genes
    all_aliases = [gene] + aliases
    data = {
        "variant": variant,
        "gene": gene,
        "gene_aliases": all_aliases
    }
    file_path = os.path.join(folder_name, "data.json")
    with open(file_path, "w") as file:
        json.dump(data, file)

    #update log file
    log_file = open(folder_name+"/log.txt", "a")
    log_file.write("The file 'data.json' has been updated with the gene aliases.\n")
    log_file.close()


######################################################################
###  (2) function to get all the syntactic variations from synvar  ###
######################################################################
    

def get_synvar(folder):
    #Get the info from data.json
    with open(folder+"/data.json") as json_file:
        data = json.load(json_file)

    #Get the gene list from the json file
    gene_list=data['gene_aliases']
    variant=data['variant']
    syntactic_variation_list=[variant]    

    #For each gene in the list, get the syntactic variations
    for gene in gene_list:
        url = "http://goldorak.hesge.ch/synvar/generate/literature/fromMutation?variant="+variant+"&ref="+gene
        response = requests.get(url)

        if response.status_code == 200:            
            data_synvar = response.text  
            root = ElementTree.fromstring(data_synvar)
            
            #Get all the syntactic variations, if they are not in the list yet            
            syntactic_variations = root.findall(".//syntactic-variation-list/syntactic-variation")
            
            # Add each variation to the list           
            for variation in syntactic_variations:
                variation_clean=variation.text.replace("  ", " ").replace("\n", "")                 
                if variation_clean not in syntactic_variation_list:
                    syntactic_variation_list.append(variation_clean)

    
    # New data to add to the json file
    data["syntactic_variations"] = syntactic_variation_list
    
    with open(folder + "/data.json", "w") as json_file:
        json.dump(data, json_file, indent=4)

    #update log file
    log_file = open(folder+"/log.txt", "a")
    log_file.write("The file 'data.json' has been updated with the syntactic variations of the variant.\n")
    log_file.close()

##############################################################
###  (3) function to get all the information from litvar2  ###
##############################################################

def get_litvar2(folder):    
    #Open the variant folder and load data.json
    with open(folder+"/data.json") as json_file:
        data = json.load(json_file)

    #Get the gene list from the json file
    original_variant=data['variant']
    original_gene=data['gene']
    gene_list=data['gene_aliases']
    syntactic_variations=data['syntactic_variations']

    pmids_list_litvar=[]
    pmcids_list_litvar=[]    
    all_pmids_count_litvar=0
    all_pmcids_count_litvar=0

    #for each syntactic variation, get the information from litvar2 autocomplete
    for variant in syntactic_variations:        

        url = "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/?query="+variant
        response = requests.get(url)

        if response.status_code == 200:
            data_litvar2 = response.json()         
                
            if len(data_litvar2)>0:                
                #update log file
                log_file = open(folder+"/log.txt", "a")
                log_file.write(variant + " variant found in litvar2 autocomplete\n")
                log_file.close()
                
                #For each variant, get the rsid, gene and pmids_count from litvar2 autocomplete
                for i in range(len(data_litvar2)): 
                    
                    #initialize variables
                    rsid_litvar=""
                    gene_list_litvar=[]
                    pmids_count_litvar=0 
                    pmcids_count_litvar=0                      

                    if 'gene' in data_litvar2[i]:
                        gene_list_litvar=data_litvar2[i]['gene'] 
                    if 'pmids_count' in data_litvar2[i]:
                        pmids_count_litvar=data_litvar2[i]['pmids_count'] 
                    if 'pmcids_count' in data_litvar2[i]:
                        pmcids_count_litvar=data_litvar2[i]['pmcids_count']  
                    if 'rsid' in data_litvar2[i]:
                        rsid_litvar=data_litvar2[i]['rsid']
                                      
                    #if gene is in the list all_gene_list
                    if any(gene in gene_list_litvar for gene in gene_list): 
                        #if pmid_count>0, get the pmids and pmcids from litvar2
                        if(int(pmids_count_litvar)>0):
                            #Get the pmids from litvar2
                            url3 = "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/get/litvar%40"+rsid_litvar+"%23%23/publications"
                            response3 = requests.get(url3)
                            if response3.status_code == 200:
                                data_litvar_pub = response3.json()
                                pmids_list_litvar.append(data_litvar_pub['pmids'])
                                pmcids_list_litvar.append(data_litvar_pub['pmcids'])                           
                               
                    
    #Save the new data in data.json
    flattened_pmids_list_litvar = list(set(item for sublist in pmids_list_litvar for item in sublist))
    flattened_pmcids_list_litvar = list(set(item for sublist in pmcids_list_litvar for item in sublist))
    all_pmids_count_litvar=len(flattened_pmids_list_litvar)
    all_pmcids_count_litvar=len(flattened_pmcids_list_litvar)
    data = {
        "variant": original_variant,
        "gene": original_gene,
        "gene_aliases": gene_list,
        "syntactic_variations": syntactic_variations,
        "pmids_count": all_pmids_count_litvar,
        "pmcids_count": all_pmcids_count_litvar,
        "pmids": flattened_pmids_list_litvar,
        "pmcids": flattened_pmcids_list_litvar
    }

    with open(folder+"/data.json", "w") as outfile:
        json.dump(data, outfile, indent=4)

    #update log file
    log_file = open(folder+"/log.txt", "a")
    log_file.write("The file 'data.json' has been updated with the information from litvar2.\n")
    log_file.close()


#########################################################
###  (4) function to get the synonyms of the disease  ###
#########################################################

def get_disease_synonims(folder, api_key, disease):
    #Search the disease in mesh
    handle = Entrez.esearch(db="mesh", term=disease, api_key=api_key)
    record = Entrez.read(handle)
    handle.close()

    ids = record["IdList"]

    #for id in id_list, download the synonyms of the disease from mesh
    synonyms = []

    if len(ids)>0:
        for id in ids:
            handle = Entrez.esummary(db="mesh", id=id, api_key=api_key)
            data=Entrez.read(handle) 
            handle.close()            
            synonyms.extend(data[0]["DS_MeshTerms"])            
        
        #Append synonyms in data.json from folder
        with open(folder+"/data.json") as json_file:
            data = json.load(json_file)
            data["disease"] = disease
            data['disease_synonyms'] = synonyms

        with open(folder+"/data.json", 'w') as outfile:
            json.dump(data, outfile, indent=4)
        
        #update log file
        log_file = open(folder+"/log.txt", "a")
        log_file.write("Synonyms of disease "+disease+" added to data.json\n")
        log_file.close()
    else:
        #update log file
        log_file = open(folder+"/log.txt", "a")
        log_file.write("There are no synonyms of disease "+disease+" to add to data.json\n")
        log_file.close()


###########################################################
###  (5) function to get tagged entities from pubtator  ###
###########################################################

def get_pubtator(folder):
    #Open data.json     
    with open(folder+"/data.json") as json_file:
        data = json.load(json_file)

    #Get the pmids_count from the json file
    pmids=[]
    pmcids=[]
    if 'pmids' in data:
        pmids=data['pmids']
    if 'pmcids' in data:
        pmcids=data['pmcids']    

    #Download the annotated publications from pubtator
    if len(pmids)>0:
        #Pubtator maximum number of pmids is 100
        #If len is greater than 100, split the list in chunks of 100
        if len(pmids)>100:
            pmids_chunks = [pmids[x:x+100] for x in range(0, len(pmids), 100)]
            for pmids_chunk in pmids_chunks:
                url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml?pmids="+",".join(list(map(str, pmids_chunk)))+"&concepts=gene,disease,mutation"
                r = requests.get(url)
                #extract the annotations from the xml file
                tree = ElementTree.fromstring(r.content)
                for document in tree.findall('document'):
                    pmid=document.find('id').text
                    if not os.path.exists(folder+"/publications/"+pmid):
                        os.makedirs(folder+"/publications/"+pmid)
                    with open(folder+"/publications/"+pmid+"/pubtator_annotations.json", "w") as f:
                        annotations=[]
                        for passage in document.findall('passage'):
                            for annotation in passage.findall('annotation'): 
                                annotations.append({
                                    "text": annotation.find('text').text,
                                    "infons": {
                                        "type": annotation.find('infons').find('type').text,
                                        "identifier": annotation.find('infons').find('identifier').text
                                    },
                                    "locations": {
                                        "offset": annotation.find('location').attrib['offset'],
                                        "length": annotation.find('location').attrib['length']
                                    }
                                })
                        json.dump(annotations, f, indent=4)
        else:
            url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml?pmids="+",".join(list(map(str, pmids)))+"&concepts=gene,disease,mutation"
            r = requests.get(url)
            #extract the annotations from the xml file
            tree = ElementTree.fromstring(r.content)
            for document in tree.findall('document'):
                pmid=document.find('id').text
                if not os.path.exists(folder+"/publications/"+pmid):
                    os.makedirs(folder+"/publications/"+pmid)
                with open(folder+"/publications/"+pmid+"/pubtator_annotations.json", "w") as f:
                    annotations=[]
                    for passage in document.findall('passage'):
                        for annotation in passage.findall('annotation'):       
                            annotations.append({
                                "text": annotation.find('text').text,
                                "location": annotation.find('location').attrib,
                                "infons": [infon_data for infon in annotation.findall('infon') for infon_data in [{
                                    "key": infon.get('key'),
                                    "value": infon.text
                                }]]
                            })
                    json.dump(annotations, f)

        #update log file
        log_file = open(folder+"/log.txt", "a")
        log_file.write("Annotations for "+pmid+" downloaded.\n")
        log_file.close()        

    if len(pmcids)>0:
        #if length of pmids is greater than 0, download the annotated publications from pubtator
        if len(pmcids)>100:
            pmcids_chunks = [pmcids[x:x+100] for x in range(0, len(pmcids), 100)]
            for pmcids_chunk in pmcids_chunks:
                url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml?pmcids="+",".join(pmcids_chunk)+"&concepts=gene,disease,mutation"
                r = requests.get(url)
                #extract the annotations from the xml file
                tree = ElementTree.fromstring(r.content)
                for document in tree.findall('document'):
                    pmcid=document.find('id').text
                    if not os.path.exists(folder+"/publications/PMC"+pmcid):
                        os.makedirs(folder+"/publications/PMC"+pmcid)
                    with open(folder+"/publications/PMC"+pmcid+"/pubtator_annotations.json", "w") as f:
                        annotations=[]
                        for passage in document.findall('passage'):
                            for annotation in passage.findall('annotation'):     
                                annotations.append({
                                    "text": annotation.find('text').text,
                                    "infons": {
                                        "type": annotation.find('infons').find('type').text,
                                        "identifier": annotation.find('infons').find('identifier').text
                                    },
                                    "locations": {
                                        "offset": annotation.find('location').attrib['offset'],
                                        "length": annotation.find('location').attrib['length']
                                    }
                                })
                        json.dump(annotations, f, indent=4)
        else:
            url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml?pmcids="+",".join(pmcids)+"&concepts=gene,disease,mutation"            
            r = requests.get(url)
            #extract the annotations from the xml file
            tree = ElementTree.fromstring(r.content)
            for document in tree.findall('document'):
                pmcid=document.find('id').text            
                if not os.path.exists(folder+"/publications/PMC"+pmcid):
                    os.makedirs(folder+"/publications/PMC"+pmcid)
                with open(folder+"/publications/PMC"+pmcid+"/pubtator_annotations.json", "w") as f:
                    annotations=[]
                    for passage in document.findall('passage'):
                        for annotation in passage.findall('annotation'):     
                            annotations.append({
                                "text": annotation.find('text').text,
                                "location": annotation.find('location').attrib,
                                "infons": [infon_data for infon in annotation.findall('infon') for infon_data in [{
                                    "key": infon.get('key'),
                                    "value": infon.text
                                }]]
                            })
                    json.dump(annotations, f)

        #update log file
        log_file = open(folder+"/log.txt", "a")
        log_file.write("Annotations for PMC"+pmcid+" downloaded.\n")
        log_file.close()


############################################################
###  (6) function to filter the articles by the disease  ###
############################################################

def filter_by_disease(folder, api_key):
    results=[]
    #Open data.json and get the disease synonyms
    with open(folder+"/data.json") as json_file:
        data = json.load(json_file)
        if "disease_synonyms" in data:
            disease_synonyms = data['disease_synonyms']

            #If there are publications downloaded
            if os.path.exists(folder+"/publications") and len(os.listdir(folder+"/publications"))>0:
                #For each publication folder in publications folder, open pubtator_annotations.json and check if the disease synonyms are in the annotations
                for publication in os.listdir(folder+"/publications"):
                    print("Checking publication "+publication)
                    for file in os.listdir(folder+"/publications/"+publication):        
                        if file.endswith("pubtator_annotations.json"):
                            with open(folder+"/publications/"+publication+"/"+file) as json_file:  
                                cheked_diseases = []              
                                data = json.load(json_file)                
                                for annotation in data: 
                                    annotated_disease_synonyms=[]
                                    #If the annotation is a disease
                                    if "infons" in annotation and len(annotation["infons"])>0 and "value" in annotation["infons"][0]:                       
                                        text=annotation["infons"][0]["value"]   
                                        if "MESH" in text:
                                            print("Cheking annotation "+text)
                                            text_disease = text.split("MESH:")[1] 
                                            if text_disease in cheked_diseases:
                                                print("Disease "+text_disease+" already checked")
                                                break
                                            else:
                                                cheked_diseases.append(text_disease)                          
                                                handle = Entrez.esearch(db="mesh", term=text_disease+"[MeSH Unique ID]", api_key=api_key)
                                                record = Entrez.read(handle)
                                                handle.close()
                                                ids = record["IdList"]
                                                handle = Entrez.esummary(db="mesh", id=ids, api_key=api_key)
                                                data=Entrez.read(handle) 
                                                handle.close()    

                                                annotated_disease_synonyms.extend(data[0]["DS_MeshTerms"])  

                                                print("Cheking disease synonyms "+str(annotated_disease_synonyms))     

                                                #If any disease synonym is in the annotations, print the publication and the disease
                                                if any(disease_synonym in annotated_disease_synonyms for disease_synonym in disease_synonyms):
                                                        print("Publication "+publication+" contains the disease")
                                                        #append publication, disease and annotation to results
                                                        results.append({"publication":publication,"disease":annotation["text"],"annotation":text})
                                                        break
                                                else:
                                                    print("MESH not in the synonyms")
                                                    break

            #Save results
            with open("results/"+folder+"/filter.json", 'w') as outfile:
                json.dump(outfile, indent=4)
                print("Results saved in "+folder+"/filter.json")
        
                                

#######################
###  Main function  ###
#######################

def main():
    Entrez.email = "[PUT YOUR EMAIL HERE]"  
    api_key = "[PUT YOUR API KEY HERE]"  
    file_name = "[PUT YOUR FILE NAME HERE].csv"
 
    #Open test.csv file containing on each row a variant, a gene and a disease. Separated by ';'. The first row is the header.
    with open(file_name) as f:
        lines = f.readlines()
        for line in lines[1:]:
            variant=line.split(';')[0]
            gene=line.split(';')[1]
            disease=line.split(';')[2].rstrip('\n')
            print("Starting job for variant: "+variant+", gene: "+gene+" and disease: "+disease+" \n")

            #get the current date and hour to create the folder name
            now = datetime.datetime.now()
            current_date = now.strftime("%Y-%m-%d")
            current_hour = now.strftime("%H-%M-%S")

            #clean the variant and gene names to avoid problems with the folder name
            variant_cleaned = ''.join(e for e in variant if e.isalnum())
            gene_cleaned = ''.join(e for e in gene if e.isalnum())
            folder_name = f"results/{current_date}-{current_hour}_{variant_cleaned}-{gene_cleaned}"

            #create the folder
            os.mkdir(folder_name)

            #create a log file inside the folder
            log_file = open(folder_name+"/log.txt", "w")
            log_file.write("Log file for variant: "+variant+", gene: "+gene+" and disease: "+disease+" \n")
            log_file.close()

            #Get the aliases of the gene
            print("Getting the aliases of the gene...")
            get_aliases(folder_name, gene, variant, api_key)

            #Get the syntactic variations from synvar
            print("Getting the syntactic variations from synvar...")
            get_synvar(folder_name)

            #Get the information from litvar2
            print("Getting the information from litvar2...")
            get_litvar2(folder_name)            

            #Get synonyms of the disease
            print("Getting synonyms of the disease...")
            get_disease_synonims(folder_name, api_key, disease)

            #Get the tagged entities from pubtator
            print("Getting the tagged entities from pubtator...")
            get_pubtator(folder_name)

            #Filter the publications by disease
            print("Filtering the publications by disease...")
            filter_by_disease(folder_name, api_key)

            #update log file, append mode
            log_file = open(folder_name+"/log.txt", "a")
            log_file.write("Job finished!")
            log_file.close()

            print("Job finished!\n")

if __name__ == "__main__":
    main()
            

