from datetime import date
import os
import glob

month_year_stamp = date.today().strftime("%b_%Y")
species_assembly_counts = 'species_counts_' + month_year_stamp + '.txt'

configfile: "config.yaml"

CATEGORY = ['major', 'large', 'interm']

rule target:
    input:
        expand("{category}_species_metadata", category=CATEGORY)

rule species_count:
    output:
        species_assembly_counts      
    shell:
        "python assemblydb_entrez_query.py {config[email]} {config[api_key]}"

rule retrieve_id_list:
    input:
        species_assembly_counts
    output:
        directory(expand("id_list_{category}", category=CATEGORY))
    shell:
        "python idlist_retriever.py {config[email]} {config[api_key]} {input}"

rule obtain_metadata:
    input:
        "id_list_{category}"
    output:
        directory("{category}_species_metadata")
    run:
        for i in range(0, len(CATEGORY)):
            id_file_list = glob.glob('id_list_' + CATEGORY[i] + '/*_idlist.txt')
            for index in range(0, len(id_file_list)):
                shell("python metadata_print_fileindex.py {config[email]} {config[api_key]} id_list_{wildcards.category} {wildcards.category}_species_metadata {index}")