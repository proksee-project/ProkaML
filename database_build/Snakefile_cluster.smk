from datetime import date
import os
import glob

month_year_stamp = date.today().strftime("%b_%Y")
species_assembly_counts = 'species_counts_' + month_year_stamp + '.txt'

configfile: "config.yaml"

CATEGORY = ['interm', 'large', 'major']

rule target:
    input:
        "command.sh"

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
        "python idlist_retriever_categorical.py {config[email]} {config[api_key]} {input}"

rule cluster_script:
    input:
        expand("id_list_{category}", category=CATEGORY)
    output:
        "command.sh"
    shell:
        "python create_cluster_jobscripts.py {config[email]} {config[api_key]} {config[output_directory]} {config[partition]}"
