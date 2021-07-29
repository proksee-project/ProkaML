from datetime import date
import os
import glob

month_year_stamp = date.today().strftime("%b_%Y")
species_assembly_counts = 'species_counts_' + month_year_stamp + '.txt'

configfile: "config.yaml"

rule target:
    input:
        "species_metadata"

rule species_count:
    output:
        species_assembly_counts      
    shell:
        "python assemblydb_entrez_query.py {config[email]} {config[api_key]}"

checkpoint retrieve_id_list:
    input:
        species_assembly_counts
    output:
        directory("id_list")
    shell:
        "python idlist_retriever.py {config[email]} {config[api_key]} {input}"

def metadata_input(wildcards):
    checkpoint_output = checkpoints.retrieve_id_list.get(**wildcards).output[0]
    return expand("id_list/{i}.txt",
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)

rule obtain_metadata:
    input:
        metadata_input
    output:
        directory("species_metadata")
    shell:
        "python metadata_print_filename.py {config[email]} {config[api_key]} {wildcards.input} {output}"
