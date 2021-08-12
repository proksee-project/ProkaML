from datetime import date
import os
import glob

month_year_stamp = date.today().strftime("%b_%Y")
species_assembly_counts = 'species_assemblycounts_' + month_year_stamp + '.txt'
species_assembly_counts_taxonomy = 'species_assemblycounts_' + month_year_stamp + '_taxonomy.txt'
configfile: "config.yaml"

rule target:
    input:
        "species_metadata"

rule species_count:
    #conda: config["conda_environment"]
    output:
        species_assembly_counts      
    shell:
        "python get_species_assemblies.py {config[email]} {config[api_key]}"

rule append_taxonomical_kingdom:
    input:
        species_assembly_counts
    output:
        species_assembly_counts_taxonomy
    shell:
        "python get_taxonomical_kingdom.py {config[email]} {config[api_key]}"

checkpoint retrieve_id_list:
    input:
        species_assembly_counts_taxonomy
    output:
        directory("id_list")
    shell:
        "python retrieve_idlist.py {config[email]} {config[api_key]} {input}"

def get_idlist_filepath(wildcards):
    checkpoint_output = checkpoints.retrieve_id_list.get(**wildcards).output[0]
    return expand("id_list/{i}.txt", i=glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)

rule obtain_metadata:
    input:
        get_idlist_filepath
    output:
        directory("species_metadata")
    run:
        for file in input:
            command = "python get_genomic_metadata_filename.py {config[email]} {config[api_key]} {file}"
            shell(command)
