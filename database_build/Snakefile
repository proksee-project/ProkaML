'''
Copyright:

University of Manitoba & National Microbiology Laboratory, Canada, 2021

Written by: Arnab Saha Mandal

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
'''

from datetime import date
import os
import glob

month_year_stamp = date.today().strftime("%b_%Y")
species_assembly_counts = 'species_assemblycounts_' + month_year_stamp + '.txt'
species_assembly_counts_taxonomy = 'species_assemblycounts_' + month_year_stamp + '_taxonomy.txt'
configfile: "config.yaml"
integrated_metadata_file = 'well_represented_species_metadata.txt'

rule target:
    input:
        integrated_metadata_file

rule species_count:
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

rule obtain_metadata:
    input:
        "id_list/{i}_idlist.txt"
    output:
        "species_metadata/{i}_metadata.txt"
    shell:
        "python get_genomic_metadata_filename.py {config[email]} {config[api_key]} {input}"

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.retrieve_id_list.get(**wildcards).output[0]
    return expand("species_metadata/{i}_metadata.txt", i=glob_wildcards(os.path.join(checkpoint_output, "{i}_idlist.txt")).i)

rule concatenate_metadata:
    input:
        aggregate_input
    output:
        integrated_metadata_file
    shell:
        "python concatenate_metadata.py"
