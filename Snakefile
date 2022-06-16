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

import os
import glob
from os import environ

configfile: "config.yaml"

rule target:
    input:
        "species_well_represented_assemblyQC.joblib.gz"


checkpoint get_assembly_UIDs:
    conda: 
        "dependencies.yaml"
    output:
        directory("database_build/temporary_outputs/entrez_id_list")
    shell:
        "python database_build/get_assembly_UIDs.py {config[email]} {config[api_key]}"


rule obtain_entrez_metadata:
    conda: 
        "dependencies.yaml"
    input:
        "database_build/temporary_outputs/entrez_id_list/Assembly_UID_chunk{i}.txt"
    output:
        "database_build/temporary_outputs/entrez_species_metadata/Assembly_UID_chunk{i}_metadata.txt"
    shell:
        "python database_build/get_entrez_metadata.py {config[email]} {config[api_key]} {input}"


def log_entrez_metadata_report(wildcards):
    ck_output = checkpoints.get_assembly_UIDs.get(**wildcards).output[0]
    return expand("database_build/temporary_outputs/log_entrez_metadata_chunk{i}.txt", 
        i=glob_wildcards(os.path.join(ck_output,"Assembly_UID_chunk{i}.txt")).i)


rule log_entrez_metadata:
    conda: 
        "dependencies.yaml"
    input:
        log_entrez_metadata_report
    shell:
        "python database_build/log_entrez_metadata.py"


def entrez_metadata_input(wildcards):
    ck_output = checkpoints.get_assembly_UIDs.get(**wildcards).output[0]
    return expand("database_build/temporary_outputs/entrez_species_metadata/Assembly_UID_chunk{i}_metadata.txt",
        i=glob_wildcards(os.path.join(ck_output,"Assembly_UID_chunk{i}.txt")).i)


checkpoint get_species_counts:
    conda: 
        "dependencies.yaml"
    input:
        entrez_metadata_input
    output:
        directory("database_build/temporary_outputs/species_reorganized_metadata")
    shell:
        "python database_build/get_species_counts.py {config[email]} {config[api_key]}"


rule append_gc_content:
    conda: 
        "dependencies.yaml"
    input:
        "database_build/temporary_outputs/species_reorganized_metadata/{j}_metadata.txt"
    output:
        "database_build/temporary_outputs/species_reorganized_metadata_gc/{j}_metadata_gc.txt"
    shell:
        "python database_build/append_gc_content.py {config[email]} {config[api_key]} {input}"


def log_gc_content(wildcards):
    ck_output = checkpoints.get_species_counts.get(**wildcards).output[0]
    return expand("database_build/{j}_log_gc.txt", j=glob_wildcards(os.path.join(ck_output,"{j}_metadata.txt")).j)


rule log_gc_content:
    conda: 
        "dependencies.yaml"
    input:
        log_gc_content
    shell:
        "python database_build/log_gc_content.py"


rule concatenate_metadata:
    conda: 
        "dependencies.yaml"
    input:
        log_gc_content
    output:
        "well_represented_species_metadata.txt"
    shell:
        "python database_build/concatenate_metadata.py"


rule preprocess_metadata:
    conda: 
        "dependencies.yaml"
    input:
        "well_represented_species_metadata.txt"
    output:
        "well_represented_species_metadata_normalized.txt"
    shell:
        "python preprocessing/cmd_preprocess_normalize.py"


rule generate_machine_learning_model:
    conda: 
        "dependencies.yaml"
    input:
        "well_represented_species_metadata_normalized.txt"
    output:
        "species_well_represented_assemblyQC.joblib.gz"
    shell:
        "python machine_learning/cmd_build_apply_model.py"
