"""
Copyright Government of Canada 2021

Written by:

Eric Marinier
    National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this work except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
"""

"""
This script builds a tab-separated mapping file to be used with ONE corresponding Mash sketch file. The script
maps NCBI accession IDs found in the 'query comment' field of the output of 'mash info' (of a Mash sketch file)
to their NCBI taxonomy ID, their ranked lineage, and their full lineage.

The mapping file will have a format similar to the following:

[Accession ID]    [Taxonomy ID]    [Ranked Lineage (multiple columns)]    [Full Name Lineage (one column)]
NC_005100    10116    Rattus norvegicus    Rattus   Muridae Rodentia   [...]   Organisms; Eukaryota; [...]

The following files will be required:

mash_info: the output from running 'mash info' on a Mash sketch; the file produced by this program is
    expected to work ONLY with this Mash sketch file
ranked_lineage: a NCBI file mapping NCBI taxonomy IDs to their 'ranked' lineage
full_name_lineage: a NCBI file mapping NCBI taxonomy IDs to their 'full' lineage
merged: a NCBI file mapping old NCBI accession IDs to their new, merged NCBI accession IDs
(multiple) mapping files: a list of NCBI files mapping NCBI accession IDs to NCBI taxonomy IDs

As of 2021-06-23, the following files were used when building a mapping file:

mash_info: output of 'mash info -t' on this sketch: https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh

ranked_lineage: rankedlineage.dmp in https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

full_name_lineage: fullnamelineage.dmp in https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

merged: merged.dmp in https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

mapping files:
    - https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_wgs.accession2taxid.gz
    - https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz
    - https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    - https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
"""

import sys
import argparse

MASH_INFO_FILENAME = "mash_info_filename"
RANKED_LINEAGE_FILENAME = "ranked_lineage_filename"
FULL_NAME_LINEAGE_FILENAME = "full_name_lineage_filename"
MERGED_FILENAME = "merged_filename"
MAPPING_FILENAMES = "mapping_filenames"


def build_mapping_file():
    """
    Builds a tab-separated mapping file. The file maps NCBI accession IDs to its associated tax ID
    and taxonomy lineage information.

    POST
        If successful, tab-separated information mapping:

        [(Old/New) NCBI Accession ID] -> [TaxID, Ranked Lineage, Full Name Lineage]

        will be written to standard output. For example:

        [Accession ID]    [Taxonomy ID]    [Ranked Lineage]    [Full Name Lineage]
        NC_005100    10116    Rattus norvegicus    Rattus   Muridae Rodentia   [...]   Organisms; Eukaryota; [...]
    """

    # Argument parsing:
    args = parse_arguments()
    parameters = vars(args)

    mash_info_filename = parameters.get(MASH_INFO_FILENAME)
    ranked_lineage_filename = parameters.get(RANKED_LINEAGE_FILENAME)
    full_name_lineage_filename = parameters.get(FULL_NAME_LINEAGE_FILENAME)
    merged_filename = parameters.get(MERGED_FILENAME)
    mapping_filenames = parameters.get(MAPPING_FILENAMES)

    # Build dictionary mapping [NCBI Accession ID] -> [NCBI Taxonomy ID]
    accessions = get_accession_IDs(mash_info_filename)
    add_taxonomy_to_accessions(accessions, mapping_filenames)

    # Build dictionary mapping [NCBI Taxnomy ID] -> [Ranked Lineage, Full Name Lineage]
    lineages = map_taxonomy_to_lineage(ranked_lineage_filename)
    append_full_name_lineage(lineages, full_name_lineage_filename)

    # Build dictionary mapping [Old NCBI Accession ID] -> [New NCBI Accession ID]
    # This is because some accession IDs are no longer used.
    merged_ids = build_merged_IDs(merged_filename)

    # Output the file mapping [(Old/New) NCBI Accession ID] -> [TaxID, Ranked Lineage, Full Name Lineage]
    output_mapping(accessions, lineages, merged_ids)


def parse_arguments():
    """
    Parses the command-line arguments are returns a populated Namespace object with argument strings as
    attributes.

    RETURNS:
        args (argparse.Namespace): an Namespace object with argument strings as attributes
    """

    parser = argparse.ArgumentParser(description='Builds an NCBI ID-to-taxonomy mapping file.')

    parser.add_argument('-mi', dest=MASH_INFO_FILENAME, type=str, required=True,
                        help='The output of "mash info" for the Mash sketch file.')
    parser.add_argument('-r', dest=RANKED_LINEAGE_FILENAME, type=str, required=True,
                        help='NCBI ranked lineage file mapping IDs to lineages.')
    parser.add_argument('-f', dest=FULL_NAME_LINEAGE_FILENAME, type=str, required=True,
                        help='NCBI full name lineage file mapping IDs to full name lineages.')
    parser.add_argument('-m', dest=MERGED_FILENAME, type=str, required=True,
                        help='NCBI file mapping old IDs to new, merged IDs.')
    parser.add_argument('-maps', dest=MAPPING_FILENAMES, type=str, required=True, nargs='+',
                        help='NCBI files mapping Accession IDs to taxonomy IDs.')

    args = parser.parse_args()
    
    return args

def get_accession_IDs(mash_info_filename):
    """
    Get all NCBI accession IDs that are observed in the output of 'mash info'. The format of the 'mash info' file
    needs to be very particular:

    1000    143726002    GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz    [1870 seqs] NC_004354.4 Drosophila [...]

    That is, the accession needs to appear in the 'query comment' field of 'mash info' output. It may optionally be preceeded
    by a "[# seqs]" substring.

    PARAMETERS:
        mash_info_filename (str): the name of the file containing the output of 'mash_info'

    RETURNS:
        accessions (dict(str->None)): a dictionary containing only NCBI accession IDs as keys, all mapping to None
    """

    accessions = {}

    with open(mash_info_filename) as f:

        next(f) # Skip header

        line = next(f, None)

        while line:
            line = line.strip()

            tokens = line.split("\t")
            string = tokens[3]

            # Check for optional "[# seqs]":
            if string.startswith("["):
                string = string.split("] ")[1]

            # Grab the Accession ID:
            accession = string.split(" ")[0]
            accession = accession.split(".")[0] # Remove the version (NC_000001.1 -> NC_000001)
            accessions[accession] = None

            line = next(f, None)

    return accessions

def add_taxonomy_to_accessions(accessions, mapping_filenames):
    """
    Adds NCBI taxonomy IDs to a NCBI accession IDs dictionary by iterating over several files mapping accession IDs to
    taxonomy IDs and adding matches when found.

    PARAMETERS:
        accessions (dict(str->None)): a dictionary with NCBI accession IDs as keys; any values may be overwritten. This
            dictionary will be editted.

    POST
        The 'accessions' dictionary will be modified to have accessions added as values for each taxonomy ID key.
    """

    for filename in mapping_filenames:
        with open(filename) as f:

            next(f) # skip header

            for line in f:

                tokens = line.split()
                accession = tokens[0]
                taxid = tokens[2]

                if accession in accessions:

                    accessions[accession] = taxid

def map_taxonomy_to_lineage(ranked_lineage_filename):
    """
    Maps taxnomy IDs to ranked lineage information. In particular, builds a dictionary of the following form:

    [Taxonomy ID] -> [Taxonomy Name, Species, Genus, Family, Order, Class, Phylum, Kingdom, Superkingdom]

    PARAMETERS:
        ranked_lineage_filename (str): the name of the file mapping NCBI taxonomy IDs to their full ranked lineage

    RETURNS:
        lineages (dict(int->list(str))): a dictionary mapping taxonmy IDs to ranked lineage information
    """

    lineages = {}

    with open(ranked_lineage_filename) as f:

        for line in f:

            tokens = line.split("|")

            taxid = tokens[0].strip()

            tax_name = tokens[1].strip() if len(tokens[1].strip()) > 0 else "-"
            species = tokens[2].strip() if len(tokens[2].strip()) > 0 else "-"
            genus = tokens[3].strip() if len(tokens[3].strip()) > 0 else "-"
            family = tokens[4].strip() if len(tokens[4].strip()) > 0 else "-"
            order = tokens[5].strip() if len(tokens[5].strip()) > 0 else "-"
            clas = tokens[6].strip() if len(tokens[6].strip()) > 0 else "-"
            phylum = tokens[7].strip() if len(tokens[7].strip()) > 0 else "-"
            kingdom = tokens[8].strip() if len(tokens[8].strip()) > 0 else "-"
            superkingdom = tokens[9].strip() if len(tokens[9].strip()) > 0 else "-"

            lineages[taxid] = [tax_name, species, genus, family, order, clas, phylum, kingdom, superkingdom]

    return lineages

def append_full_name_lineage(lineages, full_name_lineage_filename):
    """
    Appends the 'full name lineage' to the lineage dictionary mapping taxonomy IDs to ranked lineage information. The
    dictionary will have the following form after complection:

    [Taxonomy ID] -> [Taxonomy Name, Species, Genus, Family, Order, Class, Phylum, Kingdom, Superkingdom, Full Name]

    The 'full name lineage' tends to look like the following:

    "cellular organisms; Eukaryota; Opisthokonta; Metazoa; [...]; Monotremata; Ornithorhynchidae; Ornithorhynchus"

    PARAMETERS:
        lineages (dict(int->list(str))): a dictionary mapping taxonmy IDs to ranked lineage information; this dictionary
            will be modified after execution
        full_name_lineage_filename (str): the name of the file containing taxonomy ID to full name lineage mapping
            information

    POST
        The 'lineages' parameter will have full name lineages appended to the lineage information for each taxonomy ID.
    """
    
    with open(full_name_lineage_filename) as f:

        for line in f:

            tokens = line.split("|")

            taxid = tokens[0].strip()

            full_lineage = tokens[2].strip()
            full_lineage = full_lineage[:-1]  # remove the trailing ";" character

            lineages[taxid].append(full_lineage)


def build_merged_IDs(merged_filename):
    """
    Builds a dictionary mapping old, deprecated NCBI accession IDs to new NCBI accession IDs. This is (usually?) because
    the old IDs have been merged together into a new ID.

    PARAMETERS:
        merged_filename (str): the name of the file mapping old NCBI accession IDs to new NCBI accession IDs

    RETURN:
        merged_ids (dict(str->str)): a dictionary mapping old IDs to new IDs
    """

    merged_ids = {}

    with open(merged_filename) as f:

        for line in f:

            tokens = line.split("|")

            old_id = tokens[0].strip()
            new_id = tokens[1].strip()

            merged_ids[old_id] = new_id

    return merged_ids

def output_mapping(accessions, lineages, merged_ids):
    """
    Iterates through various dictionaries to output, line-by-line, and writes the mapping of NCBI accession ID to
    its associated taxnomic lineage information.

    PARAMETERS:
        accessions (dict(str->str)): dictionary mapping NCBI accession IDs to NCBI taxonomy IDs
        lineages (dict(int->list(str))): dictionary mapping taxonomy IDs to ranked and full lineage information
        merged_ids (dict(str->str)): dictionary mapping old NCBI accession IDs to new NCBI accession IDs

    POST
        If successful, tab-separated information mapping:

        [(Old/New) NCBI Accession ID] -> [TaxID, Ranked Lineage, Full Name Lineage]

        will be written to standard output. For example:

        [Accession ID]    [Taxonomy ID]    [Ranked Lineage]    [Full Name Lineage]
        NC_005100    10116    Rattus norvegicus    Rattus   Muridae Rodentia   [...]   Organisms; Eukaryota; [...]
    """

    for accession in accessions:

        taxid = accessions[accession]

        if taxid in lineages:
            lineage = lineages[taxid]

            output = ""
            output += accession + "\t"
            output += taxid + "\t"
            output += "\t".join(lineage)

        elif taxid in merged_ids:

            new_id = merged_ids[taxid]
            lineage = lineages[new_id]

            output = ""
            output += accession + "\t"
            output += taxid + "\t"
            output += "\t".join(lineage)

        else:
            output = "FAILURE: " + str(taxid)

        print(output)

# Run the program:
build_mapping_file()

