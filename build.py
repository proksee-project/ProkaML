"""
Copyright Government of Canada 2023

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

import sys
import pandas
import numpy
import csv

DELIMITER = "\t"


def calculate_staistics(dataframe):
    """
    Calculates the summary statistics for the dataframe. This includes various percentiles for the N50, number of
    contigs, L50, and assembly length. Basically, it takes assembly statistics for many assemblies of the same species
    and summarizes it one row, showing the calculated percentiles.

    PARAMETERS:
        dataframe (pandas.DataFrame): dataframe where each row represents an assembly and each column represents an
                                      assembly summary statistic (name, N50, number of contigs, etc.).

    RETURNS:
        statistics (string[]): a summarized report of the percentiles for the entire dataframe, in the following
                               format:

                               [species, counts,
                               n50 5%tile, n50 20%tile, n50 50%tile, n50 80%tile, n50 95%tile,
                               contigs 5%tile, contigs 20%tile, contigs 50%tile, contigs 80%tile, contigs 95%tile,
                               l50 5%tile, l50 20%tile, l50 50%tile, l50 80%tile, l50 95%tile,
                               length 5%tile, length 20%tile, length 50%tile, length 80%tile, length 95%tile]
    """

    SPECIES = "Organism Name"
    N50 = "ContigN50"
    NUM_CONTIGS = "Contig count"
    L50 = "ContigL50"
    LENGTH = "Total length"

    statistics = []

    species = dataframe[SPECIES].tolist()[0]
    n50_list = dataframe[N50].tolist()
    num_contigs_list = dataframe[NUM_CONTIGS].tolist()
    l50_list = dataframe[L50].tolist()
    length_list = dataframe[LENGTH].tolist()

    statistics.append(species)

    # Observed counts for the given species:
    statistics.append(len(dataframe[SPECIES].tolist()))

    statistics.append(numpy.percentile(n50_list, 5))
    statistics.append(numpy.percentile(n50_list, 20))
    statistics.append(numpy.percentile(n50_list, 50))
    statistics.append(numpy.percentile(n50_list, 80))
    statistics.append(numpy.percentile(n50_list, 95))

    statistics.append(numpy.percentile(num_contigs_list, 5))
    statistics.append(numpy.percentile(num_contigs_list, 20))
    statistics.append(numpy.percentile(num_contigs_list, 50))
    statistics.append(numpy.percentile(num_contigs_list, 80))
    statistics.append(numpy.percentile(num_contigs_list, 95))

    statistics.append(numpy.percentile(l50_list, 5))
    statistics.append(numpy.percentile(l50_list, 20))
    statistics.append(numpy.percentile(l50_list, 50))
    statistics.append(numpy.percentile(l50_list, 80))
    statistics.append(numpy.percentile(l50_list, 95))

    statistics.append(numpy.percentile(length_list, 5))
    statistics.append(numpy.percentile(length_list, 20))
    statistics.append(numpy.percentile(length_list, 50))
    statistics.append(numpy.percentile(length_list, 80))
    statistics.append(numpy.percentile(length_list, 95))

    statistics = [statistics[0]] + [statistics[1]] + ['%.2f' % elem for elem in statistics[2:]]

    return statistics


if __name__ == '__main__':
    """
    This script summarizes assembly statistics by species into percentiles.

    usage:
        python build.py [assembly statistics metadata file, ...]

    examples:
        python build.py Staphylococcus_aureus_metadata.txt Campylobacter_coli_metadata.txt
        python build.py major_species_metadata/*

    The script takes as arguments multiple assembly statistics metadata files, each of which
    is formatted in the following way:

    assembly metadata file column format:

    Organism Name, Organism Detail, Strain/Isolate, Assembly Name, Genbank Accession, Refseq Accession, \
    Genome Coverage, Submission Date, Last Update Date, Refseq Exclusion Reason, ContigN50, Contig count, \
    ContigL50, Total length, Assembly Method, Sequencing Technology

    The first row is expected to be a header containing the above column headers.

    The script will write out a file called "assembly_metadata_summary.csv" containing the summarized
    overview of all passed files. This file will automatically be overwritten if it already exists.
    """

    SPECIES = 0
    REFSEQ_ACCESSION = 5
    N50 = 10
    NUM_CONTIGS = 11
    L50 = 12
    LENGTH = 13

    ASSEMBLY_METHOD = 14
    TECHNOLOGY = 15

    COLS = [SPECIES, REFSEQ_ACCESSION, N50, NUM_CONTIGS, L50, LENGTH, ASSEMBLY_METHOD, TECHNOLOGY]

    species_statistics = []

    for file_location in sys.argv[1:]:

        print(file_location)
        dataframe = pandas.read_csv(file_location, delimiter=DELIMITER)
        dataframe = dataframe.iloc[:, COLS]
        species = dataframe["Organism Name"].tolist()[0]

        # Filter for sequencing platform
        searchfor = ["Illumina"]
        # searchfor = ["Illumina", "454", "Ion Torrent", "IonTorrent",
        #             "SOLiD", "Sanger", "Solexa", "BGI", "CompleteGenomics"]
        dataframe = dataframe[dataframe['Sequencing Technology'].str.contains(
                              "|".join(searchfor), na=False, case=False)]

        # Check for all NA columns:
        # i.e. If there is no data for an entire column of this particular species,
        # we won't be able to calculate the percentile for it, so skip the species
        # and go to the next.
        bad_columns = dataframe.columns[dataframe.isna().all()].tolist()
        if len(bad_columns) > 0:
            continue

        # Filter by only included in RefSeq
        dataframe = dataframe[dataframe['Refseq Accession'].str.contains("_", na=False, case=False)]

        dataframe = dataframe.iloc[:, [0, 2, 3, 4, 5]]

        # Check for no rows:
        # i.e. If there are no rows for the species, then there is no data (only a header row),
        # and we won't be able to calculate anything, so skip the species and go to the next.
        if dataframe.count == 0:
            continue

        statistics = calculate_staistics(dataframe)
        species_statistics.append(statistics)

    species_statistics = sorted(species_statistics, key=lambda x: x[0])

    with open("assembly_metadata_summary.csv", "w") as output:
        output.write("species,count,n50_05,n50_20,n50_50,n50_80,n50_95," +
                     "contig_count_05,contig_count_20,contig_count_50,contig_count_80,contig_count_95," +
                     "l50_05,l50_20,l50_50,l50_80,l50_95,length_05,length_20,length_50,length_80,length_95\n")
        writer = csv.writer(output)
        writer.writerows(species_statistics)
