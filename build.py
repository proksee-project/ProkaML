import sys
import pandas
import numpy
import csv

DELIMITER = "\t"

def calculate_staistics(dataframe):

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
        searchfor = ["Illumina", "454", "Ion Torrent", "IonTorrent", "SOLiD", "Sanger", "Solexa", "BGI", "CompleteGenomics"]
        dataframe = dataframe[dataframe['Sequencing Technology'].str.contains("|".join(searchfor), na=False, case=False)]

        # Check for all NA columns:
        bad_columns = dataframe.columns[dataframe.isna().all()].tolist()
        if len(bad_columns) > 0:
            #species_statistics.append([species])
            continue

        # Filter by only included in RefSeq
        dataframe = dataframe[dataframe['Refseq Accession'].str.contains("_", na=False, case=False)]

        dataframe = dataframe.iloc[:, [0, 2, 3, 4, 5]]

        # Check for no rows:
        if dataframe.count == 0:
            #species_statistics.append([species])
            continue

        statistics = calculate_staistics(dataframe)
        species_statistics.append(statistics)

    species_statistics = sorted(species_statistics, key=lambda x: x[0])

    with open("output.csv", "w") as output:
        output.write("species,count,n50_05,n50_20,n50_50,n50_80,n50_95,contig_count_05,contig_count_20,contig_count_50,contig_count_80,contig_count_95,l50_05,l50_20,l50_50,l50_80,l50_95,length_05,length_20,length_50,length_80,length_95\n")
        writer = csv.writer(output)
        writer.writerows(species_statistics)

