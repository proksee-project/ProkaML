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

    statistics.append(numpy.percentile(n50_list, 5))
    statistics.append(numpy.percentile(n50_list, 20))
    statistics.append(numpy.percentile(n50_list, 80))
    statistics.append(numpy.percentile(n50_list, 95))

    statistics.append(numpy.percentile(num_contigs_list, 5))
    statistics.append(numpy.percentile(num_contigs_list, 20))
    statistics.append(numpy.percentile(num_contigs_list, 80))
    statistics.append(numpy.percentile(num_contigs_list, 95))

    statistics.append(numpy.percentile(l50_list, 5))
    statistics.append(numpy.percentile(l50_list, 20))
    statistics.append(numpy.percentile(l50_list, 80))
    statistics.append(numpy.percentile(l50_list, 95))

    statistics.append(numpy.percentile(length_list, 5))
    statistics.append(numpy.percentile(length_list, 20))
    statistics.append(numpy.percentile(length_list, 80))
    statistics.append(numpy.percentile(length_list, 95))

    statistics = [statistics[0]] + ['%.2f' % elem for elem in statistics[1:]]

    return statistics

if __name__ == '__main__':

    SPECIES = 0
    N50 = 10
    NUM_CONTIGS = 11
    L50 = 12
    LENGTH = 13

    ASSEMBLY_METHOD = 14
    TECHNOLOGY = 15

    COLS = [SPECIES, N50, NUM_CONTIGS, L50, LENGTH, ASSEMBLY_METHOD, TECHNOLOGY]

    species_statistics = []

    for file_location in sys.argv[1:]:

        print(file_location)
        dataframe = pandas.read_csv(file_location, delimiter=DELIMITER)
        dataframe = dataframe.iloc[:, COLS]

        searchfor = ["Illumina", "454", "Ion Torrent", "IonTorrent"]
        dataframe = dataframe[dataframe['Sequencing Technology'].str.contains("|".join(searchfor), na=False, case=False)]
        dataframe = dataframe.iloc[:, :5]

        statistics = calculate_staistics(dataframe)
        species_statistics.append(statistics)

    species_statistics = sorted(species_statistics, key=lambda x: x[0])

    with open("output.csv", "w") as output:
        output.write("species,n50_05,n50_20,n50_80,n50_95,contig_count_05,contig_count_20,contig_count_80,contig_count_95,l50_05,l50_20,l50_80,l50_95,length_05,length_20,length_80,length_95\n")
        writer = csv.writer(output)
        writer.writerows(species_statistics)

