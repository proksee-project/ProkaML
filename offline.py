import sys

mash_info_filename = sys.argv[1]
ranked_lineage_filename = sys.argv[2]
merged_filename = sys.argv[3]

accessions = {}

with open(mash_info_filename) as f:

    next(f)

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

for filename in sys.argv[4:]:
    with open(filename) as f:

        next(f) # skip header

        for line in f:

            tokens = line.split()

            accession = tokens[0]
            taxid = tokens[2]

            if accession in accessions:

                accessions[accession] = taxid


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

merged_ids = {}

with open(merged_filename) as f:

    for line in f:

        tokens = line.split("|")

        old_id = tokens[0].strip()
        new_id = tokens[1].strip()

        merged_ids[old_id] = new_id

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

