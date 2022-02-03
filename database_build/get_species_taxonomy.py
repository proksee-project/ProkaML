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

import pandas as pd
from Bio import Entrez
from datetime import date
import argparse
import re
import constants as const


class Taxonomy():

    def __init__(self, species, email, api_key):
        self.species = species
        self.email = email
        self.api_key = api_key

    def get_taxonomy_record(self):

        Entrez.email = self.email
        Entrez.api_key = self.api_key

        for attempts in range(1, 4):
            try:
                handle = Entrez.esearch(db=const.TAXONOMY_DATABASE, term=self.species)

            except Exception:
                record = {}

            else:
                record = Entrez.read(handle)
                if record:
                    break

        return record


    def get_taxonomy_dict(self, record):

        taxonomy_dict = {}
        if record and len(record[const.RECORD_IDLIST]) > 0:
            for attempts in range(1, 4):        
                try:
                    efetch = Entrez.efetch(db=const.TAXONOMY_DATABASE, 
                                        id=record[const.RECORD_IDLIST],
                                        retmode=const.TAXONOMY_RETMODE)

                except Exception:
                    pass

                else:
                    taxon_summary = Entrez.read(efetch)
                    for i in range(len(taxon_summary[0][const.LINEAGE])):
                        taxonomy_dict[taxon_summary[0][const.LINEAGE][i][const.RANK]] = \
                            taxon_summary[0][const.LINEAGE][i][const.SCIENTIFIC_NAME]

                    break

        return taxonomy_dict


    def get_taxonomy_kingdom(self, taxonomy_dict):

        if taxonomy_dict and const.SUPERKINGDOM in taxonomy_dict:
            kingdom = taxonomy_dict[const.SUPERKINGDOM]

        else:
            kingdom = const.NA

        return kingdom


    def get_taxonomy_phylum(self, taxonomy_dict):

        if taxonomy_dict and const.PHYLUM in taxonomy_dict:
            phylum = taxonomy_dict[const.PHYLUM]

        else:
            phylum = const.NA

        return phylum


    def get_taxonomy_class(self, taxonomy_dict):

        if taxonomy_dict and const.CLASS in taxonomy_dict:
            taxonomy_class = taxonomy_dict[const.CLASS]

        else:
            taxonomy_class = const.NA

        return taxonomy_class


    def get_taxonomy_order(self, taxonomy_dict):

        if taxonomy_dict and const.ORDER in taxonomy_dict:
            order = taxonomy_dict[const.ORDER]

        else:
            order = const.NA

        return order


    def get_taxonomy_family(self, taxonomy_dict):

        if taxonomy_dict and const.FAMILY in taxonomy_dict:
            family = taxonomy_dict[const.FAMILY]

        else:
            family = const.NA

        return family


    def get_taxonomy_genus(self, taxonomy_dict):

        if taxonomy_dict and const.GENUS in taxonomy_dict:
            genus = taxonomy_dict[const.GENUS]

        else:
            genus = const.NA

        return genus


    def get_full_taxonomy(self):
        taxonomy_record = self.get_taxonomy_record()
        taxonomy_dict = self.get_taxonomy_dict(taxonomy_record)
        kingdom = self.get_taxonomy_kingdom(taxonomy_dict)
        phylum = self.get_taxonomy_phylum(taxonomy_dict)
        taxonomy_class = self.get_taxonomy_class(taxonomy_dict)
        order = self.get_taxonomy_order(taxonomy_dict)
        family = self.get_taxonomy_family(taxonomy_dict)
        genus = self.get_taxonomy_genus(taxonomy_dict)

        return [kingdom, phylum, taxonomy_class, order, family, genus]
