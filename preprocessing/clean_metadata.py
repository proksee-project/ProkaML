'''
Copyright:
University of Manitoba & National Microbiology Laboratory, Canada, 2020
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

'''
This program processes NCBI annotated metadata to calculate median of log of microbial specific
species specific genomic attributes
'''

import pandas as pd
import numpy as np
import os
import sys
import re 

class OrganizeMetadata():

	def __init__(self):
		pass

	def organize_assembly_method(self, dataframe):
		dataframe = dataframe.replace({'Assembly Method':r'.*[Aa]5.*$'},{'Assembly Method':'A5-miseq'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^A[Bb][Yy][Ss]{2}.*$'},{'Assembly Method':'ABYSS'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^[Aa][Ll]{2}[Pp][Aa][Tt][Hh][Ss].*$'},{'Assembly Method':'allpaths'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^AMO[Ss]cmp.*$'},{'Assembly Method':'AMOScmp'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'.*[Bb]owtie.*$'},{'Assembly Method':'Bowtie'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^CA.*$'},{'Assembly Method':'Celera'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Celera.*$'},{'Assembly Method':'Celera'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'.*[Cc][Aa][Nn][Uu].*'},{'Assembly Method':'Canu'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^[Cc][Ll][Cc].*$'},{'Assembly Method':'CLC'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'.*DNA[Ss][Tt][Aa][Rr].*$'},{'Assembly Method':'DNASTAR'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Seq[Mm]an.*$'},{'Assembly Method':'DNASTAR'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Ed.*?na.*$'},{'Assembly Method':'Edena'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'.*Falcon.*$'},{'Assembly Method':'Falcon'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Flye.*$'},{'Assembly Method':'Flye'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Geneious.*$'},{'Assembly Method':'Geneious'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^[Gg][Ss].*$'},{'Assembly Method':'GS'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Roche\sGS.*$'},{'Assembly Method':'GS'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'.*[Hh][Gg][Aa][Pp].*$'},{'Assembly Method':'HGAP'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'.*Hierarchical\sGenome.*$'},{'Assembly Method':'HGAP'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'.*[Ii][Dd][Bb][Aa].*$'},{'Assembly Method':'IDBA_UD'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^INNU[Cc][Aa].*$'},{'Assembly Method':'INNUca'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^MaSu[Rr]CA.*'},{'Assembly Method':'MaSuRCA'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^[Mm][Ee][Gg][Aa][Hh][Ii][Tt].*$'},{'Assembly Method':'Megahit'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^MetaBAT.*$'},{'Assembly Method':'MetaBAT'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^[Mm][Ii][Rr][Aa].*$'},{'Assembly Method':'MIRA'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^454\sNewbler.*$'},{'Assembly Method':'Newbler'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^[Nn]ewbler.*$'},{'Assembly Method':'Newbler'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Roche\sNewbler.*$'},{'Assembly Method':'Newbler'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^pilon.*$'},{'Assembly Method':'pilon'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^ProkaryoteAssembly.*$'},{'Assembly Method':'ProkaryoteAssembly'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Ray.*$'},{'Assembly Method':'Ray'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^[Ss]hovill.*$'},{'Assembly Method':'Shovill'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^[Ss][Kk][Ee][Ss][Aa].*$'},{'Assembly Method':'SKESA'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'.*SMRT.*$'},{'Assembly Method':'SMRT'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^S[Oo][Aa][Pp].*$'},{'Assembly Method':'SOAPdenovo'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'.*[Ss][Pp][Aa][Dd][Ee].*$'},{'Assembly Method':'SPAdes'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^[Uu]ni[Cc]y.*$'},{'Assembly Method':'unicycler'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'.*[Vv]el.+?t.*$'},{'Assembly Method':'Velvet'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^[Ww]gs.*$'},{'Assembly Method':'WgsAssembler'},regex=True)

		dataframe = dataframe.replace({'Assembly Method':r'^ARGO.*$'},{'Assembly Method':'ARGO'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^PATRIC.*$'},{'Assembly Method':'PATRIC'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Platanus.*$'},{'Assembly Method':'Platanus'},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^QUAST.*$'},{'Assembly Method':'QUAST'},regex=True)

		dataframe = dataframe.replace({'Assembly Method':r'^Allora.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Arachne.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Artemis.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^As\s.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^BOI.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^breseq.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^BWA.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^CAP3.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^CISA.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Clover.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^custom.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^De.*?[Nn]ovo.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^FALCON.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Galaxy.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^GeneStudio.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Illumina.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'.*in-house.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^IonTorrent.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Kbase.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Lasergene.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Mauve.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Microbe.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^[Mm]inia.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^MIX.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^MyPro.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^nanopolish.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^NextG[Ee][Nn]e.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^NovoAlign.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^other.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Parallel.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^PRJNA.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^prokka.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Shasta.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^SoftGenetics.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^SSPACE.*'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^TRIMMOMATIC.*'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Turing.*'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^Unknown.*$'},{'Assembly Method':np.nan},regex=True)
		dataframe = dataframe.replace({'Assembly Method':r'^VAAL.*$'},{'Assembly Method':np.nan},regex=True)
		
		return dataframe

	def organize_sequencing_technology(self, dataframe):
		dataframe = dataframe.replace({'Sequencing Technology':r'.*454.*$'},{'Sequencing Technology':'454'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^BGI.*$'},{'Sequencing Technology':'BGI'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^Complete\sGenomics.*$'},{'Sequencing Technology':'CompleteGenomics'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^DNB-Seq.*$'},{'Sequencing Technology':'BGI'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^[Ii].*?[Aa].*$'},{'Sequencing Technology':'Illumina'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^[HMhm][Ii][Ss][Ee][Qq].*$'},{'Sequencing Technology':'Illumina'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^[Nn]ext[Ss][Ee][Qq].*$'},{'Sequencing Technology':'Illumina'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^[Nn]ova[Ss][Ee][Qq].*$'},{'Sequencing Technology':'Illumina'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^Nextera.*$'},{'Sequencing Technology':'Illumina'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^llumina.*$'},{'Sequencing Technology':'Illumina'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^[Ii][Oo][Nn].*$'},{'Sequencing Technology':'IonTorrent'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^Oxford\sNanopore.*$'},{'Sequencing Technology':'OxfordNanopore'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'.*Min[Ii][Oo][Nn].*$'},{'Sequencing Technology':'OxfordNanopore'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^[Pp]ac[Bb]io.*$'},{'Sequencing Technology':'PacBio'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^Sanger.*$'},{'Sequencing Technology':'Sanger'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'^Solexa.*$'},{'Sequencing Technology':'Solexa'},regex=True)
		dataframe = dataframe.replace({'Sequencing Technology':r'.*[Ss][Oo][Ll]i[Dd].*$'},{'Sequencing Technology':'SOLiD'},regex=True)

		dataframe = dataframe.replace({'Sequencing Technology':r'.*Others.*$'},{'Sequencing Technology':np.nan},regex=True)

		return dataframe

	def subset_long_read(self, dataframe):
		long_read_index = dataframe[(dataframe['Sequencing Technology'] == 'OxfordNanopore') | \
                  					(dataframe['Sequencing Technology'] == 'PacBio') | \
                  					(dataframe['Assembly Method'] == 'Canu') | \
									(dataframe['Assembly Method'] == 'Falcon') | \
									(dataframe['Assembly Method'] == 'Flye') | \
									(dataframe['Assembly Method'] == 'HGAP') | \
									(dataframe['Assembly Method'] == 'SMRT') | \
									(dataframe['Assembly Method'] == 'pilon')].index

		return long_read_index
