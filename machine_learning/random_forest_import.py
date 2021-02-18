import pandas as pd
import numpy as np
import re
import os
from pathlib import Path 
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import confusion_matrix
import joblib

rfc = RandomForestClassifier()

START_DIR = Path(__file__).resolve().parents[1]
filepath = '{}/gc_content_calculate/major_species_normalized.txt'.format(str(START_DIR))

def assembly_method_categorical(x):    
    x = x.replace({'Assembly Method':r'.*[Aa]5.*$'},{'Assembly Method':'A5-miseq'},regex=True)
    x = x.replace({'Assembly Method':r'^A[Bb][Yy][Ss]{2}.*$'},{'Assembly Method':'ABYSS'},regex=True)
    x = x.replace({'Assembly Method':r'^[Aa][Ll]{2}[Pp][Aa][Tt][Hh][Ss].*$'},{'Assembly Method':'allpaths'},regex=True)
    x = x.replace({'Assembly Method':r'^AMO[Ss]cmp.*$'},{'Assembly Method':'AMOScmp'},regex=True)
    x = x.replace({'Assembly Method':r'.*[Bb]owtie.*$'},{'Assembly Method':'Bowtie'},regex=True)
    x = x.replace({'Assembly Method':r'^CA.*$'},{'Assembly Method':'Celera'},regex=True)
    x = x.replace({'Assembly Method':r'^Celera.*$'},{'Assembly Method':'Celera'},regex=True)
    x = x.replace({'Assembly Method':r'.*[Cc][Aa][Nn][Uu].*'},{'Assembly Method':'Canu'},regex=True)
    x = x.replace({'Assembly Method':r'^[Cc][Ll][Cc].*$'},{'Assembly Method':'CLC'},regex=True)
    x = x.replace({'Assembly Method':r'.*DNA[Ss][Tt][Aa][Rr].*$'},{'Assembly Method':'DNASTAR'},regex=True)
    x = x.replace({'Assembly Method':r'^Seq[Mm]an.*$'},{'Assembly Method':'DNASTAR'},regex=True)
    x = x.replace({'Assembly Method':r'^Ed.*?na.*$'},{'Assembly Method':'Edena'},regex=True)
    x = x.replace({'Assembly Method':r'^Flye.*$'},{'Assembly Method':'Flye'},regex=True)
    x = x.replace({'Assembly Method':r'^Geneious.*$'},{'Assembly Method':'Geneious'},regex=True)
    x = x.replace({'Assembly Method':r'^[Gg][Ss].*$'},{'Assembly Method':'GS'},regex=True)
    x = x.replace({'Assembly Method':r'^Roche\sGS.*$'},{'Assembly Method':'GS'},regex=True)
    x = x.replace({'Assembly Method':r'.*[Hh][Gg][Aa][Pp].*$'},{'Assembly Method':'HGAP'},regex=True)
    x = x.replace({'Assembly Method':r'.*Hierarchical\sGenome.*$'},{'Assembly Method':'HGAP'},regex=True)
    x = x.replace({'Assembly Method':r'.*[Ii][Dd][Bb][Aa].*$'},{'Assembly Method':'IDBA_UD'},regex=True)
    x = x.replace({'Assembly Method':r'^INNU[Cc][Aa].*$'},{'Assembly Method':'INNUca'},regex=True)
    x = x.replace({'Assembly Method':r'^MaSu[Rr]CA.*'},{'Assembly Method':'MaSuRCA'},regex=True)
    x = x.replace({'Assembly Method':r'^[Mm][Ee][Gg][Aa][Hh][Ii][Tt].*$'},{'Assembly Method':'Megahit'},regex=True)
    x = x.replace({'Assembly Method':r'^MetaBAT.*$'},{'Assembly Method':'MetaBAT'},regex=True)
    x = x.replace({'Assembly Method':r'^[Mm][Ii][Rr][Aa].*$'},{'Assembly Method':'MIRA'},regex=True)
    x = x.replace({'Assembly Method':r'^454\sNewbler.*$'},{'Assembly Method':'Newbler'},regex=True)
    x = x.replace({'Assembly Method':r'^[Nn]ewbler.*$'},{'Assembly Method':'Newbler'},regex=True)
    x = x.replace({'Assembly Method':r'^Roche\sNewbler.*$'},{'Assembly Method':'Newbler'},regex=True)
    x = x.replace({'Assembly Method':r'^pilon.*$'},{'Assembly Method':'pilon'},regex=True)
    x = x.replace({'Assembly Method':r'^ProkaryoteAssembly.*$'},{'Assembly Method':'ProkaryoteAssembly'},regex=True)
    x = x.replace({'Assembly Method':r'^Ray.*$'},{'Assembly Method':'Ray'},regex=True)
    x = x.replace({'Assembly Method':r'^[Ss]hovill.*$'},{'Assembly Method':'Shovill'},regex=True)
    x = x.replace({'Assembly Method':r'^[Ss][Kk][Ee][Ss][Aa].*$'},{'Assembly Method':'SKESA'},regex=True)
    x = x.replace({'Assembly Method':r'.*SMRT.*$'},{'Assembly Method':'SMRT'},regex=True)
    x = x.replace({'Assembly Method':r'^S[Oo][Aa][Pp].*$'},{'Assembly Method':'SOAPdenovo'},regex=True)
    x = x.replace({'Assembly Method':r'.*[Ss][Pp][Aa][Dd][Ee].*$'},{'Assembly Method':'SPAdes'},regex=True)
    x = x.replace({'Assembly Method':r'^[Uu]ni[Cc]y.*$'},{'Assembly Method':'unicycler'},regex=True)
    x = x.replace({'Assembly Method':r'.*[Vv]el.+?t.*$'},{'Assembly Method':'Velvet'},regex=True)
    x = x.replace({'Assembly Method':r'^[Ww]gs.*$'},{'Assembly Method':'WgsAssembler'},regex=True)

    x = x.replace({'Assembly Method':r'^ARGO.*$'},{'Assembly Method':'ARGO'},regex=True)
    x = x.replace({'Assembly Method':r'^PATRIC.*$'},{'Assembly Method':'PATRIC'},regex=True)
    x = x.replace({'Assembly Method':r'^Platanus.*$'},{'Assembly Method':'Platanus'},regex=True)
    x = x.replace({'Assembly Method':r'^QUAST.*$'},{'Assembly Method':'QUAST'},regex=True)

    x = x.replace({'Assembly Method':r'^Allora.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Arachne.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Artemis.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^As\s.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^BOI.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^breseq.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^BWA.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^CAP3.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^CISA.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Clover.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^custom.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^De.*?[Nn]ovo.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^FALCON.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Galaxy.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^GeneStudio.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Illumina.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'.*in-house.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^IonTorrent.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Kbase.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Lasergene.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Mauve.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Microbe.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^[Mm]inia.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^MIX.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^MyPro.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^nanopolish.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^NextG[Ee][Nn]e.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^NovoAlign.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^other.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Parallel.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^PRJNA.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^prokka.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Shasta.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^SoftGenetics.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^SSPACE.*'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^TRIMMOMATIC.*'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Turing.*'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^Unknown.*$'},{'Assembly Method':np.nan},regex=True)
    x = x.replace({'Assembly Method':r'^VAAL.*$'},{'Assembly Method':np.nan},regex=True)
    
    return x

def seq_plat_categorical(x):
    x = x.replace({'Sequencing Technology':r'.*454.*$'},{'Sequencing Technology':'454'},regex=True)
    x = x.replace({'Sequencing Technology':r'^BGI.*$'},{'Sequencing Technology':'BGI'},regex=True)
    x = x.replace({'Sequencing Technology':r'^Complete\sGenomics.*$'},{'Sequencing Technology':'CompleteGenomics'},regex=True)
    x = x.replace({'Sequencing Technology':r'^DNB-Seq.*$'},{'Sequencing Technology':'BGI'},regex=True)
    x = x.replace({'Sequencing Technology':r'^[Ii].*?[Aa].*$'},{'Sequencing Technology':'Illumina'},regex=True)
    x = x.replace({'Sequencing Technology':r'^[HMhm][Ii][Ss][Ee][Qq].*$'},{'Sequencing Technology':'Illumina'},regex=True)
    x = x.replace({'Sequencing Technology':r'^[Nn]ext[Ss][Ee][Qq].*$'},{'Sequencing Technology':'Illumina'},regex=True)
    x = x.replace({'Sequencing Technology':r'^[Nn]ova[Ss][Ee][Qq].*$'},{'Sequencing Technology':'Illumina'},regex=True)
    x = x.replace({'Sequencing Technology':r'^Nextera.*$'},{'Sequencing Technology':'Illumina'},regex=True)
    x = x.replace({'Sequencing Technology':r'^llumina.*$'},{'Sequencing Technology':'Illumina'},regex=True)
    x = x.replace({'Sequencing Technology':r'^[Ii][Oo][Nn].*$'},{'Sequencing Technology':'IonTorrent'},regex=True)
    x = x.replace({'Sequencing Technology':r'^Oxford\sNanopore.*$'},{'Sequencing Technology':'OxfordNanopore'},regex=True)
    x = x.replace({'Sequencing Technology':r'.*Min[Ii][Oo][Nn].*$'},{'Sequencing Technology':'OxfordNanopore'},regex=True)
    x = x.replace({'Sequencing Technology':r'^[Pp]ac[Bb]io.*$'},{'Sequencing Technology':'PacBio'},regex=True)
    x = x.replace({'Sequencing Technology':r'^Sanger.*$'},{'Sequencing Technology':'Sanger'},regex=True)
    x = x.replace({'Sequencing Technology':r'^Solexa.*$'},{'Sequencing Technology':'Solexa'},regex=True)
    x = x.replace({'Sequencing Technology':r'.*[Ss][Oo][Ll]i[Dd].*$'},{'Sequencing Technology':'SOLiD'},regex=True)

    x = x.replace({'Sequencing Technology':r'.*Others.*$'},{'Sequencing Technology':np.nan},regex=True)
    
    return x

#Creating refseq inclusion and exclusion (minus multi-isolates) dataframes
df = pd.read_table(filepath)
refseqY_df = df[df['Refseq Accession'].notnull()]
refseqN_df = df[df['Refseq Accession'].isnull()]
refseqY_df = refseqY_df.assign(datalabel=1)
refseqN_further = refseqN_df[~refseqN_df['Refseq Exclusion Reason'].str.contains('multi-isolate|\[\]')]
refseqN_further = refseqN_further.assign(datalabel=2)

#Data cleaning of assembly methods and sequencing platforms
refseqY_df = assembly_method_categorical(refseqY_df)
refseqN_further = assembly_method_categorical(refseqN_further)
refseqY_df = seq_plat_categorical(refseqY_df)
refseqN_further = seq_plat_categorical(refseqN_further)

#Creating dataframes of inclusion and exclusion sets with columns of interest
inclusion_set = refseqY_df[['logn50_norm','logcontigcount_norm','logl50_norm','logtotlen_norm','Assembly Method','Sequencing Technology','GC_norm','datalabel']]
exclusion_set = refseqN_further[['logn50_norm','logcontigcount_norm','logl50_norm','logtotlen_norm','Assembly Method','Sequencing Technology','GC_norm','datalabel']]

#Identifying index of long read associated assemblers and sequencing technologies
long_read_index_incl = inclusion_set[(inclusion_set['Sequencing Technology'] == 'OxfordNanopore') | \
                                (inclusion_set['Sequencing Technology'] == 'PacBio') | \
                                (inclusion_set['Assembly Method'] == 'Canu') | \
                                (inclusion_set['Assembly Method'] == 'Flye') | \
                                (inclusion_set['Assembly Method'] == 'HGAP') | \
                                (inclusion_set['Assembly Method'] == 'SMRT') | \
                                (inclusion_set['Assembly Method'] == 'pilon')].index

long_read_index_excl = exclusion_set[(exclusion_set['Sequencing Technology'] == 'OxfordNanopore') | \
                                (exclusion_set['Sequencing Technology'] == 'PacBio') | \
                                (exclusion_set['Assembly Method'] == 'Canu') | \
                                (exclusion_set['Assembly Method'] == 'Flye') | \
                                (exclusion_set['Assembly Method'] == 'HGAP') | \
                                (exclusion_set['Assembly Method'] == 'SMRT') | \
                                (exclusion_set['Assembly Method'] == 'pilon')].index

#Excluding rows containing long read data
inclusion_set.drop(long_read_index_incl, inplace=True)
exclusion_set.drop(long_read_index_excl, inplace=True)

#Subsetting dataframe to contain only numerical attributes
inclusion_set_numerical = inclusion_set[['logn50_norm','logcontigcount_norm','logl50_norm','logtotlen_norm','GC_norm','datalabel']]
exclusion_set_numerical = exclusion_set[['logn50_norm','logcontigcount_norm','logl50_norm','logtotlen_norm','GC_norm','datalabel']]
frame = [inclusion_set_numerical, exclusion_set_numerical]
integrated_frame = pd.concat(frame, ignore_index=True)
inclusion = integrated_frame.loc[integrated_frame['datalabel']==1]
exclusion = integrated_frame.loc[integrated_frame['datalabel']==2]

#Subsamping inclusion frame to equalize exclusion frame
inclusion_balanced = inclusion.sample(n=1841)
balanced = [inclusion_balanced, exclusion]
integrated_bal = pd.concat(balanced, ignore_index=True)


#Defining data matrix and label vector
X_bal = integrated_bal.drop('datalabel', axis=1)
Y_bal = integrated_bal['datalabel']

#Some GCnormalized values are missing, so we perform median imputation
from sklearn.impute import SimpleImputer
imp = SimpleImputer(missing_values=np.nan, strategy='median')

imp = imp.fit(X_bal)
X_bal_imp = imp.transform(X_bal)


#training and importing random forest model with coverage included
rfc.fit(X_bal_imp, Y_bal)
joblib.dump(rfc, "random_forest_n50_contigcount_l50_totlen_gccontent.joblib")
