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

import sys
from gc_content import CalculateGCContent

"""
This program takes a spreadsheet of NCBI assemblies metadata, computes genomic attributes of interest for every assembly
and writes the genomic attributes as additional columns within the spreadsheet.
Currently, the only calculated attribute (in addition to NCBI derived attributes) is overall GC content,
which is the fraction of G or C bases of all nucleotide bases of an assembly.
"""

email = sys.argv[1]
api_key = sys.argv[2]
filename = sys.argv[3]
output_dir = sys.argv[4]

calculate_gc_content = CalculateGCContent(email, api_key, filename, output_dir)
calculate_gc_content.append_gc_content()
