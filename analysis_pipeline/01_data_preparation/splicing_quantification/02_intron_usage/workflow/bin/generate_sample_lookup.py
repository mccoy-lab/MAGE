#!/usr/bin/env python3

import sys
import os
import glob
import pandas as pd

leafcutter_dir = sys.argv[1]
metadata_file = sys.argv[2]

# Read metadata into a pandas dataframe
metadata = pd.read_csv(metadata_file, sep='\t')
metadata_dict = dict(zip(metadata["internal_libraryID"], metadata["sample_kgpID"]))

# Process files
for filepath in glob.glob(os.path.join(leafcutter_dir, '*.sorted.gz')):
    internal_id = os.path.basename(filepath).split('.')[0]
    tgp_id = metadata_dict.get(internal_id, '')
    print(f"{internal_id}\t{internal_id}")
