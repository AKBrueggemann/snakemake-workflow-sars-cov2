'# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.' 
'# Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)' 
'# This file may not be copied, modified, or distributed' 
'# except according to those terms.
'
'Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
' 
'Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
' 
'This file may not be copied, modified, or distributed
' 
'except according to those terms.
'
import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd


def get_read_statistics(sm_input, sm_output):
    all_read_counts = []
    for input in sm_input:
        with open(input) as in_file:
            all_read_counts.append(next(in_file).rstrip())
    all_read_counts = pd.Series(all_read_counts, dtype=int)

    with open(sm_output, "w") as f:
        print("Length of all reads:\n", file=f)
        print(all_read_counts.describe().apply(lambda x: format(x, "f")), file=f)


get_read_statistics(snakemake.input, snakemake.output[0])
