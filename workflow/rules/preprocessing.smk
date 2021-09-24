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
rule update_sample:
    input:
        "config/pep/samples.csv",
    log:
        "logs/sample_update/preprocessing/sample_csv_update.txt",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/update-sample-sheet.py"
