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
sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.output[0], "w") as out:
    print(*snakemake.params.mixtures, sep="\n", file=out)
