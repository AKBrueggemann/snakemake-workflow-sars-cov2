'# Copyright ' + str(datetime.datetime.now().year)  + ' Thomas Battenfeld, Alexander Thomas, Johannes Köster.' 
'# Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)' 
'# This file may not be copied, modified, or distributed' 
'# except according to those terms.
'
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

import pysam

sars_cov2_id, _ = snakemake.params.reference_genome[0].split(".", 1)


def is_sars_cov2(record, mate=False):
    if mate:
        return record.next_reference_name.startswith(sars_cov2_id)
    else:
        return record.reference_name.startswith(sars_cov2_id)


with pysam.AlignmentFile(snakemake.input[0], "rb") as inbam:
    with pysam.AlignmentFile(snakemake.output[0], "wb", template=inbam) as outbam:
        for record in inbam:
            assert record.is_paired
            if (
                (record.is_unmapped and record.mate_is_unmapped)
                or (is_sars_cov2(record) and record.mate_is_unmapped)
                or (is_sars_cov2(record, mate=True) and record.is_unmapped)
                or (is_sars_cov2(record) and is_sars_cov2(record, mate=True))
            ):
                outbam.write(record)
