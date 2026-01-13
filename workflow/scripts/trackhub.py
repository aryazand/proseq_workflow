import trackhub
from snakemake.script import snakemake
import os
from pathlib import Path

hub, genomes_file, genome, trackdb = trackhub.default_hub(
    hub_name=snakemake.params.hub_name,
    defaultPos=snakemake.params.defaultPos,
    short_label=snakemake.params.short_label,
    long_label=snakemake.params.long_label,
    genome=snakemake.params.genome,
    email=snakemake.params.email,
)

# Add TSRs track
tsr_track = trackhub.Track(
    name="TSRs",
    tracktype="bigBed",
    bigDataUrl=os.path.relpath(snakemake.input.tsrs, snakemake.output.trackdb),
    shortLabel="TSRs",
    longLabel="Transcription Start Regions",
    visibility="dense",
)

trackdb.add_tracks(tsr_track)

# Loop through bigwig files in snakemake.input.bw and add to trackhub
for bw in snakemake.input.bw:
    bw_basename=os.path.basename(bw)
    bw_name=os.path.splitext(bw_basename)[0]

    bw_track = trackhub.Track(
        name=bw_name,
        tracktype="bigWig",
        bigDataUrl=os.path.relpath(bw, snakemake.output.trackdb),
        shortLabel=bw_name,
        longLabel=bw_name,
        visibility="full",
    )

    trackdb.add_tracks(bw_track)

trackhub.upload.stage_hub(hub, staging=snakemake.output.dir)