import trackhub
from snakemake.script import snakemake
import os

#################
# Initiate Hub
#################
hub = trackhub.Hub(
    hub=snakemake.config["ucsc_trackhub"]["hub_name"],
    short_label=snakemake.config["ucsc_trackhub"]["short_label"],
    long_label=snakemake.config["ucsc_trackhub"]["long_label"],
    email=snakemake.config["ucsc_trackhub"]["email"]
)

#########################################
# Add genomes for genome.txt file
#########################################

for assembly_name, assembly_data in snakemake.config["ucsc_trackhub"]["genomes"].items():

    genome = trackhub.Assembly(
        genome=assembly_name,
        twobit_file=os.path.abspath(snakemake.input.genome_2bit),
        organism=assembly_data["organism"],
        defaultPos=assembly_data["defaultPos"],
        scientificName=assembly_data["scientificName"],
        description=assembly_data["description"],
        html_string=assembly_data["description"],
        orderKey=4800
    )

    genomes_file = trackhub.GenomesFile()
    hub.add_genomes_file(genomes_file)

    # Add TrackDb 
    # we also need to create a trackDb and add it to the genome
    trackdb = trackhub.TrackDb()
    genome.add_trackdb(trackdb)

    # add the genome to the genomes file here:
    genomes_file.add_genome(genome)

    #######################
    # Add data to trackdb.txt
    #######################

    # Add genome model
    genome_model = trackhub.Track(
        name=assembly_data["annotation"]["track_name"],
        tracktype="bigGenePred",
        source=os.path.abspath(snakemake.input.genome_genePred),
        shortLabel=assembly_data["annotation"]["short_label"],
        longLabel=assembly_data["annotation"]["long_label"],
        visibility="pack",
    )

    trackdb.add_tracks(genome_model)

    # Add TSRs track
    tsr_track = trackhub.Track(
        name="TSRs",
        tracktype="bigBed",
        source=os.path.abspath(snakemake.input.tsrs),
        shortLabel="TSRs",
        longLabel="Transcription Start Regions",
        visibility="dense",
    )

    trackdb.add_tracks(tsr_track)

    # Loop through bigwig files in snakemake.input.bw and add to trackhub
    for bw in snakemake.input.plus_bw:
        bw_basename=os.path.basename(bw)
        bw_name=os.path.splitext(bw_basename)[0]

        bw_track = trackhub.Track(
            name=bw_name,
            tracktype="bigWig",
            source=os.path.abspath(bw),
            shortLabel=bw_name,
            longLabel=bw_name,
            visibility="full",
            color="0,0,100"
        )

        trackdb.add_tracks(bw_track)

    for bw in snakemake.input.minus_bw:
        bw_basename=os.path.basename(bw)
        bw_name=os.path.splitext(bw_basename)[0]

        bw_track = trackhub.Track(
            name=bw_name,
            tracktype="bigWig",
            source=os.path.abspath(bw),
            shortLabel=bw_name,
            longLabel=bw_name,
            visibility="full",
            negateValues="on",
            color="113,35,124"
        )

        trackdb.add_tracks(bw_track)

    trackhub.upload.stage_hub(hub, staging=snakemake.output.dir)