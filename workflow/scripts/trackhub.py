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
    # Add genome model to trackdb.txt
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

    #######################
    # Add TSRs to trackdb.txt
    #######################

    tsr_track = trackhub.Track(
        name="TSRs",
        tracktype="bigBed",
        source=os.path.abspath(snakemake.input.tsrs),
        shortLabel="TSRs",
        longLabel="Transcription Start Regions",
        visibility="dense",
    )

    trackdb.add_tracks(tsr_track)

    #######################
    # Add BigWig to trackdb.txt
    #######################

    # Loop through bigwig files in snakemake.input.bw and add to trackhub
    for bw in snakemake.input.fiveprime_plus_bw:
        bw_basename=os.path.basename(bw)
        bw_name=os.path.splitext(bw_basename)[0]

        bw_track = trackhub.Track(
            name=bw_name + "_fiveprime",
            tracktype="bigWig",
            source=os.path.abspath(bw),
            shortLabel=bw_name + ("(5')"),
            longLabel=bw_name + ("(5' ends)"),
            visibility="full",
            autoScale="on",
            color=snakemake.config["ucsc_trackhub"]["process_bw"]["plus_color"]
        )

        trackdb.add_tracks(bw_track)

    for bw in snakemake.input.fiveprime_minus_bw:
        bw_basename=os.path.basename(bw)
        bw_name=os.path.splitext(bw_basename)[0]

        bw_track = trackhub.Track(
            name=bw_name + "_fiveprime",
            tracktype="bigWig",
            source=os.path.abspath(bw),
            shortLabel=bw_name + ("(5')"),
            longLabel=bw_name + ("(5' ends)"),
            visibility="full",
            autoScale="on",
            negateValues=snakemake.config["ucsc_trackhub"]["process_bw"]["negateValues_for_minus_strand"],
            color=snakemake.config["ucsc_trackhub"]["process_bw"]["minus_color"]
        )

        trackdb.add_tracks(bw_track)

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
            autoScale="on",
            color=snakemake.config["ucsc_trackhub"]["process_bw"]["plus_color"]
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
            autoScale="on",
            negateValues=snakemake.config["ucsc_trackhub"]["process_bw"]["negateValues_for_minus_strand"],
            color=snakemake.config["ucsc_trackhub"]["process_bw"]["minus_color"]
        )

        trackdb.add_tracks(bw_track)

    #######################
    # Add groups 
    #######################

    group_name=snakemake.config["ucsc_trackhub"]["group_name"]
    group_label=snakemake.config["ucsc_trackhub"]["group_label"]

    grouping = trackhub.groups.GroupDefinition(
        name=group_name,
        label=group_label,
        priority=1,
        default_is_closed=False)

    groups_file = trackhub.groups.GroupsFile([grouping])
    genome.add_groups(groups_file)

    # We can now add the "group" parameter to all the children of the trackDb
    for track in trackdb.children:
        track.add_params(group=group_name)


    #######################
    # Stage Trackhub 
    #######################
    
    trackhub.upload.stage_hub(hub, staging=snakemake.output.dir)