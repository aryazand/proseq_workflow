import trackhub
from snakemake.script import snakemake
import os

#################
# Initiate Hub
#################
hub = trackhub.Hub(
    hub=snakemake.config["ucsc_trackhub"]["hub_file"]["hub_name"],
    short_label=snakemake.config["ucsc_trackhub"]["hub_file"]["short_label"],
    long_label=snakemake.config["ucsc_trackhub"]["hub_file"]["long_label"],
    email=snakemake.config["ucsc_trackhub"]["hub_file"]["email"]
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
        html_string=assembly_data["htmlDocumentation"],
        orderKey=assembly_data["orderKey"]
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

    # Creat Annotation group
    annotation_group_name = snakemake.config["ucsc_trackhub"]["hub_file"]["hub_name"] + "_annotations"
    annotation_group_label = snakemake.config["ucsc_trackhub"]["hub_file"]["short_label"] + " Annotations"

    annotation_group = trackhub.groups.GroupDefinition(
        name=annotation_group_name,
        label=annotation_group_label,
        priority=1,
        default_is_closed=False)
    
    # Add genome model
    genome_model = trackhub.Track(
        name=assembly_data["trackDb"]["annotation"]["track_name"],
        tracktype="bigGenePred",
        source=os.path.abspath(snakemake.input.genome_genePred),
        shortLabel=assembly_data["trackDb"]["annotation"]["short_label"],
        longLabel=assembly_data["trackDb"]["annotation"]["long_label"],
        visibility="pack",
    )

    genome_model.add_params(group=annotation_group_name)
    trackdb.add_tracks(genome_model)

    # Add Stringtie Annotation tracks

    for st in snakemake.input.stringtie:

        st_basename=os.path.basename(st)
        st_name=os.path.splitext(st_basename)[0]

        stringtie_track = trackhub.Track(
            name=st_name + "_transcripts",
            tracktype="bigBed",
            source=os.path.abspath(st),
            shortLabel=st_name + " transcripts",
            longLabel=st_name + " transcripts predicted by stringtie",
            visibility="pack",
        )

        stringtie_track.add_params(group=annotation_group_name)
        trackdb.add_tracks(stringtie_track)
    
    #######################
    # Add TSRs to trackdb.txt
    #######################

    # Create TSR group
    tsrs_group_name = snakemake.config["ucsc_trackhub"]["hub_file"]["hub_name"] + "_tsrs"
    tsrs_group_label = snakemake.config["ucsc_trackhub"]["hub_file"]["short_label"] + " TSRs"

    tsrs_group = trackhub.groups.GroupDefinition(
        name=tsrs_group_name,
        label=tsrs_group_label,
        priority=2,
        default_is_closed=False)
    
    for tsr in snakemake.input.tsrs:

        tsr_basename=os.path.basename(tsr)
        tsr_name=os.path.splitext(tsr_basename)[0]

        tsr_track = trackhub.Track(
            name=tsr_name + "_tsrs",
            tracktype="bigBed",
            source=os.path.abspath(tsr),
            shortLabel="TSRs in " + tsr_name,
            longLabel="TSRs detected in " + tsr_name,
            visibility="pack",
            spectrum="on"
        )

        tsr_track.add_params(group=tsrs_group_name)
        trackdb.add_tracks(tsr_track)

    #######################
    # Add BigWig to trackdb.txt
    #######################

    bw_group_name = snakemake.config["ucsc_trackhub"]["hub_file"]["hub_name"] + "_bw"
    bw_group_label = snakemake.config["ucsc_trackhub"]["hub_file"]["short_label"] + " BigWigs"

    bw_group = trackhub.groups.GroupDefinition(
        name=bw_group_name,
        label=bw_group_label,
        priority=3,
        default_is_closed=False)
    
    fiveprime_bw_group_name = snakemake.config["ucsc_trackhub"]["hub_file"]["hub_name"] + "_fiveprime_bw"
    fiveprime_bw_group_label = snakemake.config["ucsc_trackhub"]["hub_file"]["short_label"] + "5' ends BigWigs"

    fiveprime_bw_group = trackhub.groups.GroupDefinition(
        name=fiveprime_bw_group_name,
        label=fiveprime_bw_group_label,
        priority=4,
        default_is_closed=False)
    
    threeprime_bw_group_name = snakemake.config["ucsc_trackhub"]["hub_file"]["hub_name"] + "_threprime__bw"
    threeprime_bw_group_label = snakemake.config["ucsc_trackhub"]["hub_file"]["short_label"] + "3' ends BigWigs"

    threeprime_bw_group = trackhub.groups.GroupDefinition(
        name=threeprime_bw_group_name,
        label=threeprime_bw_group_label,
        priority=5,
        default_is_closed=False)
    
    
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
            maxHeightPixels="100:50:8",
            color=assembly_data["trackDb"]["bw"]["plus_color"]
        )
        
        trackdb.add_tracks(bw_track)
        bw_track.add_params(group=fiveprime_bw_group_name)


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
            maxHeightPixels="100:50:8",
            negateValues=assembly_data["trackDb"]["bw"]["negateValues_for_minus_strand"],
            color=assembly_data["trackDb"]["bw"]["minus_color"]
        )

        trackdb.add_tracks(bw_track)
        bw_track.add_params(group=fiveprime_bw_group_name)
        
    for bw in snakemake.input.threeprime_plus_bw:
        bw_basename=os.path.basename(bw)
        bw_name=os.path.splitext(bw_basename)[0]

        bw_track = trackhub.Track(
            name=bw_name + "_threeprime",
            tracktype="bigWig",
            source=os.path.abspath(bw),
            shortLabel=bw_name + ("(3')"),
            longLabel=bw_name + ("(3' ends)"),
            visibility="full",
            autoScale="on",
            maxHeightPixels="100:50:8",
            color=assembly_data["trackDb"]["bw"]["plus_color"]
        )
        
        trackdb.add_tracks(bw_track)
        bw_track.add_params(group=threeprime_bw_group_name)


    for bw in snakemake.input.threeprime_minus_bw:
        bw_basename=os.path.basename(bw)
        bw_name=os.path.splitext(bw_basename)[0]

        bw_track = trackhub.Track(
            name=bw_name + "_threeprime",
            tracktype="bigWig",
            source=os.path.abspath(bw),
            shortLabel=bw_name + ("(3')"),
            longLabel=bw_name + ("(3' ends)"),
            visibility="full",
            autoScale="on",
            maxHeightPixels="100:50:8",
            negateValues=assembly_data["trackDb"]["bw"]["negateValues_for_minus_strand"],
            color=assembly_data["trackDb"]["bw"]["minus_color"]
        )

        trackdb.add_tracks(bw_track)
        bw_track.add_params(group=threeprime_bw_group_name)

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
            maxHeightPixels="100:50:8",
            color=assembly_data["trackDb"]["bw"]["plus_color"]
        )
        
        trackdb.add_tracks(bw_track)
        bw_track.add_params(group=bw_group_name)

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
            maxHeightPixels="100:50:8",
            negateValues=assembly_data["trackDb"]["bw"]["negateValues_for_minus_strand"],
            color=assembly_data["trackDb"]["bw"]["minus_color"]
        )

        trackdb.add_tracks(bw_track)
        bw_track.add_params(group=bw_group_name)

    #######################
    # Define Group File
    ####################### 

    # Add group file
    groups_file = trackhub.groups.GroupsFile([annotation_group, tsrs_group, bw_group])
    genome.add_groups(groups_file)

    #######################
    # Stage Trackhub 
    #######################

    trackhub.upload.stage_hub(hub, staging=snakemake.output.dir)