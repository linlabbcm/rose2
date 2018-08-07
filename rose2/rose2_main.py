#!/usr/bin/env python

# PROGRAM TO STITCH TOGETHER REGIONS TO FORM ENHANCERS, MAP READ DENSITY TO STITCHED REGIONS, AND
# RANK ENHANCERS BY READ DENSITY TO DISCOVER SUPER-ENHANCERS

import argparse
import sys
import time
import os

from rose2.definitions import ROOT_DIR
from rose2.resources import genemapper, rose2_utils, utils


def parse_args(args=None):
    """Argument parser."""
    parser = argparse.ArgumentParser(
        usage=(
            "rose2 [options]"
            " -g [GENOME]"
            " -i [INPUT_REGION_GFF]"
            " -r [RANKBY_BAM_FILE]"
            " -o [OUTPUT_FOLDER]"
            " [OPTIONAL_FLAGS]"
        )
    )

    # Required flags
    parser.add_argument("-i", "--i", dest="input", default=None, type=str,
                        help=(
                            "Enter a .gff or .bed file of binding sites used to make enhancers"
                        ), required=True)
    parser.add_argument("-r", "--rankby", dest="rankby", default=None, type=str,
                        help="bamfile to rank enhancer by", required=True)
    parser.add_argument("-o", "--out", dest="out", default=None, type=str,
                        help="Enter an output folder", required=True)
    parser.add_argument("-g", "--genome", dest="genome", default=None, type=str,
                        help="Enter the genome build (MM9,MM8,HG18,HG19)", required=True)

    # Optional flags
    parser.add_argument("-b", "--bams", dest="bams", default=None, type=str,
                        help="Enter a comma separated list of additional bam files to map to",
                        required=False)
    parser.add_argument("-c", "--control", dest="control", default='', type=str,
                        help="bamfile to rank enhancer by", required=False)
    parser.add_argument("-s", "--stitch", dest="stitch", default=None, type=str,
                        help=(
                            "Enter a max linking distance for stitching. Default will determine "
                            "optimal stitching parameter"
                        ), required=False)
    parser.add_argument("-t", "--tss", dest="tss", default=0, type=int,
                        help="Enter a distance from TSS to exclude. 0 = no TSS exclusion",
                        required=False)
    parser.add_argument("--mask", dest="mask", default=None, type=str,
                        help=(
                            "Mask a set of regions from analysis. Provide a .bed or .gff of "
                            "masking regions"
                        ), required=False)

    return parser.parse_args(args)


def rose(input_file, rankby, output_folder, genome, bams=None, control='', stitch=None, tss=0,
         mask_file=None):
    """ROSE2 main function."""
    debug = False

    # Making the out folder if it doesn't exist
    out_folder = utils.format_folder(output_folder, True)

    # Figuring out folder schema
    gff_folder = utils.format_folder(out_folder + 'gff/', True)
    mapped_folder = utils.format_folder(out_folder + 'mappedGFF/', True)

    # Getting input file
    if input_file.split('.')[-1] == 'bed':
        # Converting a BED file
        input_gff_name = input_file.split('/')[-1][0:-4]
        input_gff_file = '{}{}.gff'.format(gff_folder, input_gff_name)
        utils.bed_to_gff(input_file, input_gff_file)
    elif input_file.split('.')[-1] == 'gff':
        # Copy the input GFF to the GFF folder
        input_gff_file = input_file
        os.system('cp {} {}'.format(input_gff_file, gff_folder))
    else:
        print('WARNING: INPUT FILE DOES NOT END IN .gff or .bed. ASSUMING .gff FILE FORMAT')
        # Copy the input GFF to the GFF folder
        input_gff_file = input_file
        os.system('cp {} {}'.format(input_gff_file, gff_folder))

    # Getting the list of BAM files to process
    bam_file_list = [rankby]
    if control:
        bam_file_list.append(control)

    if bams:
        bam_file_list += bams.split(',')

    for bam in bam_file_list:
        if not os.path.isfile('{}.bai'.format(bam)):
            print('INDEX FILE FOR {} IS MISSING'.format(bam))
            sys.exit()

    # Optional args

    # Stitch parameter
    stitch_window = '' if not stitch else int(stitch)

    # TSS options
    tss_window = int(tss)
    remove_tss = True if not tss_window else False

    # Getting the Bound region file used to define enhancers
    print('USING {} AS THE INPUT GFF'.format(input_gff_file))
    input_name = input_gff_file.split('/')[-1].split('.')[0]

    # Getting the genome
    print('USING {} AS THE GENOME'.format(genome))

    # Getting the correct annot file
    annotation_path = '{}/annotation'.format(ROOT_DIR)
    genome_dict = {
        'HG18': '{}/hg18_refseq.ucsc'.format(annotation_path),
        'MM9': '{}/mm9_refseq.ucsc'.format(annotation_path),
        'HG19': '{}/hg19_refseq.ucsc'.format(annotation_path),
        'MM8': '{}/mm8_refseq.ucsc'.format(annotation_path),
        'MM10': '{}/mm10_refseq.ucsc'.format(annotation_path),
        'RN4': '{}/rn4_refseq.ucsc'.format(annotation_path),
        'RN6': '{}/rn6_refseq.ucsc'.format(annotation_path),
    }

    annot_file = genome_dict[genome.upper()]

    # Get chroms found in the bams
    print('GETTING CHROMS IN BAMFILES')
    bam_chrom_list = rose2_utils.get_bam_chrom_list(bam_file_list)
    print("USING THE FOLLOWING CHROMS")
    print(bam_chrom_list)

    # Loading in the GFF and filtering by chrom
    print('LOADING AND FILTERING THE GFF')
    input_gff = rose2_utils.filter_gff(input_gff_file, bam_chrom_list)

    # Loading in the bound region reference collection
    print('LOADING IN GFF REGIONS')
    reference_collection = utils.gff_to_locus_collection(input_gff)

    print('STARTING WITH {} INPUT REGIONS'.format(len(reference_collection)))
    print('CHECKING REFERENCE COLLECTION:')
    rose2_utils.check_ref_collection(reference_collection)

    # Masking reference collection
    # See if there's a mask
    if mask_file:
        print('USING MASK FILE {}'.format(mask_file))
        # if it's a bed file
        if mask_file.split('.')[-1].upper() == 'BED':
            mask_gff = utils.bed_to_gff(mask_file)
        elif mask_file.split('.')[-1].upper() == 'GFF':
            mask_gff = utils.parse_table(mask_file, '\t')
        else:
            print("MASK MUST BE A .gff or .bed FILE")

        mask_collection = utils.gff_to_locus_collection(mask_gff)
        print('LOADING {} MASK REGIONS'.format(len(mask_collection)))

        # Now mask the reference loci
        reference_loci = reference_collection.get_loci()
        filtered_loci = [
            locus for locus in reference_loci if not mask_collection.get_overlap(locus, 'both')
        ]
        print(
            "FILTERED OUT {} LOCI THAT WERE MASKED IN {}"
            "".format(len(reference_loci) - len(filtered_loci), mask_file)
        )
        reference_collection = utils.LocusCollection(filtered_loci, 50)

    # Now stitch regions
    print('STITCHING REGIONS TOGETHER')
    stitched_collection, debug_output, stitch_window = rose2_utils.region_stitching(
        reference_collection,
        input_name,
        out_folder,
        stitch_window,
        tss_window,
        annot_file,
        remove_tss,
    )

    # Now make a stitched collection GFF
    print('MAKING GFF FROM STITCHED COLLECTION')
    stitched_gff = utils.locus_collection_to_gff(stitched_collection)
    # Making sure start/stop ordering are correct
    for line in stitched_gff:
        start = int(line[3])
        stop = int(line[4])
        if start > stop:
            line[3] = stop
            line[4] = start

    print(stitch_window)
    print(type(stitch_window))
    if not remove_tss:
        stitched_gff_file = '{}{}_{}KB_STITCHED.gff'.format(
            gff_folder,
            input_name,
            str(stitch_window / 1000),
        )
        stitched_gff_name = '{}_{}KB_STITCHED'.format(input_name, str(stitch_window / 1000))
        debug_out_file = '{}{}_{}KB_STITCHED.debug'.format(
            gff_folder,
            input_name,
            str(stitch_window / 1000),
        )
    else:
        stitched_gff_file = '{}{}_{}KB_STITCHED_TSS_DISTAL.gff'.format(
            gff_folder,
            input_name,
            str(stitch_window / 1000),
        )
        stitched_gff_name = '{}_{}KB_STITCHED_TSS_DISTAL'.format(
            input_name,
            str(stitch_window / 1000),
        )
        debug_out_file = '{}{}_{}KB_STITCHED_TSS_DISTAL.debug'.format(
            gff_folder,
            input_name,
            str(stitch_window / 1000),
        )

    # Writing debug output to disk
    if debug:
        print('WRITING DEBUG OUTPUT TO DISK AS {}'.format(debug_out_file))
        utils.unparse_table(debug_output, debug_out_file, '\t')

    # Write the GFF to disk
    print('WRITING STITCHED GFF TO DISK AS {}'.format(stitched_gff_file))
    utils.unparse_table(stitched_gff, stitched_gff_file, '\t')

    # Setting up the overall output file
    output_file_1 = out_folder + stitched_gff_name + '_ENHANCER_REGION_MAP.txt'
    print('OUTPUT WILL BE WRITTEN TO  {}'.format(output_file_1))

    # Try to use the bamliquidatior_path.py script on cluster, otherwise, failover to local
    # (in path), otherwise fail
    bamliquidator_path = 'bamliquidator_batch'

    bam_file_list_unique = list(bam_file_list)
    bam_file_list_unique = utils.uniquify(bam_file_list_unique)
    # Prevent redundant mapping
    print("MAPPING TO THE FOLLOWING BAMS:")
    print(bam_file_list_unique)
    for bam_file in bam_file_list_unique:

        bam_file_name = bam_file.split('/')[-1]

        # Mapping to the stitched GFF
        mapped_out_1_folder = '{}{}_{}_MAPPED'.format(
            mapped_folder,
            stitched_gff_name,
            bam_file_name,
        )
        mapped_out_1_file = '{}{}_{}_MAPPED/matrix.txt'.format(
            mapped_folder,
            stitched_gff_name,
            bam_file_name,
        )
        if utils.check_output(mapped_out_1_file, 0.2, 0.2):
            print("FOUND {} MAPPING DATA FOR BAM: {}".format(stitched_gff_file, mapped_out_1_file))
        else:
            cmd1 = "{} --sense . -e 200 --match_bamToGFF -r {} -o {} {}".format(
                bamliquidator_path,
                stitched_gff_file,
                mapped_out_1_folder,
                bam_file,
            )
            print(cmd1)

            os.system(cmd1)
            if utils.check_output(mapped_out_1_file, 0.2, 5):
                print(
                    "SUCCESSFULLY MAPPED TO {} FROM BAM: {}"
                    "".format(stitched_gff_file, bam_file_name)
                )
            else:
                print(
                    "ERROR: FAILED TO MAP {} FROM BAM: {}"
                    "".format(stitched_gff_file, bam_file_name)
                )
                sys.exit()

    print('BAM MAPPING COMPLETED NOW MAPPING DATA TO REGIONS')
    # Calculate density by region
    # TODO: Need to fix this function to account for different outputs of liquidator
    rose2_utils.map_collection(
        stitched_collection,
        reference_collection,
        bam_file_list,
        mapped_folder,
        output_file_1,
        ref_name=stitched_gff_name,
    )

    print('CALLING AND PLOTTING SUPER-ENHANCERS')

    control_name = control.split('/')[-1] if control else 'NONE'

    cmd = 'ROSE2_callSuper.R {} {} {} {}'.format(
        out_folder,
        output_file_1,
        input_name,
        control_name,
    )

    print(cmd)

    os.system(cmd)

    # Calling the gene mapper
    time.sleep(5)
    tables = [
        "_SuperEnhancers.table.txt",
        "_StretchEnhancers.table.txt",
        "_SuperStretchEnhancers.table.txt",
    ]
    for table in tables:
        table_file = "{}{}".format(input_name, table)
        genemapper.map(
            os.path.join(out_folder, table_file),
            genome,
            rankby,
            control,
        )


def main(args=None):
    """Main run call."""

    # Retrieving flags
    args = parse_args(args)

    rose(
        args.input,
        args.rankby,
        args.out,
        args.genome,
        bams=args.bams,
        control=args.control,
        stitch=args.stitch,
        tss=args.tss,
        mask_file=args.mask,
    )
