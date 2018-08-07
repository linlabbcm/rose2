"""ROSE2 UTILITY FUNCTIONS."""

# Functions require samtools

import copy
import numpy
import subprocess
import sys

from collections import defaultdict
from rose2.resources import utils

# ==================================================================
# =====================HELPER FUNCTIONS=============================
# ==================================================================


def check_ref_collection(reference_collection):
    """Makes sure the names of all loci in the reference collection are unique."""
    names_list = [locus.id for locus in reference_collection.get_loci()]

    if len(names_list) != len(utils.uniquify(names_list)):
        print("ERROR: REGIONS HAVE NON-UNIQUE IDENTIFIERS")
        sys.exit()
    else:
        print("REFERENCE COLLECTION PASSES QC")
        return


def get_bam_chrom_list(bam_file_list):
    """Gets the consensus list of chromosomes mapped by the bams."""
    # Start w/ the first bam
    cmd = 'samtools idxstats {}'.format(bam_file_list[0])
    idx_stats = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    idx_stats = idx_stats.communicate()
    final_chrom_list = [
        line.split('\t')[0] for line in idx_stats[0].decode('utf-8').split('\n')[0:-2]
    ]

    # Now go through each additional bam
    for bam_file in bam_file_list:
        cmd = 'samtools idxstats {}'.format(bam_file)
        idx_stats = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        idx_stats = idx_stats.communicate()
        chrom_list = [
            line.split('\t')[0] for line in idx_stats[0].decode('utf-8').split('\n')[0:-2]
        ]
        final_chrom_list = [chrom for chrom in final_chrom_list if chrom_list.count(chrom)]

    return utils.uniquify(final_chrom_list)


def filter_gff(gff_file, chrom_list):
    """Takes in a gff and filters out all lines that don't belong to a chrom in the chrom_list."""
    gff = utils.parse_table(gff_file, '\t')
    filtered_gff = []
    exclude_list = []
    for line in gff:
        if chrom_list.count(line[0]) == 1:
            filtered_gff.append(line)
        else:
            exclude_list.append(line[0])

    exclude_list = utils.uniquify(exclude_list)
    if exclude_list:
        print("EXCLUDED GFF REGIONS FROM THE FALLING CHROMS: {}".format(','.join(exclude_list)))

    return filtered_gff


# ==================================================================
# =====================REGION STITCHING=============================
# ==================================================================

def optimize_stitching(locus_collection, name, out_folder, step_size=500):
    """Takes a locus collection and starts writing out stitching stats at step sized intervals."""
    max_stitch = 5000  # set a hard wired match stitching parameter

    stitch_table = [
        [
            'STEP',
            'NUM_REGIONS',
            'TOTAL_CONSTIT',
            'TOTAL_REGION',
            'MEAN_CONSTIT',
            'MEDIAN_CONSTIT',
            'MEAN_REGION',
            'MEDIAN_REGION',
            'MEAN_STITCH_FRACTION',
            'MEDIAN_STITCH_FRACTION',
        ]
    ]
    # First consolidate the collection
    locus_collection = locus_collection.stitch_collection(stitch_window=0)
    total_constit = sum([locus.len() for locus in locus_collection.get_loci()])
    step = 0
    while step <= max_stitch:
        print("Getting stitch stats for {} (bp)".format(step))
        stitch_collection = locus_collection.stitch_collection(stitch_window=step)
        num_regions = len(stitch_collection)
        stitch_loci = stitch_collection.get_loci()
        region_lengths = [locus.len() for locus in stitch_loci]
        total_region = sum(region_lengths)
        constit_lengths = []
        for locus in stitch_loci:

            constit_loci = locus_collection.get_overlap(locus)
            constit_lengths.append(sum([locus.len() for locus in constit_loci]))

        mean_constit = round(numpy.mean(constit_lengths), 2)
        median_constit = round(numpy.median(constit_lengths), 2)

        mean_region = round(numpy.mean(region_lengths), 2)
        median_region = round(numpy.median(region_lengths), 2)

        stitch_fractions = [
            float(constit_length) / float(region_length)
            for constit_length, region_length in zip(constit_lengths, region_lengths)
        ]
        mean_stitch_fraction = round(numpy.mean(stitch_fractions), 2)
        median_stitch_fraction = round(numpy.median(stitch_fractions), 2)

        new_line = [
            step,
            num_regions,
            total_constit,
            total_region,
            mean_constit,
            median_constit,
            mean_region,
            median_region,
            mean_stitch_fraction,
            median_stitch_fraction,
        ]

        stitch_table.append(new_line)

        step += step_size

    # Write the stitch table to disk
    stitch_param_file = '{}{}_stitch_params.tmp'.format(out_folder, name)
    utils.unparse_table(stitch_table, stitch_param_file, '\t')
    # Call the rscript
    r_cmd = 'ROSE2_stitchOpt.R {} {} {}'.format(stitch_param_file, out_folder, name)
    print(r_cmd)
    # Get back the stitch parameter
    r_output = subprocess.Popen(r_cmd, stdout=subprocess.PIPE, shell=True)
    r_output_test = r_output.communicate()

    print(r_output_test)

    stitch_param = r_output_test[0].decode('utf-8').split('\n')[2]
    try:
        stitch_param = int(stitch_param)
    except ValueError:
        print("INVALID STITCHING PARAMETER. STITCHING OPTIMIZATION FAILED")
        sys.exit()

    return stitch_param


def region_stitching(reference_collection, name, out_folder, stitch_window, tss_window,
                     annot_file, remove_tss=True):
    """Region stitching."""
    print('PERFORMING REGION STITCHING')
    # First have to turn bound region file into a locus collection

    # Need to make sure this names correctly... Each region should have a unique name
    # reference_collection

    debug_output = []
    # Filter out all bound regions that overlap the TSS of an ACTIVE GENE
    if remove_tss:
        print('REMOVING TSS FROM REGIONS USING AN EXCLUSION WINDOW OF {}BP'.format(tss_window))

        # First make a locus collection of TSS
        start_dict = utils.make_start_dict(annot_file)

        # Now make_TSS loci for active genes
        remove_ticker = 0
        # This loop makes a locus centered around +/- tss_window of transcribed genes then adds it
        # to the list tss_loci
        tss_loci = []
        for gene_id in start_dict.keys():
            tss_loci.append(utils.make_tss_locus(gene_id, start_dict, tss_window, tss_window))

        # This turns the tss_loci list into a LocusCollection
        # 50 is the internal parameter for LocusCollection and doesn't really matter
        tss_collection = utils.LocusCollection(tss_loci, 50)

        # Gives all the loci in reference_collection
        bound_loci = reference_collection.get_loci()

        # This loop will check if each bound region is contained by the TSS exclusion zone
        # This will drop out a lot of the promoter only regions that are tiny
        # Typical exclusion window is around 2kb
        for locus in bound_loci:
            if tss_collection.get_containers(locus, 'both'):
                # If true, the bound locus overlaps an active gene
                reference_collection.remove(locus)
                debug_output.append([str(locus), locus.id, 'CONTAINED'])
                remove_ticker += 1
        print('REMOVED {} LOCI BECAUSE THEY WERE CONTAINED BY A TSS'.format(remove_ticker))

    # Reference_collection is now all enriched region loci that don't overlap an active TSS

    if not stitch_window:
        print('DETERMINING OPTIMUM STITCHING PARAMTER')
        opt_collection = copy.deepcopy(reference_collection)
        stitch_window = optimize_stitching(opt_collection, name, out_folder, step_size=500)
    print('USING A STITCHING PARAMETER OF {}'.format(stitch_window))
    stitched_collection = reference_collection.stitch_collection(stitch_window, 'both')

    if remove_tss:
        # Now replace any stitched region that overlap 2 distinct genes with the original loci
        # that were there
        fixed_loci = []
        tss_loci = []
        for gene_id in start_dict.keys():
            tss_loci.append(utils.make_tss_locus(gene_id, start_dict, 50, 50))

        # This turns the tss_loci list into a LocusCollection
        # 50 is the internal parameter for LocusCollection and doesn't really matter
        tss_collection = utils.LocusCollection(tss_loci, 50)
        remove_ticker = 0
        original_ticker = 0
        for stitched_locus in stitched_collection.get_loci():
            overlapping_tss_loci = tss_collection.get_overlap(stitched_locus, 'both')
            tss_names = [start_dict[tss_locus.id]['name'] for tss_locus in overlapping_tss_loci]
            tss_names = utils.uniquify(tss_names)
            if len(tss_names) > 2:
                # Stitched_collection.remove(stitched_locus)
                original_loci = reference_collection.get_overlap(stitched_locus, 'both')
                original_ticker += len(original_loci)
                fixed_loci += original_loci
                debug_output.append([str(stitched_locus), stitched_locus.id, 'MULTIPLE_TSS'])
                remove_ticker += 1
            else:
                fixed_loci.append(stitched_locus)

        print(
            'REMOVED {} STITCHED LOCI BECAUSE THEY OVERLAPPED MULTIPLE TSSs'.format(remove_ticker)
        )
        print('ADDED BACK {} ORIGINAL LOCI'.format(original_ticker))
        fixed_collection = utils.LocusCollection(fixed_loci, 50)
        return fixed_collection, debug_output, stitch_window

    else:
        return stitched_collection, debug_output, stitch_window


# ==================================================================
# =====================REGION LINKING MAPPING=======================
# ==================================================================

def map_collection(stitched_collection, reference_collection, bam_file_list, mapped_folder,
                   output, ref_name):
    """Makes a table of factor density in a stitched locus.

    Also ranks table by number of loci stitched together.

    """
    print('FORMATTING TABLE')
    loci = stitched_collection.get_loci()

    locus_table = [['REGION_ID', 'CHROM', 'START', 'STOP', 'NUM_LOCI', 'CONSTITUENT_SIZE']]

    loci_len_list = []

    # Strip out any that are in chrY
    for locus in list(loci):
        if locus.chr == 'chrY':
            loci.remove(locus)

    for locus in loci:
        loci_len_list.append(locus.len())

    len_order = utils.order(loci_len_list, decreasing=True)
    ticker = 0
    for i in len_order:
        ticker += 1
        if not ticker % 1000:
            print(ticker)
        locus = loci[i]

        # First get the size of the enriched regions within the stitched locus
        ref_enrich_size = 0
        ref_overlapping_loci = reference_collection.get_overlap(locus, 'both')
        for ref_locus in ref_overlapping_loci:
            ref_enrich_size += ref_locus.len()

        try:
            stitch_count = int(locus.id.split('_')[0])
        except ValueError:
            stitch_count = 1
        coords = [int(x) for x in locus.coords()]

        locus_table.append(
            [locus.id, locus.chr, min(coords), max(coords), stitch_count, ref_enrich_size]
        )

    print('GETTING MAPPED DATA')
    print("USING A BAM FILE LIST:")
    print(bam_file_list)
    for bam_file in bam_file_list:

        bam_file_name = bam_file.split('/')[-1]

        print('GETTING MAPPING DATA FOR  {}'.format(bam_file))
        # Assumes standard convention for naming enriched region gffs

        # Opening up the mapped GFF
        print('OPENING {}{}_{}_MAPPED/matrix.txt'.format(mapped_folder, ref_name, bam_file_name))

        mapped_gff = utils.parse_table(
            '{}{}_{}_MAPPED/matrix.txt'.format(mapped_folder, ref_name, bam_file_name),
            '\t',
        )

        signal_dict = defaultdict(float)
        print('MAKING SIGNAL DICT FOR {}'.format(bam_file))
        mapped_loci = []
        for line in mapped_gff[1:]:

            chrom = line[1].split('(')[0]
            start = int(line[1].split(':')[-1].split('-')[0])
            end = int(line[1].split(':')[-1].split('-')[1])
            mapped_loci.append(utils.Locus(chrom, start, end, '.', line[0]))
            try:
                signal_dict[line[0]] = float(line[2]) * (abs(end - start))
            except ValueError:
                print('WARNING NO SIGNAL FOR LINE:')
                print(line)
                continue

        mapped_collection = utils.LocusCollection(mapped_loci, 500)
        locus_table[0].append(bam_file_name)

        for i in range(1, len(locus_table)):
            signal = 0.0
            line = locus_table[i]
            line_locus = utils.Locus(line[1], line[2], line[3], '.')
            overlapping_regions = mapped_collection.get_overlap(line_locus, sense='both')
            for region in overlapping_regions:
                signal += signal_dict[region.id]
            locus_table[i].append(signal)

    utils.unparse_table(locus_table, output, '\t')
