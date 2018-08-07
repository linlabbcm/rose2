"""SET OF GENERAL UTILITY FUNCTIONS FOR SEQ DATA."""

# Locus and LocusCollection classes were generously provided by Graham Ruby and Charles Lin
# Additional functions are found from online sources and are noted in the comments

import os
import gzip
import sys
import time

# ==================================================================
# ==========================I/O FUNCTIONS===========================
# ==================================================================


def open_plus(file_name, mode='r'):
    """Open function that can handle gzip files."""
    if file_name.split('.')[-1] == 'gz':
        return gzip.open(file_name, mode + 'b')
    return open(file_name, mode)


def parse_table(file_name, sep, header=False, excel=False):
    """Opens standard delimited files and outputs a nested list."""
    with open_plus(file_name) as infile:
        line = infile.readline()

        # First figure out if this table came from excel and has \r for line breaks
        if line.count('\r') > 0:
            excel = True

    with open_plus(file_name) as infile:
        if header:
            header = infile.readline()  # disposes of the header

        table = []
        if excel:
            lines_raw = infile.readlines()
            for line in lines_raw:
                line = line.replace('\r', '')
                line = line.rstrip()
                table.append(line.split(sep))
        else:
            for line in infile:
                line = line.rstrip().split(sep)
                table.append(line)

        return table


def unparse_table(table, output, sep):
    """Transforms a 'parse_table' output into a text file."""
    with open_plus(output, 'w') as fh_out:
        if not sep:
            for i in table:
                fh_out.write(str(i) + '\n')
        else:
            for line in table:
                fh_out.write(sep.join([str(x) for x in line]) + '\n')


def format_bed(bed, output=''):
    """Formats a bed file from UCSC or MACS into a WUSTL gateway compatible bed."""
    new_bed = []

    if isinstance(bed, str):
        bed = parse_table(bed, '\t')

    index_ticker = 1
    for line in bed:
        new_line = line[0:4]
        try:
            strand = line[5]
        except IndexError:
            strand = '.'
        new_line += [index_ticker, strand]
        index_ticker += 1
        new_bed.append(new_line)

    if output:
        unparse_table(new_bed, output, '\t')
    else:
        return new_bed


def bed_to_gff(bed, output=''):
    """Turns a bed into a gff file."""
    if isinstance(bed, str):
        bed = parse_table(bed, '\t')

    bed = format_bed(bed)

    # Determine if this is a long bed or a short bed
    if len(bed[0]) == 6:  # this is a full format bed
        bed_style = 'long'
    elif len(bed[0]) == 5:  # this is the medium  length bed with strand
        bed_style = 'medium'
    elif len(bed[0]) == 3:  # this is the minimum length bed
        bed_style = 'short'
    else:
        print('this is probably not actually a bed')
        print(bed[0])
        print('exiting now because the bed is sad')
        sys.exit()
    print('this bed has {} columns and is a {} bed'.format(len(bed[0]), bed_style))

    gff = []
    for line in bed:
        if bed_style == 'long':
            gff_line = [
                line[0],
                line[3],
                '',
                line[1],
                line[2],
                line[4],
                line[5],
                '',
                line[3],
            ]
        elif bed_style == 'medium':
            gff_line = [line[0], '', '', line[1], line[2], '', line[4], '', '']
        else:
            gff_line = [line[0], '', '', line[1], line[2], '', '.', '', '']
        gff.append(gff_line)

    if output:
        unparse_table(gff, output, '\t')
    else:
        return gff


def format_folder(folder_name, create=False):
    """Checks if a folder exists and returns a bool.

    If 'create=True' it creates a folder if it doesn't exist already.

    """
    if folder_name[-1] != '/':
        folder_name += '/'

    try:
        os.listdir(folder_name)
        return folder_name
    except OSError:
        print('Folder {} does not exist'.format(folder_name))
        if create:
            os.mkdir(folder_name)
            return folder_name
    return False


def check_output(file_name, wait_time=1, time_out=30):
    """Checks for the presence of a file every N minutes.

    If it exists, returns True.
    Default is 1 minute with a max timeOut of 30 minutes.

    """
    wait_time = int(wait_time * 60)

    time_out = int(time_out * 60)

    max_ticker = time_out / wait_time
    ticker = 0

    file_exists = False
    while not file_exists:
        try:
            size1 = os.stat(file_name).st_size
            time.sleep(.5)
            size2 = os.stat(file_name).st_size
            if size1 == size2:
                file_exists = True
            else:
                time.sleep(wait_time)
                ticker += 1
        except OSError:
            time.sleep(wait_time)
            ticker += 1
        if ticker == max_ticker:
            break

    time.sleep(.1)
    if file_exists:
        return True
    else:
        print('OPERATION TIMED OUT. FILE {} NOT FOUND'.format(file_name))
        return False


def get_parent_folder(input_file):
    """Returns the parent folder for any file."""
    parent_folder = '/'.join(input_file.split('/')[:-1]) + '/'
    if parent_folder == '/':
        return './'
    else:
        return parent_folder


# ==================================================================
# ===================ANNOTATION FUNCTIONS===========================
# ==================================================================


def make_start_dict(annot_file, gene_list=[]):
    """Makes a dictionary keyed by refseq ID.

    It contains information about chrom/start/stop/strand/common name.

    """
    if isinstance(gene_list, str):
        gene_list = parse_table(gene_list, '\t')
        gene_list = [line[0] for line in gene_list]

    if annot_file.upper().count('REFSEQ') == 1:
        refseq_table, refseq_dict = import_refseq(annot_file)
        if not gene_list:
            gene_list = refseq_dict.keys()

        start_dict = {}
        for gene in gene_list:
            if gene not in refseq_dict:
                continue
            start_dict[gene] = {}
            start_dict[gene]['sense'] = refseq_table[refseq_dict[gene][0]][3]
            start_dict[gene]['chr'] = refseq_table[refseq_dict[gene][0]][2]
            start_dict[gene]['start'] = get_tsss([gene], refseq_table, refseq_dict)
            if start_dict[gene]['sense'] == '+':
                start_dict[gene]['end'] = [int(refseq_table[refseq_dict[gene][0]][5])]
            else:
                start_dict[gene]['end'] = [int(refseq_table[refseq_dict[gene][0]][4])]
            start_dict[gene]['name'] = refseq_table[refseq_dict[gene][0]][12]

    return start_dict


def get_tsss(gene_list, refseq_table, refseq_dict):
    """Generic function to get the TSS of any gene."""
    if not gene_list:
        refseq = refseq_table
    else:
        refseq = refseq_from_key(gene_list, refseq_dict, refseq_table)
    tss = []
    for line in refseq:
        if line[3] == '+':
            tss.append(line[4])
        if line[3] == '-':
            tss.append(line[5])
    tss = map(int, tss)

    return list(tss)


def import_refseq(refseq_file, return_multiples=False):
    """Opens up a refseq file downloaded by UCSC.

    Takes in a refseq table and makes a refseq table and a refseq dictionary for keying the table.

    """
    refseq_table = parse_table(refseq_file, '\t')
    refseq_dict = {}
    ticker = 1
    for line in refseq_table[1:]:
        if line[1] in refseq_dict:
            refseq_dict[line[1]].append(ticker)
        else:
            refseq_dict[line[1]] = [ticker]
        ticker = ticker + 1

    multiples = []
    for i in refseq_dict:
        if len(refseq_dict[i]) > 1:
            multiples.append(i)

    if return_multiples:
        return refseq_table, refseq_dict, multiples
    else:
        return refseq_table, refseq_dict


def refseq_from_key(refseq_key_list, refseq_dict, refseq_table):
    """Function that grabs refseq lines from refseq IDs."""
    type_refseq = []
    for name in refseq_key_list:
        if name in refseq_dict:
            type_refseq.append(refseq_table[refseq_dict[name][0]])
    return type_refseq


def make_transcript_collection(annot_file, up_search, down_search, window=500, gene_list=[]):
    """Makes a LocusCollection w/ each transcript as a locus.

    Takes in either a refseqfile or an ensemblGFF.

    """
    if annot_file.upper().count('REFSEQ') == 1:
        refseq_table, refseq_dict = import_refseq(annot_file)
        locus_list = []
        ticker = 0
        if not gene_list:
            gene_list = refseq_dict.keys()

        for line in refseq_table[1:]:
            if line[1] in gene_list:
                if line[3] == '-':
                    locus = Locus(
                        line[2],
                        int(line[4]) - down_search,
                        int(line[5]) + up_search,
                        line[3],
                        line[1],
                    )
                else:
                    locus = Locus(
                        line[2],
                        int(line[4]) - up_search,
                        int(line[5]) + down_search,
                        line[3],
                        line[1],
                    )
                locus_list.append(locus)
                ticker = ticker + 1
                if not ticker % 1000:
                    print(ticker)

    trans_collection = LocusCollection(locus_list, window)
    return trans_collection


# ==================================================================
# ========================LOCUS INSTANCE============================
# ==================================================================

# Locus and LocusCollection instances courtesy of Graham Ruby

class Locus(object):
    """Standard locus class for tracking genomic loci."""
    __chr_dict = dict()
    __sense_dict = {'+': '+', '-': '-', '.': '.'}

    def __init__(self, chr, start, end, sense, ID='', score=0):
        """Initialize attributes."""
        coords = sorted([int(start), int(end)])
        # This method for assigning chromosome should help avoid storage of redundant strings.
        if chr not in self.__chr_dict:
            self.__chr_dict[chr] = chr
        self.chr = self.__chr_dict[chr]
        self.sense = self.__sense_dict[sense]
        self.start = int(coords[0])
        self.end = int(coords[1])
        self.id = ID
        self.score = score

    def len(self):
        """get locus length."""
        return self.end - self.start + 1

    def get_antisense_locus(self):
        """Get antisense locus."""
        if self.sense == '.':
            return self
        switch = {'+': '-', '-': '+'}
        return Locus(self.chr, self.start, self.end, switch[self.sense])

    def coords(self):
        """Returns a sorted list of the coordinates."""
        return [self.start, self.end]

    def overlaps(self, other_locus):
        """Find if loci overlap.

        True if two loci share any coordinates in common.

        """
        if self.chr != other_locus.chr:
            return False
        if not(self.sense == '.' or
               other_locus.sense == '.' or
               self.sense == other_locus.sense):
            return False
        if self.start > other_locus.end or other_locus.start > self.end:
            return False
        return True

    def contains(self, other_locus):
        """Find if locus contains other locus.

        True if all the nucleotides of the given locus overlap with the self locus.

        """
        if self.chr != other_locus.chr:
            return False
        if not(self.sense == '.' or
               other_locus.sense == '.' or
               self.sense == other_locus.sense):
            return False
        if self.start > other_locus.start or other_locus.end > self.end:
            return False
        return True

    def overlaps_antisense(self, other_locus):
        """Same as overlaps, but considers the opposite strand."""
        return self.get_antisense_locus().overlaps(other_locus)

    def contains_antisense(self, other_locus):
        """Same as contains, but considers the opposite strand."""
        return self.get_antisense_locus().contains(other_locus)

    def __hash__(self):
        """Define hash."""
        return self.start + self.end

    def __eq__(self, other):
        """Define Locus objects equality."""
        if self.__class__ != other.__class__:
            return False
        if self.chr != other.chr:
            return False
        if self.start != other.start:
            return False
        if self.end != other.end:
            return False
        if self.sense != other.sense:
            return False
        return True

    def __ne__(self, other):
        """Opposite to equality."""
        return not self.__eq__(other)

    def __str__(self):
        """Locus string representation."""
        return self.chr + '(' + self.sense + '):' + '-'.join(map(str, self.coords()))


class LocusCollection(object):
    """A collection of locus objects used for querying large sets of loci."""
    def __init__(self, loci, window_size=50):
        """Initialize attributes."""
        # Top-level keys are chr, then strand, no space
        self.__chr_to_coord_to_loci = dict()
        self.__loci = dict()
        self.__win_size = window_size
        for lcs in loci:
            self.__add_locus(lcs)

    def __add_locus(self, lcs):
        """Add locus to the collection."""
        if lcs not in self.__loci:
            self.__loci[lcs] = None
            if lcs.sense == '.':
                chr_key_list = [lcs.chr + '+', lcs.chr + '-']
            else:
                chr_key_list = [lcs.chr + lcs.sense]
            for chr_key in chr_key_list:
                if chr_key not in self.__chr_to_coord_to_loci:
                    self.__chr_to_coord_to_loci[chr_key] = dict()
                for num in self.__get_key_range(lcs):
                    if num not in self.__chr_to_coord_to_loci[chr_key]:
                        self.__chr_to_coord_to_loci[chr_key][num] = []
                    self.__chr_to_coord_to_loci[chr_key][num].append(lcs)

    def __get_key_range(self, locus):
        """Get key range."""
        start = locus.start / self.__win_size
        end = locus.end / self.__win_size + 1  # add 1 because of the range
        return range(int(start), int(end))

    def __len__(self):
        """Define length."""
        return len(self.__loci)

    def append(self, new):
        """Append locus."""
        self.__add_locus(new)

    def has_locus(self, locus):
        """Check if locus in the loci dict"""
        return locus in self.__loci

    def remove(self, old):
        """Remove locus from the collection."""
        if old not in self.__loci:
            raise ValueError("Requested locus isn't in collection.")
        del self.__loci[old]
        if old.sense == '.':
            sense_list = ['+', '-']
        else:
            sense_list = [old.sense]
        for k in self.__get_key_range(old):
            for sense in sense_list:
                self.__chr_to_coord_to_loci[old.chr + sense][k].remove(old)

    def get_loci(self):
        """Get loci."""
        return self.__loci.keys()

    def __subset_helper(self, locus, sense):
        """"""
        sense = sense.lower()
        if ['sense', 'antisense', 'both'].count(sense) != 1:
            raise ValueError("Sense command invalid: '{}'.".format(sense))
        matches = dict()
        senses = ['+', '-']
        if locus.sense == '.' or sense == 'both':
            def lamb(s): return True
        elif sense == 'sense':
            def lamb(s): return bool(s == locus.sense)
        elif sense == 'antisense':
            def lamb(s): return bool(s != locus.sense)
        else:
            raise ValueError("Sense value was inappropriate: '{}'.".format(sense))

        for s in filter(lamb, senses):
            chr_key = locus.chr + s
            if chr_key in self.__chr_to_coord_to_loci:
                for n in self.__get_key_range(locus):
                    if n in self.__chr_to_coord_to_loci[chr_key]:
                        for lcs in self.__chr_to_coord_to_loci[chr_key][n]:
                            matches[lcs] = None

        return matches.keys()

    def get_overlap(self, locus, sense='sense'):
        """Returns all members of the collection that overlap the locus.

        :param sense: can be 'sense' (default), 'antisense', or 'both'.

        """
        matches = self.__subset_helper(locus, sense)
        # Now, get rid of the ones that don't really overlap
        real_matches = dict()
        if sense == 'sense' or sense == 'both':
            for i in filter(lambda lcs: lcs.overlaps(locus), matches):
                real_matches[i] = None
        if sense == 'antisense' or sense == 'both':
            for i in filter(lambda lcs: lcs.overlaps_antisense(locus), matches):
                real_matches[i] = None
        return real_matches.keys()

    def get_containers(self, locus, sense='sense'):
        """Returns all members of the collection that contain the locus.

        Sense can be 'sense' (default), 'antisense', or 'both'.

        """
        matches = self.__subset_helper(locus, sense)
        # Now, get rid of the ones that don't really overlap
        real_matches = dict()
        if sense == 'sense' or sense == 'both':
            for i in filter(lambda lcs: lcs.contains(locus), matches):
                real_matches[i] = None
        if sense == 'antisense' or sense == 'both':
            for i in filter(lambda lcs: lcs.contains_antisense(locus), matches):
                real_matches[i] = None
        return real_matches.keys()

    def stitch_collection(self, stitch_window=1, sense='both'):
        """Reduces the collection by stitching together overlapping loci.

        Returns a new collection. Initializing stichWindow to 1 helps collect directly adjacent
        loci.

        :param sense: can be 'sense' (default), 'antisense', or 'both'.

        """
        locus_list = self.get_loci()
        old_collection = LocusCollection(locus_list, 500)

        stitched_collection = LocusCollection([], 500)

        for locus in locus_list:
            if old_collection.has_locus(locus):
                old_collection.remove(locus)
                overlapping_loci = old_collection.get_overlap(
                    Locus(
                        locus.chr,
                        locus.start - stitch_window,
                        locus.end + stitch_window,
                        locus.sense,
                        locus.id,
                    ),
                    sense,
                )

                stitch_ticker = 1
                while overlapping_loci:
                    stitch_ticker += len(overlapping_loci)
                    overlap_coords = locus.coords()

                    for overlapping_locus in overlapping_loci:
                        overlap_coords += overlapping_locus.coords()
                        old_collection.remove(overlapping_locus)
                    overlap_coords = [int(x) for x in overlap_coords]
                    if sense == 'both':
                        locus = Locus(
                            locus.chr,
                            min(overlap_coords),
                            max(overlap_coords),
                            '.',
                            locus.id,
                        )
                    else:
                        locus = Locus(
                            locus.chr,
                            min(overlap_coords),
                            max(overlap_coords),
                            locus.sense,
                            locus.id,
                        )
                    overlapping_loci = old_collection.get_overlap(
                        Locus(
                            locus.chr,
                            locus.start - stitch_window,
                            locus.end + stitch_window,
                            locus.sense,
                        ),
                        sense
                    )
                locus.id = '{}_{}_lociStitched'.format(stitch_ticker, locus.id)

                stitched_collection.append(locus)

            else:
                continue

        return stitched_collection


# ==================================================================
# ========================LOCUS FUNCTIONS===========================
# ==================================================================

def locus_collection_to_gff(locus_collection):
    """Get a gff format list from a LocusCollection object."""
    loci_list = locus_collection.get_loci()
    gff = []
    for locus in loci_list:
        new_line = [
            locus.chr,
            locus.id,
            '',
            locus.coords()[0],
            locus.coords()[1],
            '',
            locus.sense,
            '',
            locus.id,
        ]
        gff.append(new_line)
    return gff


def gff_to_locus_collection(gff, window=500):
    """Opens up a gff file and turns it into a LocusCollection instance."""
    loci_list = []
    if isinstance(gff, str):
        gff = parse_table(gff, '\t')

    for line in gff:
        # Use line[1] as the locus ID.  If that is empty use line[8]
        if line[1]:
            name = line[1]
        elif line[8]:
            name = line[8]
        else:
            name = '{}:{}:{}-{}'.format(line[0], line[6], line[3], line[4])

        loci_list.append(Locus(line[0], line[3], line[4], line[6], name))
    return LocusCollection(loci_list, window)


def make_tss_locus(gene, start_dict, upstream, downstream):
    """Given a start_dict, make a locus for any gene's TSS w/ upstream and downstream windows."""
    start = start_dict[gene]['start'][0]
    if start_dict[gene]['sense'] == '-':
        return Locus(start_dict[gene]['chr'], start - downstream, start + upstream, '-', gene)
    else:
        return Locus(start_dict[gene]['chr'], start - upstream, start + downstream, '+', gene)


def make_search_locus(locus, up_search, down_search):
    """Takes a locus and expands it by a fixed upstream/downstream amount.

    Spits out the new larger locus.

    """
    if locus.sense == '-':
        search_locus = Locus(
            locus.chr,
            locus.start - down_search,
            locus.end + up_search,
            locus.sense,
            locus.id,
        )
    else:
        search_locus = Locus(
            locus.chr,
            locus.start - up_search,
            locus.end + down_search,
            locus.sense,
            locus.id,
        )
    return search_locus


# ==================================================================
# ========================MISC FUNCTIONS============================
# ==================================================================

# uniquify function
# by Peter Bengtsson
# Used under a creative commons license
# sourced from  here: http://www.peterbe.com/plog/uniqifiers-benchmark

def uniquify(seq, idfun=None):
    """Makes a list unique.

    Order preserving.

    """
    seen = {}
    result = []
    for item in seq:
        if idfun is None:
            marker = item
        else:
            marker = idfun(item)
        if marker in seen:
            continue
        seen[marker] = 1
        result.append(item)
    return result


# Taken from http://code.activestate.com/recipes/491268/

def order(x, none_is_last=True, decreasing=False):
    """Returns the ordering of the elements of x.

    The list [ x[j] for j in order(x) ] is a sorted version of x.

    Missing values in x are indicated by None. If none_is_last is true, then missing values are
    ordered to be at the end. Otherwise, they are ordered at the beginning.

    """
    omit_none = False
    if none_is_last is None:
        none_is_last = True
        omit_none = True

    n = len(x)
    ix = range(n)
    if None not in x:
        ix.sort(reverse=decreasing, key=lambda j: x[j])
    else:
        # Handle None values properly.
        def key(i, x=x):
            elem = x[i]
            # Valid values are True or False only.
            if decreasing == none_is_last:
                return not(elem is None), elem
            else:
                return elem is None, elem
        ix = range(n)
        ix.sort(key=key, reverse=decreasing)

    if omit_none:
        n = len(x)
        for i in range(n-1, -1, -1):
            if x[ix[i]] is None:
                n -= 1
        return ix[:n]
    return ix
