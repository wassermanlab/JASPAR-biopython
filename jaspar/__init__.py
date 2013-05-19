from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna as dna
import re
import math

import sys
#sys.path.append("/homed/home/dave/devel/biopython-1.61")
sys.path.append("/home/anthony/PostDoc/JASPAR2013/Biopython_package/biopython-master/")
from Bio_dev import motifs


class Motif(motifs.Motif):
    """
    A subclass of Motif used in parsing JASPAR files.

    Represents a JASPAR profile as a Motif with additional metadata information
    if available

    """
    def __init__(self, matrix_id, name, alphabet=dna, instances=None,
                 counts=None, collection=None, tf_class=None, tf_family=None,
                 species=None, tax_group=None, acc=None, data_type=None,
                 medline=None, pazar_id=None, comment=None):
        """
        Construct a JASPAR Motif instance.

        """
        motifs.Motif.__init__(self, alphabet, instances, counts)
        self.name = name
        # TODO we may want to store the base ID and version in separate
        # attributes or just provide methods to split/return the matrix ID
        # into base_id and version
        self.matrix_id = matrix_id
        self.collection = collection
        self.tf_class = tf_class
        self.tf_family = tf_family
        self.species = species      # May have multiple so species is a list.
                                    # The species are actually specified as
                                    # taxonomy IDs.
        self.tax_group = tax_group
        self.acc = acc              # May have multiple so acc is a list.
        self.data_type = data_type
        self.medline = medline
        self.pazar_id = pazar_id
        self.comment = comment

    def __str__(self):
        """
        Return a string represention of the JASPAR profile. We choose to
        provide only the filled metadata information.

        """
        tf_name_str = "TF name\t{0}\n".format(self.name)
        matrix_id_str = "Matrix ID\t{0}\n".format(self.matrix_id)
        the_string = "".join([tf_name_str, matrix_id_str])
        if self.collection:
            collection_str = "Collection\t{0}\n".format(self.collection)
            the_string = "".join([the_string, collection_str])
        if self.tf_class:
            tf_class_str = "TF class\t{0}\n".format(self.tf_class)
            the_string = "".join([the_string, tf_class_str])
        if self.tf_family:
            tf_family_str = "TF family\t{0}\n".format(self.tf_family)
            the_string = "".join([the_string, tf_family_str])
        if self.species:
            species_str = "Species\t{0}\n".format(",".join(self.species))
            the_string = "".join([the_string, species_str])
        if self.tax_group:
            tax_group_str = "Taxonomic group\t{0}\n".format(self.tax_group)
            the_string = "".join([the_string, tax_group_str])
        if self.acc:
            acc_str = "Accession\t{0}\n".format(self.acc)
            the_string = "".join([the_string, acc_str])
        if self.data_type:
            data_type_str = "Data type used\t{0}\n".format(self.data_type)
            the_string = "".join([the_string, data_type_str])
        if self.medline:
            medline_str = "Medline\t{0}\n".format(self.medline)
            the_string = "".join([the_string, medline_str])
        if self.pazar_id:
            pazar_id_str = "PAZAR ID\t{0}\n".format(self.pazar_id)
            the_string = "".join([the_string, pazar_id_str])
        if self.comment:
            comment_str = "Comments\t{0}\n".format(self.comment)
            the_string = "".join([the_string, comment_str])
        matrix_str = "Matrix:\n{0}\n\n".format(self.counts)
        the_string = "".join([the_string, matrix_str])
        return the_string

    def __hash__(self):
        """
        Return the hash key corresponding to the JASPAR profile

        :note: We assume the unicity of matrix IDs

        """
        return self.matrix_id.__hash__()

    def __eq__(self, other):
        return self.matrix_id == other.matrix_id

    def ic(self):
        """
        Compute the total information content.
        XXX This really belongs in the base Motif class

        """

        pwm = self.pwm
        alphabet = self.alphabet
        background = self.background
        if background:
            background = dict(background)
        else:
            background = dict.fromkeys(sorted(self.alphabet.letters), 1.0)

        total = sum(background.values())
        for l in alphabet.letters:
            background[l] /= total

        ic = 0
        for i in range(self.length):
            for l in alphabet.letters:
                p = pwm[l][i]
                b = background[l]
                if b > 0:
                    if p > 0:
                        ic += p * math.log(p / b, 2)

        return ic

    def gc_content(self):
        """
        Compute the GC content.
        XXX This really belongs in the base Motif class
        """
        counts = self.counts
        alphabet = self.alphabet
        gc_total = 0
        total = 0
        for i in xrange(self.length):
            for letter in alphabet.letters:
                if letter == 'C' or letter == 'G':
                    gc_total += counts[letter][i]

                total += counts[letter][i]

        return float(gc_total) / total


class Record(list):
    """
    Represents a list of jaspar motifs

    Attribute:
        o version: The JASPAR version used

    """

    def __init__(self):
        self.version = None

    def __str__(self):
        return "\n".join([str(the_motif) for the_motif in self])

    def to_dict(self):
        """
        Return the list of matrices as a dictionnary of matrices

        """

        dic = {}
        for motif in self:
            dic[motif.matrix_id] = motif
        return dic


def read(handle, format):
    """
    Read motif(s) from a file in one of several different JASPAR formats.
    Call the appropriate routine based on the format passed.
    """

    format = format.lower()
    if format == "pfm":
        record = _read_pfm(handle)
        return record
    elif format == "sites":
        record = _read_sites(handle)
        return record
    elif format == "jaspar":
        record = _read_jaspar(handle)
        return record
    else:
        raise ValueError("Unknown format %s" % format)


def write(motif):
    """Returns the pfm representation of the motif
    """
    letters = "ACGT"
    counts = motif.counts
    lines = []
    for letter in letters:
        terms = map(str, counts[letter])
        line = "\t".join(terms) + "\n"
        lines.append(line)
    # Finished; glue the lines together
    text = "".join(lines)
    return text


def _read_pfm(handle):
    """
    Reads the motif from a JASPAR .pfm file
    """

    alphabet = dna
    counts = {}

    letters = "ACGT"
    for letter, line in zip(letters, handle):
        words = line.split()
        #if there is a letter in the beginning, ignore it
        if words[0] == letter:
            words = words[1:]
        counts[letter] = map(float, words)

    motif = Motif(matrix_id=None, name=None, alphabet=alphabet, counts=counts)
    motif.mask = "*" * motif.length
    record = Record()
    record.append(motif)

    return record


def _read_sites(handle):
    """
    Reads the motif from JASPAR .sites file
    """

    alphabet = dna
    instances = []

    for line in handle:
        if not line.startswith(">"):
            break
        # line contains the header ">...."
        # now read the actual sequence
        line = handle.next()
        instance = ""
        for c in line.strip():
            if c == c.upper():
                instance += c
        instance = Seq(instance, alphabet)
        instances.append(instance)

    instances = motifs.Instances(instances, alphabet)
    motif = Motif(
        matrix_id=None, name=None, alphabet=alphabet, instances=instances
    )
    motif.mask = "*" * motif.length
    record = Record()
    record.append(motif)

    return record


def _read_jaspar(handle):
    """
    Read motifs from a JASPAR formatted file

    Format is one or more records of the form, e.g.:
    >MA0001.1 AGL3
    A  [ 0  3 79 40 66 48 65 11 65  0 ]
    C  [94 75  4  3  1  2  5  2  3  3 ]
    G  [ 1  0  3  4  1  0  5  3 28 88 ]
    T  [ 2 19 11 50 29 47 22 81  1  6 ]

    """

    alphabet = dna
    counts = {}

    record = Record()

    head_pat = re.compile(r"^>\s*(\S+)(\s+(\S+))?")
    row_pat = re.compile(r"\s*([ACGT])\s*\[\s*(.*)\s*\]")

    identifier = None
    name = None
    row_count = 0
    for line in handle:
        line.rstrip('\r\n')

        head_match = head_pat.match(line)
        row_match = row_pat.match(line)

        if head_match:
            identifier = head_match.group(1)
            if head_match.group(2):
                name = head_match.group(2)
            else:
                name = identifier
        elif row_match:
            (letter, counts_str) = row_match.group(1, 2)

            words = counts_str.split()

            counts[letter] = map(float, words)

            row_count += 1

            if row_count == 4:
                record.append(Motif(identifier, name, alphabet=alphabet,
                                    counts=counts))

                identifier = None
                name = None
                counts = {}
                row_count = 0

    return record
