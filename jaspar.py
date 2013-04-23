from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna as dna
from math import sqrt
import re

import sys
sys.path.append("/homed/home/dave/devel/biopython-1.61")
from Bio_dev import motifs


class Motif(motifs.Motif):
    """
    A subclass of Motif used in parsing JASPAR files.

    Represents a JASPAR profile as a Motif with additional metadata information
    if available

    """
    def __init__(self, matrix_id, name, alphabet=dna, instances=None,
                 counts=None, tf_class=None, tf_family=None, species=None,
                 tax_group=None, acc=None, data_type=None, medline=None,
                 pazar_id=None, comment=None):
        """
        Construct a JASPAR Motif instance.

        """
        motifs.Motif.__init__(self, alphabet, counts=counts)
        self.name = name
        self.matrix_id = matrix_id
        # We assume a uniform distribution of the nt in the background
        self.background = dict.fromkeys(dna.letters, 0.25)
        # Number of sequences used to make the counts
        nb_seq = sum([counts[nt][0] for nt in dna.letters])
        self.pseudocounts = dict(
            (nt, self.background[nt] * sqrt(nb_seq)) for nt in dna.letters)
        self.tf_class = tf_class
        self.tf_family = tf_family
        self.species = species  # May have multiple so species is a list
        self.tax_group = tax_group
        self.acc = acc
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


def read(handle, format):
    alphabet = dna
    counts = {}
    if format == "pfm":
        # reads the motif from Jaspar .pfm file
        letters = "ACGT"
        for letter, line in zip(letters, handle):
            words = line.split()
            #if there is a letter in the beginning, ignore it
            if words[0] == letter:
                words = words[1:]
            counts[letter] = map(float, words)
        motif = Motif(
            matrix_id=None, name=None, alphabet=alphabet, counts=counts
        )
        motif.mask = "*" * motif.length
        return motif
    elif format == "sites":
        # reads the motif from Jaspar .sites file
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
        return motif
    elif format == "jaspar":
        motifs = []

        head_pat = re.compile(r"^>\s*(\S+)\s+(\S+)")
        row_pat = re.compile(r"\s*([ACGT])\s*\[\s*(.*)\s*\]")

        id = None
        name = None
        row_count = 0
        for line in handle:
            line.rstrip('\r\n')

            head_match = head_pat.match(line)
            row_match  = row_pat.match(line)

            if head_match:
                (id, name) = head_match.group(1, 2)
            elif row_match:
                (letter, counts_str) = row_match.group(1, 2)

                words = counts_str.split()

                counts[letter] = map(float, words)

                row_count += 1

                if row_count == 4:
                    motifs.append(
                        Motif(id, name, alphabet=alphabet, counts=counts)
                    )
                    
                    id = None
                    name = None
                    counts = {}
                    row_count = 0

        return motifs
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
