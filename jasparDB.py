import MySQLdb as mdb

from Bio.Alphabet.IUPAC import unambiguous_dna as dna

import sys
sys.path.append("/homed/home/dave/devel/biopython-1.61")
#sys.path.append("/home/anthony/PostDoc/JASPAR2013/Biopython_package/biopython-master/")
from Bio_dev.motifs import jaspar, matrix
from warnings import warn

DFLT_JASPAR_COLLECTION = 'CORE'

class JASPAR5(object):
    """
    Class representing a JASPAR5 DB. VERY loosely based on the perl
    TFBS::DB::JASPAR5 module.

    Note: We will only implement reading of JASPAR motifs from the DB.
    Unlike the perl module, we will not attempt to implement any methods to
    store JASPAR motifs or create a new DB at this time.
    """

    def __init__(self, host=None, name=None, user=None, password=None):
        """
        Construct a JASPAR5 instance and connect to specified DB
        """

        self.name = name
        self.host = host
        self.user = user
        self.password = password

        self.dbh = mdb.connect(host, user, password, name)

    def __str__(self):
        """
        Return a string represention of the JASPAR5 DB connection.
        
        ### TODO
            This should be in some sort of standardized format, e.g. a
            dbi connect string
        """

        str = "%s\@%s:%s" % (self.user, self.host, self.name)

        return str

    def __hash__(self):
        """
        ### TODO Determine if we even need this.
        """

        ### TODO Is this correct?
        return hash(self.__str__())

    def fetch_motif_by_id(self, id):
        """
        Fetch a single JASPAR motif from the DB by it's JASPAR matrix ID
        (e.g. 'MA0001.1').

        NOTE: The perl TFBS module allows you to specify the type of matrix to
        return (PFM, PWM, ICM) but matrices are always stored in JASAPR as
        PFMs so this does not really belong here. Once a PFM is fetched the
        pwm() and pssm() methods can be called to return the normalized and
        log-odds matrices.
        """
         
        if not id:
            # TODO Throw some sort of exception if no id is provided?
            pass

        # separate stable ID and version number
        (base_id, version) = _split_jaspar_id(id) 
        if not version:
            # if ID contains no version portion, fetch latest version by default
            version = self._fetch_latest_version(base_id)

        # fetch internal JASPAR matrix ID - also a check for validity
        int_id = self._fetch_internal_id(base_id, version)

        # fetch JASPAR motif using internal ID
        motif = self._fetch_motif_by_internal_id(int_id)

        return motif

    def fetch_motif_by_name(self, name):
        """
        Fetch a single JASPAR motif from the DB by the TF name
        (e.g. RUNX1).
        """
         
        if not name:
            # TODO Throw some sort of exception if no name is provided?
            pass

        # Name is not guaranteed to be unique. There may be more than one
        # motif with the same name. In this case, return the first motif
        # fetched from the database and provide a warning.
        # This is the way it was handled in the TFBS perl modules, but is
        # this good practice? Is it pythonesque? If the motif are from
        # different collections perhaps we should prefer the one from
        # the default JASPAR collection?
        # TODO Decide how to handle this and implement it.
        sql = "select distinct BASE_ID from MATRIX where NAME = '%s'" % name

        cur = self.dbh.cursor()
        cur.execute(sql)

        rows = cur.fetchall()
        base_ids = []
        for row in rows:
            base_ids.append(row[0])

        num_motifs = len(base_ids)
        if num_motifs > 1:
            warn("There are %d JASPAR motifs with name '%s'" % (num_motifs, name))

        return self.fetch_motif_by_id(base_ids[0])

    def fetch_motif_set(
        self, collection=DFLT_JASPAR_COLLECTION, tf_name=None, tf_class=None,
        tf_family=None, matrix_id=None, tax_group=None, species=None,
        pazar_id=None, data_type=None, medline=None, min_ic=0, min_length=0,
        min_sites=0, all=False, all_versions=False
    ):
        """
        Fetch a set of motifs based on the provided selection criteria.

        Return a set of motifs.

        The the logic in the perl TFBS modules was that if no collection
        was specified, set it to default. But it is quite possible a user
        would want to select a set of motifs across different collectons.
        In the python context the collection argument could be explicitly set
        to None.
        TODO Think about this and make a decision.
        """

        motif_set = set()

        # Fetch the internal IDs of the motifs using the criteria provided
        int_ids = self._fetch_internal_id_list(
            collection = collection,
            tf_name = tf_name,
            tf_class = tf_class,
            tf_family = tf_family,
            matrix_id = matrix_id,
            tax_group = tax_group,
            species = species,
            pazar_id = pazar_id,
            data_type = data_type,
            medline = medline,
            all = all,
            all_versions = all_versions
        )

        for int_id in int_ids:
            motif = self._fetch_motif_by_internal_id(int_id)

            pfm = motif.counts

            # Filter motifs to those with matrix IC greater than min_ic
            if min_ic:
                # TODO create a total_ic() method. This should be a
                # method of the FrequencyPositionMatrix class
                if counts.total_ic() < min_ic:
                    continue

            # Filter motifs to those with minimum length of min_length
            if min_length:
                if motif.length < min_length:
                    continue

            # TODO - Not provided in perl TFBS modules but perhaps we should
            # also have a max_length filter.

            # Filter motifs to those composed of at least this many sites.
            # The perl TFBS module assumes column sums may be different but
            # this should be strictly enforced here we will ignore this and
            # just use the first colum sum.
            if min_sites:
                num_sites = sum(
                    [motif.counts[nt][0] for nt in motif.alphabet.letters]
                )
                if  num_sites < min_sites:
                    continue

            motif_set.add(motif)

        return motif_set

    def _fetch_latest_version(self, base_id):
        """
        Get the latest version number for the given base_id,

        """

        sql = "select VERSION from MATRIX where BASE_id = '%s' order by VERSION desc limit 1" % base_id

        cur = self.dbh.cursor()
        cur.execute(sql)

        latest = cur.fetchone()[0]

        return latest

    def _fetch_internal_id(self, base_id, version):
        """
        Fetch the internal id for a base id + version. Also checks if this
        combo exists or not
        """

        sql = "select id from MATRIX where BASE_id = '%s' and VERSION = '%s'" % (base_id, version)

        cur = self.dbh.cursor()
        cur.execute(sql)

        int_id = cur.fetchone()[0]

        return int_id


    def _fetch_motif_by_internal_id(self, int_id):
        # fetch basic motif information
        sql = "select BASE_ID, VERSION, COLLECTION, NAME from MATRIX where id = %d" % int_id

        cur = self.dbh.cursor()
        cur.execute(sql)

        row = cur.fetchone()

        base_id     = row[0]
        version     = row[1]
        collection  = row[2]
        name        = row[3]

        matrix_id = "".join([base_id, '.', str(version)])

        # fetch the counts matrix
        counts = self._fetch_counts_matrix(int_id)

        # Create new JASPAR motif
        motif = jaspar.Motif(
            matrix_id, name, collection = collection, counts = counts
        )

        # fetch species
        sql = "select TAX_ID from MATRIX_SPECIES where id = %d" % int_id
        cur.execute(sql)
        tax_ids = []
        rows = cur.fetchall()
        for row in rows:
            tax_ids.append(row[0])

        motif.species = tax_ids

        # fetch protein accession numbers
        sql = "select ACC FROM MATRIX_PROTEIN where id = %d" % int_id
        cur.execute(sql)
        accs = []
        rows = cur.fetchall()
        for row in rows:
            accs.append(row[0])

        motif.acc = accs

        # fetch remaining annotation as tags from the ANNOTATION table
        sql = "select TAG, VAL from MATRIX_ANNOTATION where id = %d" % int_id
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            attr = row[0]
            val  = row[1]
            if attr == 'class':
                motif.tf_class = val
            elif attr == 'family':
                motif.tf_family = val
            elif attr == 'tax_group':
                motif.tax_group = val
            elif attr == 'type':
                motif.data_type = val
            elif attr == 'pazar_tf_id':
                motif.pazar_id = val
            elif attr == 'medline':
                motif.medline = val
            elif attr == 'comment':
                motif.comment = val
            else:
                # TODO If we were to implement additional abitrary tags
                #motif.tag(attr, val)
                pass

        return motif

    def _fetch_counts_matrix(self, int_id):
        """
        Fetch the counts matrix from the JASPAR DB by the internal ID

        Returns a Bio.motifs.matrix.GenericPositionMatrix
        """
        counts = {}
        cur = self.dbh.cursor()

        for base in dna.letters:
            base_counts = []

            cur.execute("select val from MATRIX_DATA where ID = %s and row = %s order by col", (int_id, base))

            rows = cur.fetchall()
            for row in rows:
                base_counts.append(row[0])

            counts[base] = map(float, base_counts)

        return matrix.GenericPositionMatrix(dna, counts)

    def _fetch_internal_id_list(
        self, collection=DFLT_JASPAR_COLLECTION, tf_name=None, tf_class=None,
        tf_family=None, matrix_id=None, tax_group=None, species=None,
        pazar_id=None, data_type=None, medline=None, all=False,
        all_versions=False
    ):
        """
        Fetch a list of internal JASPAR motif IDs based on various passed
        parameters which may then be used to fetch the actual matrices.

        Caller: fetch_motif_set()

        1: First catch non-tag things like collection, name and version,
           species. Make one query for these if they are named and check the
           IDs for "latest" unless requested not to. These are AND statements.

        2: Then do the rest on tag level. To be able to do this with actual
           AND statemnet in the tag table, we do an inner join query, which
           is kept separate just for convenience

        3: Intersect 1 and 2

        For the surviving matrices, the responsibility to do matrix-based
        feature filtering such as ic, number of sites etc, fall on the
        calling fetch_motif_set() method.

        TODO
        This is a direct python translation of the perl TFBS package
        _get_IDlist_by_query and needs serious cleanup/refactoring.

        The the logic in the perl TFBS modules was that if no collection
        was specified, set it to default. But it is quite possible a user
        would want to select a set of motifs across different collectons.
        In the python context the collection argument could be explicitly set
        to None.
        TODO Think about this and make a decision.
        """

        int_ids = []

        cur = self.dbh.cursor()

        # Should redo so that matrix_annotation queries are separate, with an
        # intersect in the end 

        # Special case 1: fetch ALL motifs. Highest priority.
        # Ignore all other selection arguments.
        if all:
            cur.execute("select ID from MATRIX")
            rows = cur.fetchall()
            
            for row in rows:
                int_ids.append(row[0])

            return int_ids
                

        # Special case 2: fetch specific motifs by their JASPAR IDs. This
        # has higher priority than any other except the above 'all' case.
        # Ignore all other selection arguments.
        if matrix_id:
            # These might be either stable IDs or stable_ID.version.
            # If just stable ID and if all_versions == 1, return all versions,
            # otherwise just the latest
            if all_versions:
                for id in matrix_id:
                    # ignore vesion here, this is a stupidity filter
                    (base_id, version) = _split_jaspar_id(id) 
                    cur.execute(
                        "select ID from MATRIX where BASE_ID = %s", base_id
                    )

                    rows = cur.fetchall()
                    for row in rows:
                        int_ids.append(row[0])
            else:
                # only the lastest version, or the requested version
                for id in ids:
                    (base_id, version) = _split_jaspar_id(id) 

                    if not version:
                        version = self._fetch_latest_version(base_id)

                    int_id = self._fetch_internal_id(base_id, version)

                    if int_id:
                        int_ids.append(int_id)

            return int_ids

        tables = ["MATRIX m"]
        where_clauses = []

        # Select by MATRIX.COLLECTION
        if collection:
            if isinstance(collection, list):
                # Multiple collections passed in as a list
                clause = "m.COLLECTION in ('"
                clause = "".join([clause, "','".join(collection)])
                clause = "".join([clause, "')"])
            else:
                # A single collection - typical usage
                clause = "m.COLLECTION = '%s'" % collection

            where_clauses.append(clause)

        # Select by MATRIX.NAME
        if tf_name:
            if isinstance(tf_name, list):
                # Multiple names passed in as a list
                clause = "m.NAME in ('"
                clause = "".join([clause, "','".join(tf_name)])
                clause = "".join([clause, "')"])
            else:
                # A single name
                clause = "m.NAME = '%s'" % tf_name

            where_clauses.append(clause)

        # Select by MATRIX_SPECIES.TAX_ID
        if species:
            tables.append("MATRIX_SPECIES ms")
            where_clauses.append("m.ID = ms.ID")

            # NOTE: species are numeric taxonomy IDs but stored as varchars
            # in the DB.
            if isinstance(species, list):
                # Multiple tax IDs passed in as a list
                clause = "ms.TAX_ID in ('"
                clause = "".join([clause, "','".join(str(s) for s in species)])
                clause = "".join([clause, "')"])
            else:
                # A single tax ID
                clause = "ms.TAX_ID = '%s'" % str(species)

            where_clauses.append(clause)


        #
        # Tag based selection from MATRIX_ANNOTATION
        # Differs from perl TFBS module in that the matrix class explicitly
        # has a tag attribute corresponding to the tags in the database. This
        # provides tremendous flexibility in adding new tags to the DB and
        # being able to select based on those tags with out adding new code.
        # In the JASPAR Motif class we have elected to use specific attributes
        # for the most commonly used tags and here correspondingly only allow
        # selection on these attributes.
        #
        # The attributes corresponding to the tags for which selection is
        # provided are:
        #
        #   Attribute   Tag
        #   tf_class    class
        #   tf_family   family
        #   pazar_id    pazar_tf_id
        #   medline     medline
        #   data_type   type
        #   tax_group   tax_group
        #

        # Select by MATRIX_ANNOTATION VAL="class"
        if tf_class:
            tables.append("MATRIX_ANNOTATION ma1")
            where_clauses.append("m.ID = ma1.ID")

            clause = "ma1.TAG = 'class'"
            if isinstance(tf_class, list):
                clause = "".join([clause, " and ma1.VAL in ('"])
                clause = "".join([clause, "','".join(tf_class)])
                clause = "".join([clause, "')"])
            else:
                # A single tax ID
                clause = "".join([clause, " and ma1.VAL = '%s' " % tf_class])

            where_clauses.append(clause)

        # Select by MATRIX_ANNOTATION VAL="family"
        if tf_family:
            tables.append("MATRIX_ANNOTATION ma2")
            where_clauses.append("m.ID = ma2.ID")

            clause = "ma2.TAG = 'family'"
            if isinstance(tf_family, list):
                clause = "".join([clause, "ma2.VAL in ('"])
                clause = "".join([clause, "','".join(family)])
                clause = "".join([clause, "')"])
            else:
                # A single tax ID
                clause = "".join([clause, "ma2.VAL = '%s' " % tf_family])

            where_clauses.append(clause)

        # Select by MATRIX_ANNOTATION VAL="pazar_tf_id"
        if tf_family:
            tables.append("MATRIX_ANNOTATION ma3")
            where_clauses.append("m.ID = ma3.ID")

            clause = "ma3.TAG = 'pazar_tf_id'"
            if isinstance(pazar_id, list):
                clause = "".join([clause, "ma3.VAL in ('"])
                clause = "".join([clause, "','".join(pazar_id)])
                clause = "".join([clause, "')"])
            else:
                # A single tax ID
                clause = "".join(["ma3.VAL = '%s' " % pazar_id])

            where_clauses.append(clause)

        # Select by MATRIX_ANNOTATION VAL="medline"
        if medline:
            tables.append("MATRIX_ANNOTATION ma4")
            where_clauses.append("m.ID = ma4.ID")

            clause = "ma4.TAG = 'medline'"
            if isinstance(medline, list):
                clause = "".join([clause, "ma4.VAL in ('"])
                clause = "".join([clause, "','".join(medline)])
                clause = "".join([clause, "')"])
            else:
                # A single tax ID
                clause = "".join(["ma4.VAL = '%s' " % medline])

            where_clauses.append(clause)

        # Select by MATRIX_ANNOTATION VAL="type"
        if data_type:
            tables.append("MATRIX_ANNOTATION ma5")
            where_clauses.append("m.ID = ma5.ID")

            clause = "ma5.TAG = 'data_type'"
            if isinstance(data_type, list):
                clause = "".join([clause, "ma5.VAL in ('"])
                clause = "".join([clause, "','".join(data_type)])
                clause = "".join([clause, "')"])
            else:
                # A single tax ID
                clause = "".join(["ma5.VAL = '%s' " % data_type])

            where_clauses.append(clause)

        # Select by MATRIX_ANNOTATION VAL="tax_group"
        if tax_group:
            tables.append("MATRIX_ANNOTATION ma6")
            where_clauses.append("m.ID = ma6.ID")

            clause = "ma6.TAG = 'tax_group'"
            if isinstance(tax_group, list):
                clause = "".join([clause, " and ma6.VAL in ('"])
                clause = "".join([clause, "','".join(tax_group)])
                clause = "".join([clause, "')"])
            else:
                # A single tax ID
                clause = "".join([clause, " and ma6.VAL = '%s' " % tax_group])

            where_clauses.append(clause)

        sql = "".join(["select distinct(m.ID) from ", ", ".join(tables)])

        if where_clauses:
            sql = "".join([sql, " where ", " and ".join(where_clauses)])
        
        #print "sql = %s" % sql

        cur.execute(sql)
        rows = cur.fetchall()

        for row in rows:
            id = row[0]
            if all_versions:
                int_ids.append(id)
            else:
                # is the latest version?
                if self._is_latest_version(id):
                    int_ids.append(id)

        if len(int_ids) < 1:
            warn("Warning: Zero motifs returned with current select critera")

        return int_ids

    def _is_latest_version(self, int_id):
        # Does this internal ID represened the latest version of the JASPAR
        # matrix (collapse on base ids)
        cur = self.dbh.cursor()

        cur.execute("select count(*) from MATRIX where BASE_ID = (select BASE_ID from MATRIX where ID = %s) and VERSION > (select VERSION from MATRIX where ID = %s)", (int_id, int_id))

        row = cur.fetchone();

        count = row[0]

        if count == 0:
            # no matrices with higher version ID and same base id
            return(True)

        return(False)

def _split_jaspar_id(id):
    id_split = id.split('.')

    base_id = None
    version = None
    if len(id_split) == 2:
        base_id = id_split[0]
        version = id_split[1]
    else:
        base_id = id

    return (base_id, version)
