import MySQLdb as mdb

from Bio.Alphabet.IUPAC import unambiguous_dna as dna

import sys
sys.path.append("/homed/home/dave/devel/biopython-1.61")
#sys.path.append("/home/anthony/PostDoc/JASPAR2013/Biopython_package/biopython-master/")
from Bio_dev.motifs import jaspar, matrix


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

    def fetch_profile_by_id(self, id):
        """
        Fetch a single JASPAR profile from the DB by it's JASPAR matrix id
        (e.g. 'MA0001.1').

        NOTE: The perl TFBS module allows you to specify the type of matrix to
        return (PFM, PWM, ICM) but matrices are always stored in JASAPR as
        PFMs so this does not really belong here. Once a PFM is fetched the
        pwm() and pssm() methods can be called to return the normalized and
        log-odds matrices.
        """
         
        if not id:
            # TODO Throw some sort of exception if no id is provided
            pass

        # separate stable ID and version number
        id_parts = id.split('.') 
        base_id = id_parts[0]
        version = None
        if len(id_parts) == 2:
            version = int(id_parts[1])
        else:
            # if ID contains no version portion, fetch latest version by default
            version = self._fetch_latest_version(base_id)

        #print "base_id = %s" % base_id
        #print "version = %d" % version

        # fetch internal JASPAR profile ID - also a check for validity
        int_id = self._fetch_internal_id(base_id, version)

        #print "int_id = %d" % int_id

        # fetch JASPAR profile using internal ID
        profile = self._fetch_profile_by_internal_id(int_id)

        return profile

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

        sql = "select id from MATRIX where BASE_id = '%s' and VERSION = '%s'" % (base_id, version);

        cur = self.dbh.cursor()
        cur.execute(sql)

        int_id = cur.fetchone()[0]

        return int_id


    def _fetch_profile_by_internal_id(self, int_id):

        # fetch basic profile information
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

        # Create new JASPAR profile
        profile = jaspar.Motif(
            matrix_id, name, collection = collection, counts = counts
        )

        # fetch species
        sql = "select TAX_ID from MATRIX_SPECIES where id = %d" % int_id
        cur.execute(sql)
        tax_ids = []
        rows = cur.fetchall()
        for row in rows:
            tax_ids.append(row[0])

        profile.species = tax_ids

        # fetch protein accession numbers
        sql = "select ACC FROM MATRIX_PROTEIN where id = %d" % int_id;
        cur.execute(sql)
        accs = []
        rows = cur.fetchall()
        for row in rows:
            accs.append(row[0])

        profile.acc = accs

        # fetch remaining annotation as tags from the ANNOTATION table
        sql = "select TAG, VAL from MATRIX_ANNOTATION where id = %d" % int_id
        cur.execute(sql)
        rows = cur.fetchall()
        for row in rows:
            attr = row[0]
            val  = row[1]
            if attr == 'class':
                profile.tf_class = val
            elif attr == 'family':
                profile.tf_family = val
            elif attr == 'tax_group':
                profile.tax_group = val
            elif attr == 'type':
                profile.data_type = val
            elif attr == 'pazar_tf_id':
                profile.pazar_id = val
            elif attr == 'medline':
                profile.medline = val
            elif attr == 'comment':
                profile.comment = val
            else:
                # TODO If we were to implement additional abitrary tags
                #profile.tag(attr, val)
                pass

        return profile

    def _fetch_counts_matrix(self, int_id):
        """
        Fetch the counts matrix from the JASPAR DB by the internal ID

        Returns a Bio.motifs.matrix.GenericPositionMatrix
        """

        counts = {}
        cur = self.dbh.cursor()

        #print "int_id = %d" % int_id

        for base in dna.letters:
            #print "base = '%s'" % base

            sql = "select val from MATRIX_DATA where ID = %d and row = '%s' order by col" % (int_id, base)

            base_counts = []

            cur.execute(sql)
            rows = cur.fetchall()
            for row in rows:
                base_counts.append(row[0])

            counts[base] = map(float, base_counts)

        return matrix.GenericPositionMatrix(dna, counts)
