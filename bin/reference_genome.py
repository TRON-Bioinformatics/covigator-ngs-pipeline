from pysam import FastaFile


class ReferenceGenomeReader(object):

    def __init__(self, fasta_file_path):
        """
        :type fasta_file_path: str
        """
        self.reference = FastaFile(filename=fasta_file_path)
        assert self.reference.nreferences == 1, "A reference with more than one chromosome is provided"
        self.chromosome = self.reference.references[0]

    def read_position_from_reference(self, position, length=1):
        """
        :param position: 1-based position
        :type position: int
        :rtype: str
        """
        try:
            base = self.reference.fetch(reference=self.chromosome, start=position-1, end=position + length - 1)
        except Exception as ex:
            raise ValueError("Provided chromosome and position ({}-{}) returned no base".format(
                self.chromosome, position, str(ex)))
        if len(base) == 0:
            raise ValueError("Provided chromosome and position ({}-{}) returned no base: ".format(
                self.chromosome, position))
        return base.upper()
