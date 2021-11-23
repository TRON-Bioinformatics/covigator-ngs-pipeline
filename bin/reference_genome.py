from pysam import FastaFile


class ReferenceGenomeReader(object):

    def __init__(self, fasta_file_path):
        """
        :type fasta_file_path: str
        """
        self.reference = FastaFile(filename=fasta_file_path)

    def read_position_from_reference(self, chromosome, position, length=1):
        """
        :type chromosome: str
        :param position: 1-based position
        :type position: int
        :rtype: str
        """
        try:
            base = self.reference.fetch(reference=chromosome, start=position-1, end=position + length - 1)
        except Exception as ex:
            raise ValueError("Provided chromosome and position ({}-{}) returned no base".format(
                chromosome, position, str(ex)))
        if len(base) == 0:
            raise ValueError("Provided chromosome and position ({}-{}) returned no base: ".format(
                chromosome, position))
        return base.upper()
