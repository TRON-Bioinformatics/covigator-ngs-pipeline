import unittest
from bin.tests.test_phasing import TestPhasing
from bin.tests.test_assembly_variant_caller import TestCountryParser
from bin.tests.test_reference_genome import TestReferenceGenomeReader


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(TestPhasing())
    suite.addTest(TestCountryParser())
    suite.addTest(TestReferenceGenomeReader())
    result = unittest.TextTestRunner(verbosity=2).run(suite)

    if result.wasSuccessful():
        exit(0)
    else:
        exit(1)
