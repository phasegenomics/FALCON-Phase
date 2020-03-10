#!/usr/bin/env python

# tests for preprocessing script

import os
import unittest
import bin.preprocess_diploid_asm_for_fc_phase as pp

DIR = os.path.dirname(os.path.abspath(__file__))

class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.h_tig_file = os.path.join(DIR, "collateral/h_tigs_test.fasta")
        self.p_tig_file = os.path.join(DIR, "collateral/p_tigs_test.fasta")
        self.both_file = os.path.join(DIR, "collateral/combined.fasta")

    def tearDown(self):
        pass

    def test_scrub_name(self):
        scrubbed = pp.scrub_name("000001F|arrow")
        self.assertEqual(scrubbed, "000001F")
        scrubbed = pp.scrub_name("000001F|quiver")
        self.assertEqual(scrubbed, "000001F")
        scrubbed = pp.scrub_name("000001F_000001F|fake")
        self.assertEqual(scrubbed, "000001F_000001")

    def test_is_p_tig(self):
        self.assertTrue(pp.is_p_contig("000001F"))
        self.assertFalse(pp.is_p_contig("000001F_001"))
        self.assertTrue(pp.is_p_contig("0"))
        self.assertFalse(pp.is_p_contig("0_0"))

    def test_is_h_tig(self):
        self.assertFalse(pp.is_h_contig("000001F"))
        self.assertTrue(pp.is_h_contig("000001F_001"))
        self.assertFalse(pp.is_h_contig("0"))
        self.assertTrue(pp.is_h_contig("0_0"))

    def test_parse_separated_fasta(self):
        p_tigs = pp.parse_separated_fasta(self.p_tig_file)
        h_tigs = pp.parse_separated_fasta(self.h_tig_file)
        p_names = sorted(p_tig.id for p_tig in p_tigs)
        h_names = sorted(h_tig.id for h_tig in h_tigs)
        self.assertEqual(p_names, ["000001F", "000002F", "000003F", "000004F"])
        self.assertEqual(h_names, ["000001F_001", "000002F_001", "000003F_001", "000003F_002"])

    def test_parse_combined_fasta(self):
        p_tigs, h_tigs = pp.parse_combined_fasta(self.both_file)
        p_names = sorted(p_tig.id for p_tig in p_tigs)
        h_names = sorted(h_tig.id for h_tig in h_tigs)
        self.assertEqual(p_names, ["000001F", "000002F", "000003F", "000004F"])
        self.assertEqual(h_names, ["000001F_001", "000002F_001", "000003F_001", "000003F_002"])

    def test_name_mappings(self):
        p_tigs, h_tigs = pp.parse_combined_fasta(self.both_file)
        name_mappings, omissions = pp.make_name_mappings(p_tigs, h_tigs)
        self.assertEqual(["000001F_001"], name_mappings["000001F"])
        self.assertEqual(["000003F_001", "000003F_002"], name_mappings["000003F"])
        self.assertEqual([], name_mappings["000004F"])

if __name__ == '__main__':
    unittest.main()
