#!/usr/bin/env python3

import unittest
import pathlib
import pfam_scan as ps

class Test(unittest.TestCase):

    def setUp(self):
        self.dir = pathlib.Path(__file__).parent / 'test'

    def test_read_pfam_data(self):
        data_path = self.dir / 'test.hmm.dat'
        data = ps.read_pfam_data(data_path)
        self.assertEqual(len(data), 24)
        self.assertEqual(data['Dicer_dimer'].type, 'Domain')
        self.assertEqual(data['Dicer_dimer'].clan, 'CL0196')
        self.assertEqual(data['Dicer_dimer'].ga_dom, 22.6)
        self.assertEqual(data['Dicer_dimer'].ga_seq, 22.6)
        self.assertIsNone(data['ArgoN'].clan, None)

    def test_parse_hmmscan_output(self):
        data_path = self.dir / 'test.hmm.dat'
        data = ps.read_pfam_data(data_path)
        output_path = self.dir / 'test.domtblout.txt'
        results = ps.parse_hmmscan_output(output_path, data)
        self.assertEqual(len(results['Q9SP32']), 10)
        self.assertEqual(len(results['O04379']), 7)
        self.assertEqual(len(results['O43347']), 2)
        domains = results['O04379']
        self.assertEqual(domains[0].hmm_acc, 'PF12764.10')
        self.assertEqual(domains[0].hmm_name, 'Gly-rich_Ago1')
        self.assertEqual(domains[0].evalue, 6.9e-24)
        self.assertEqual(domains[6].hmm_acc, 'PF02171.20')
        self.assertEqual(domains[6].hmm_name, 'Piwi')
        self.assertEqual(domains[6].aln_start, 677)
        self.assertEqual(domains[6].aln_end, 996)

    def test_resolve_overlapping_domains1(self):
        results = {
            'E0SP36': [
                ps.Domain(
                    63, 188, 58, 195, 'PF13401.9', 'AAA_22', 
                    'Domain', 6, 129, 137, 32.9, 7.1e-08, 1, 'CL0023'
                ), 
                ps.Domain(
                    65, 204, 65, 215, 'PF00004.32', 'AAA',
                    'Domain', 1, 120, 132, 33.3, 5.9e-08, 1, 'CL0023'
                ),
                ps.Domain(
                    312, 391, 312, 392, 'PF09079.14', 'Cdc6_C',
                    'Domain', 1, 83, 84, 83.1, 1.1e-23, 1, 'CL0123'
                )
            ]
        }
        results = ps.resolve_overlapping_domains(results)
        self.assertEqual(len(results['E0SP36']), 2)
        self.assertEqual(results['E0SP36'][0].hmm_name, 'AAA')
        self.assertEqual(results['E0SP36'][1].hmm_name, 'Cdc6_C')

    def test_resolve_overlapping_domains2(self):
        results = {
            'seq_id': [
                ps.Domain(
                    63, 188, 58, 195, 'PF13401.9', 'AAA_22', 
                    'Domain', 6, 129, 137, 32.9, 7.1e-08, 1, 'CL0023'
                ), 
                ps.Domain(
                    65, 204, 65, 215, 'PF00004.32', 'AAA',
                    'Domain', 1, 120, 132, 33.3, 5.9e-08, 1, 'CL0023'
                ),
                ps.Domain(
                    69, 210, 69, 210, 'PF00004.32', 'AAA_11',
                    'Domain', 1, 120, 132, 33.3, 2.1e-10, 1, 'CL0023'
                ),
            ]
        }
        results = ps.resolve_overlapping_domains(results)
        self.assertEqual(len(results['seq_id']), 1)
        self.assertEqual(results['seq_id'][0].hmm_name, 'AAA_11')

    def test_resolve_overlapping_domains3(self):
        results = {
            'seq_id': [
                ps.Domain(
                    63, 188, 58, 195, 'PF13401.9', 'AAA_22', 
                    'Domain', 6, 129, 137, 32.9, 7.1e-08, 1, 'CL0023'
                ), 
                ps.Domain(
                    65, 204, 65, 215, 'PF00004.32', 'AAA',
                    'Domain', 1, 120, 132, 33.3, 5.9e-08, 1, 'CL0023'
                ),
                ps.Domain(
                    69, 210, 69, 210, 'PF00004.32', 'AAA_11',
                    'Domain', 1, 120, 132, 33.3, 2.1e-10, 1, None
                ),
            ]
        }
        results = ps.resolve_overlapping_domains(results)
        self.assertEqual(len(results['seq_id']), 2)
        self.assertEqual(results['seq_id'][0].hmm_name, 'AAA')
        self.assertEqual(results['seq_id'][1].hmm_name, 'AAA_11')


if __name__ == '__main__':
    unittest.main()