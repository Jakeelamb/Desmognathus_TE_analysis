#!/usr/bin/env python3
"""
Unit tests for the superfamily proportions script.
"""

import os
import sys
import unittest
import pandas as pd
import tempfile
from pathlib import Path

# Add parent directory to path to import from the source code
parent_dir = str(Path(__file__).resolve().parent.parent)
sys.path.insert(0, parent_dir)

# Now import from the preprocessing module
from preprocessing.generate_superfamily_proportions import create_superfamily_proportions


class TestSuperfamilyProportions(unittest.TestCase):
    """Test cases for the superfamily proportions calculation."""

    def setUp(self):
        """Set up test data."""
        # Use an absolute path to find the test data
        project_root = Path(__file__).resolve().parent.parent.parent.parent
        self.test_data_path = project_root / "data/test_data/raw/sample_te_counts.csv"
        
        # Create a temporary input file with the expected format
        self.temp_dir = tempfile.mkdtemp()
        self.temp_input = os.path.join(self.temp_dir, "superfamily_breakdown.csv")
        self.temp_output = os.path.join(self.temp_dir, "superfamily_proportions.csv")
        
        # Create a sample breakdown CSV
        df = pd.DataFrame({
            'superfamily': ['DNA/hAT', 'LINE/BovB', 'LINE/L1', 'LTR/Gypsy'],
            'Desmognathus_aeneus': [142, 103, 254, 187],
            'Desmognathus_ocoee': [157, 92, 301, 162],
            'Desmognathus_orestes': [148, 98, 287, 176]
        })
        df.to_csv(self.temp_input, index=False)
    
    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)
        
    def test_create_superfamily_proportions(self):
        """Test the superfamily proportions creation function."""
        # Run the function
        result = create_superfamily_proportions(self.temp_input, self.temp_output)
        
        # Check that it succeeded
        self.assertTrue(result)
        
        # Read the output file
        result_df = pd.read_csv(self.temp_output)
        
        # Check that the result has expected columns and rows
        expected_columns = ['superfamily', 'Desmognathus_aeneus', 
                          'Desmognathus_ocoee', 'Desmognathus_orestes']
        self.assertListEqual(result_df.columns.tolist(), expected_columns)
        
        # Check that we have the right number of rows (one per superfamily)
        self.assertEqual(len(result_df), 4)
        
        # Verify each superfamily is present
        superfamilies = result_df['superfamily'].tolist()
        self.assertIn('DNA/hAT', superfamilies)
        self.assertIn('LINE/BovB', superfamilies)
        self.assertIn('LINE/L1', superfamilies)
        self.assertIn('LTR/Gypsy', superfamilies)
        
        # Check a few specific values
        # For Desmognathus_aeneus, proportions should be:
        # - DNA/hAT: 142 / (142 + 103 + 254 + 187) = 142 / 686 ≈ 0.207
        # - LINE/BovB: 103 / 686 ≈ 0.150
        # - LINE/L1: 254 / 686 ≈ 0.370
        # - LTR/Gypsy: 187 / 686 ≈ 0.273
        dna_row = result_df[result_df['superfamily'] == 'DNA/hAT'].iloc[0]
        self.assertAlmostEqual(dna_row['Desmognathus_aeneus'], 0.207, places=3)
        
        line_bovb_row = result_df[result_df['superfamily'] == 'LINE/BovB'].iloc[0]
        self.assertAlmostEqual(line_bovb_row['Desmognathus_aeneus'], 0.150, places=3)
        
        # Verify that proportions sum to 1 for each species
        for species in ['Desmognathus_aeneus', 'Desmognathus_ocoee', 'Desmognathus_orestes']:
            proportion_sum = result_df[species].sum()
            self.assertAlmostEqual(proportion_sum, 1.0, places=5)


if __name__ == '__main__':
    unittest.main() 