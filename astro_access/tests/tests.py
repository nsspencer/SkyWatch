import unittest
from astro_access import FrameInterpolator
import csv
import os


LEO_EPHEM_CSV = os.path.join(os.path.dirname(__file__), "data", "leo_sat_ephem.csv")


def read_csv(path:str) -> dict:
    # Open the CSV file in read mode
    with open(path, 'r') as file:
        # Create a CSV reader object
        csv_reader = csv.reader(file)
        
        # Iterate over each row in the CSV file
        rows = []
        for row in csv_reader:
            # Each row is a list of values representing columns
            rows.append(row)
    return rows

def get_ephem_data(path: str = LEO_EPHEM_CSV) -> dict:
    data = {}
    csv_data = read_csv(path)
    header = csv_data[0]
    
    column_data = [[] for _ in header]
    for row in csv_data[1:]:
        for index, column in enumerate(row):
            column_data[index].append(column)
    
    for heading, column in zip(header, column_data):
        data[heading] = column
        
    return data

class EphemTest(unittest.TestCase):
    
    def test_ephemeris_interpolation(self):
        ephemeis_csv = read_csv(LEO_EPHEM_CSV)
        pass


def run():
    unittest.main()