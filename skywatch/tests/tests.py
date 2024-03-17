import csv
import os

import astropy.units as u
import numpy as np
from astropy.time import Time

from skywatch.coordinates import SkyPath

LEO_EPHEM_CSV = os.path.join(os.path.dirname(__file__), "data", "leo_sat_ephem.csv")


def read_csv(path: str) -> dict:
    # Open the CSV file in read mode
    with open(path, "r") as file:
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


def get_ephem_as_skypath(path: str = LEO_EPHEM_CSV) -> SkyPath:
    leo_csv = get_ephem_data(path)
    leo_csv_times = Time(leo_csv["utc"])
    leo_csv_position = (
        np.array(
            [leo_csv["x_eci_km"], leo_csv["y_eci_km"], leo_csv["z_eci_km"]]
        ).astype(float)
        * u.km
    )
    leo_csv_velocities = np.array(
        [leo_csv["vx_eci_km_s"], leo_csv["vy_eci_km_s"], leo_csv["vz_eci_km_s"]]
    ).astype(float) * (u.km / u.s)

    sat_position = SkyPath.from_ECI(
        leo_csv_times, *leo_csv_position, leo_csv_velocities
    )
    return sat_position
