# -*- coding: utf-8 -*-
"""
Elite Data Hacks

@author: Christopher S. Josef, MD
@email: csjosef@krvmail.com
"""

import csv
from glob import glob

# enter folder with files that you want to: 1. Strip unecessary quotes 2. Make all headers lower case
path = r"C:/Users/DataSci/OneDrive - Emory University/Sepsis Calculation/Data_Export/2014/*.dsv"  # Note wildcard character added for glob().

for fname in glob(path):
    print(fname)
    with open(fname, newline="") as f:
        reader = csv.reader(
            f,
            delimiter="|",
            quotechar='"',
        )
        header = next(reader)  # Get the header row.
        header = [column.lower() for column in header]  # Lowercase the headings.
        rows = [header] + list(reader)  # Read the rest of the rows.

    with open(fname, "w", newline="") as f:
        writer = csv.writer(
            f, delimiter="|", quotechar='"', escapechar="\\", quoting=csv.QUOTE_NONE
        )
        writer.writerows(rows)  # Write new header & original rows back to file.
