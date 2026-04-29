#!/usr/bin/env python3
"""
ncbi_taxonomy.py
Usage: python ncbi_taxonomy.py accessions.txt > output.tsv
       OR edit ACCESSIONS list directly
"""

import sys
import time
from Bio import Entrez

Entrez.email = "email.address@university.edu"  # Required by NCBI

ACCESSIONS = [
    "GCA_000507305.1",
    "GCA_001939145.1",
    "GCA_003297005.1",
    "GCA_003297045.1",
    "GCA_009767595.1",
    "GCA_018327485.1",
    "GCA_905221605.1",
    "GCA_905221615.1",
    "GCA_905221625.1",
    "GCA_905221635.1",
    "GCA_905231905.1",
    "GCA_905231915.1",
    "GCA_905231925.1",
    "GCA_947184155.1",
    "GCA_001417965.1",
    "GCA_024322265.1",
    "GCF_001417965.1",
    "GCA_022179025.1",
    "GCA_900290455.1",
    "GCA_942486025.1",
    "GCA_942486035.1",
    "GCA_958299795.1",
    "GCA_958299805.1",
    "GCA_003704095.1",
    "GCA_014529365.1",
    "GCA_942486045.1",
    "GCF_003704095.1",
    "GCA_000222465.2",
    "GCA_004143615.1",
    "GCA_013753865.1",
    "GCA_014633955.1",
    "GCA_014633975.1",
    "GCA_014634005.1",
    "GCA_014634045.1",
    "GCA_014634065.1",
    "GCA_014634105.1",
    "GCA_014634125.1",
    "GCA_014634145.1",
    "GCA_014634165.1",
    "GCA_014634205.1",
    "GCA_014634225.1",
    "GCA_014634525.1",
    "GCA_014634545.1",
    "GCA_014634585.1",
    "GCA_014634605.1",
    "GCA_020536085.1",
    "GCA_025960835.1",
    "GCF_000222465.1",
    "GCF_004143615.1",
    "GCF_013753865.1"
]

def fetch_assembly_taxonomy(accessions):
    print("accession\tgenus\tspecies")
    for acc in accessions:
        # Search assembly DB to get the UID
        handle = Entrez.esearch(db="assembly", term=acc)
        record = Entrez.read(handle)
        handle.close()

        uid = record["IdList"][0]

        # Fetch assembly summary (returns JSON-like XML)
        handle = Entrez.esummary(db="assembly", id=uid, report="full")
        summary = Entrez.read(handle, validate=False)
        handle.close()

        organism = summary["DocumentSummarySet"]["DocumentSummary"][0]["SpeciesName"]
        parts = organism.split()
        genus = parts[0] if len(parts) >= 1 else "NA"
        species = parts[1] if len(parts) >= 2 else "NA"
        print(f"{acc}\t{genus}\t{species}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            accs = [line.strip() for line in f if line.strip()]
    else:
        accs = ACCESSIONS
    fetch_assembly_taxonomy(accs)
