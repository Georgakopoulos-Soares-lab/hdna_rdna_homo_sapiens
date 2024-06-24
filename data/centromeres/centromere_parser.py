import pandas as pd
from abc import abstractmethod
import subprocess
import os
from pathlib import Path
import tempfile

class CentromereParser:

    def __init__(self, accession):
        self.accession = Path(accession).resolve()

    @property
    def species(self):
        return '.'.join(self.accession.name.split('.')[:-1])

    @abstractmethod
    def parse(self):
        pass

    def extract_coverage(self):

        command = f"{self.bedtools} coverage -a {}"

class HomoSapiensCen(CentromereParser):

    def parse(self):
        df = pd.read_csv(self.accession,
                         delimiter="\t",
                         skiprows=1,
                         header=None,
                         usecols=range(4),
                         names=["seq_id", "start", "end", "compartment"],
                         )
        df.loc[:, "compartment"] = df["compartment"].str.split("_", expand=True)[0]
        return df

class PrimateCen(CentromereParser):

    def parse(self):
        df = pd.read_csv(self.accession,
                         delimiter="\t",
                         skiprows=1,
                         header=None,
                         usecols=range(4),
                         names=["seq_id", "start", "end", "compartment"],
                         )
        df.loc[:, "compartment"] = df["compartment"].str.split("(", expand=True)[0]

        lower = ['gSat', 'HSat3', 'bSat', 'cenSat', 'HSat2']
        for col in lower:
            df['compartment'] = df['compartment'].replace(col, col.lower())
        return df



if __name__ == "__main__":

    parser = HomoSapiensCen("chm13v2.0_censat_v2.0.bed")
    df = parser.parse()

    parser = PrimateCen("mGorGor1_v2.0_CenSat_v0.1.bed")
    dft = parser.parse()
    breakpoint()
