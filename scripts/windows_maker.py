import os
from dataclasses import dataclass
from typing import Optional
from pathlib import Path
import numpy as np
import pandas as pd
import logging

@dataclass
class WindowMaker:

    window_size: int = 500
    base: int = 1

    def __post_init__(self) -> None:
        logging.basicConfig(level=logging.WARNING, format="%(levelname)s:%(message)s")

    @staticmethod 
    def read_genome(genome: os.PathLike[str]) -> dict[str, int]:
        if not Path(genome).is_file():
            raise FileNotFoundError(f"Could not locate file {genome}.")
        genome_stats = {}
        with open(genome, mode="r", encoding="UTF-8") as f:
            for line in f:
                seqID, length = line.strip().split("\t")
                genome_stats.update({seqID: int(length)})
        return genome_stats

    def make_half_windows(self, df: pd.DataFrame, loci: str, genome: Optional[os.PathLike[str]] = None) -> pd.DataFrame:
        """"""
        df = df.copy()
        if genome is None:
            # print(f"User did not provide genome input. Some windows make be over the respective chromosome size.")
            genome_stats = {}
        else:
            genome_stats = WindowMaker.read_genome(genome)

        if self.base == 1:
            df["start"] = df["start"] - 1
            df["end"] = df["end"] - 1
        df["chromosomeSize"] = df["seqID"].apply(lambda seqID: genome_stats.get(seqID, np.nan))
        # if df['chromosomeSize'].isna().sum() > 0:
        # print(f"Provided chromosome sizes are incomplete. Some chromosome were not mapped to genome sizes.")
        df["chromosomeSize"] = df["chromosomeSize"].fillna(float('inf'))

        def _make_windows(df: pd.DataFrame, loci: str) -> pd.DataFrame:
            if loci != "start" and loci != "end":
                raise ValueError(f"Unknown loci {loci}.")
            if "chromosome" in df:
                df = df.rename(columns={"chromosome": "seqID"})
            if "seqID" not in df:
                raise KeyError("seqID column is not present in the dataframe.")
            if loci == "start":
                df.loc[:, "end"] = df[loci]
                df.loc[:, "start"] = np.maximum(df[loci] - self.window_size, 0)
            else:
                df.loc[:, "start"] = df[loci]
                df.loc[:, "end"] = np.minimum(df["chromosomeSize"], df[loci].copy() + self.window_size)
            columns = [col for col in df.columns if col != "start" and col != "end" and col != "seqID"]
            return df[["seqID", "start", "end"] + columns]

        expanded_df_positive = _make_windows(df.query("strand == '+'"), loci=loci)
        expanded_df_negative = _make_windows(df.query("strand == '-'"), loci="start" if loci == "end" else "end")
        df = pd.concat([expanded_df_positive, expanded_df_negative], axis=0)
        df.sort_index(ascending=True, inplace=True)
        return df.drop(columns=['chromosomeSize'])

    def make_windows(self, df: pd.DataFrame, loci: str, genome: Optional[os.PathLike[str]] = None) -> pd.DataFrame:
        """"""
        if genome is None:
            print(f"User did not provide genome input. Some windows make be over the respective chromosome size.")
            genome_stats = {}
        else:
            genome_stats = WindowMaker.read_genome(genome)

        df = df.copy()
        if self.base == 1:
            df["start"] = df["start"] - 1
            df["end"] = df["end"] - 1
        if loci != "mid":
            positive_df = df.query("strand != '-'")
            negative_df = (
                    df.query("strand == '-'")
                    .rename(
                        columns={"start": "end", 
                                 "end": "start"}
                        )
                    )
            df = pd.concat([positive_df, negative_df], axis=0)
        df["chromosomeSize"] = df["seqID"].apply(lambda seqID: genome_stats.get(seqID, np.nan))
        if df['chromosomeSize'].isna().sum() > 0:
            print(f"Provided chromosome sizes are incomplete. Some chromosome were not mapped to genome sizes.")
        df["chromosomeSize"] = df["chromosomeSize"].fillna(float('inf'))
        df.sort_index(ascending=True, inplace=True)

        def _make_windows(df: pd.DataFrame, loci: str) -> pd.DataFrame:
            if loci != "start" and loci != "end" and loci != "mid":
                raise ValueError(f"Unknown loci {loci}.")
            if "chromosome" in df:
                df = df.rename(columns={"chromosome": "seqID"})
            if "seqID" not in df:
                raise KeyError("seqID column is not present in the dataframe.")

            if loci == "start":
                df["end"] = np.minimum(df["chromosomeSize"], df[loci] + self.window_size + 1)
                df["start"] = np.maximum(df[loci] - self.window_size, 0)
            elif loci == "end":
                df["start"] = np.maximum(df[loci] - self.window_size, 0)
                df["end"] = np.minimum(df["chromosomeSize"], df[loci] + self.window_size + 1)
            else:
                df["start"] = (df["start"] + df["end"]) // 2
                df["end"] = np.minimum(df["chromosomeSize"], df["start"] + self.window_size + 1)
                df["start"] = np.maximum(df["start"] - self.window_size, 0)
            columns = [col for col in df.columns if col != "start" and col != "end" and col != "seqID"]
            df["start"] = df["start"].astype(int)
            df["end"] = df["end"].astype(int)
            return df[["seqID", "start", "end"] + columns]
        return _make_windows(df=df, loci=loci).drop(columns=["chromosomeSize"])
