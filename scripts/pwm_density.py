import numpy as np
import pandas as pd
from typing import Callable, Optional
from abc import abstractmethod 
from tqdm import tqdm
import matplotlib.pyplot as plt
from seaborn import color_palette

class PWMExtractor:

    def __init__(self) -> None:
        self.nucleotides = "agct"

    @staticmethod
    def invert(nucleotide: str) -> str:
        match nucleotide:
            case 'a':
                return 't'
            case 't':
                return 'a'
            case 'g':
                return 'c'
            case 'c':
                return 'g'
            case _:
                raise ValueError(f'Unknown nucleotide {nucleotide}.')

    def determine_template_strand(self, compartment_strand: str, motif_strand: str) -> str:
        if compartment_strand == motif_strand:
            return "non_template"
        return "template"

    def extract_template_density(self, intersect_df: pd.DataFrame, 
                                        window_size: int, 
                                        enrichment: bool = False,
                                 ) -> pd.DataFrame:
        total_counts = {
                        "template": np.zeros(2*window_size+1),
                        "non_template": np.zeros(2*window_size+1)
                    }

        total_overlap = 0
        for _, row in intersect_df.iterrows():
                compartment_strand = row['strand']
                if compartment_strand == "?":
                    continue
                start = int(row["start"])
                end = int(row["end"])
                motif_start = int(row["motif_start"])
                motif_end = int(row["motif_end"])
                motif_strand = row["motif_strand"]

                template = self.determine_template_strand(compartment_strand, motif_strand)
                temp_counts = np.zeros(2*window_size+1)

                overlap = int(row['overlap'])
                total_overlap += overlap
                origin = end - window_size - 1
                L = max(0, window_size - (origin - motif_start))
                U = min(2 * window_size + 1, window_size - (origin - motif_end))
                assert L <= U
                overlap_start = max(motif_start, start)
                overlap_end = min(motif_end, end)
                overlap_length = overlap_end - overlap_start
                assert overlap == overlap_length

                temp_counts[L:U] += 1
                if compartment_strand == "-":
                    temp_counts = temp_counts[::-1]
                total_counts[template] += temp_counts

        # honour contract
        assert total_overlap == np.sum(total_counts["template"]) + np.sum(total_counts["non_template"]), f"Overlap: {total_overlap} vs. Calculated overlap {np.sum(total_counts)}."
        if enrichment:
            name = "Enrichment"
            for typ in total_counts:
                total_counts[typ] = total_counts[typ] / np.mean(total_counts[typ])
        else:
            name = "Occurrences"
        
        for typ in total_counts:
            total_counts[typ] = pd.Series(total_counts[typ])\
                                .to_frame(name=name + f"_{typ}")
            total_counts[typ].index = total_counts[typ].index - window_size 
            if name != "Enrichment":
                total_counts[typ] = total_counts[typ].astype(int)

        template_df = pd.concat([total_counts["template"],
                                 total_counts["non_template"]
                                 ], axis=1).T
        return template_df

    @staticmethod
    def bayes_estimator(counts: int, total_counts: int, total_bins: int) -> float:
        return (counts+1)/(total_counts+total_bins)

    @staticmethod
    def expected_value(counts: int, total_counts: int, total_bins: int) -> float:
        return PWMExtractor.bayes_estimator(counts, total_counts, total_bins)

    @staticmethod
    def enrichment(counts: int, total_counts: int, total_bins: int) -> float:
        return counts/(total_counts * PWMExtractor.bayes_estimator(counts, total_counts, total_bins))

    def get_relative_positions(self, 
                               intersect_df: pd.DataFrame, 
                               window_size: int) -> pd.DataFrame:
        relative_df = []
        total_overlap = 0
        for _, row in intersect_df.iterrows():
                compartment_strand = row['strand']
                if compartment_strand == "?":
                    continue
                start = int(row['start'])
                end = int(row['end'])
                motif_start = int(row['motif_start'])
                motif_end = int(row['motif_end'])

                temp_counts = np.zeros(2*window_size+1)
                overlap = int(row['overlap'])
                total_overlap += overlap
                origin = end - window_size - 1
                L = max(0, window_size - (origin - motif_start))
                U = min(2 * window_size + 1, window_size - (origin - motif_end))

                assert L <= U
                overlap_start = max(motif_start, start)
                overlap_end = min(motif_end, end)
                overlap_length = overlap_end - overlap_start
                assert overlap == overlap_length
                assert overlap == U-L
                temp_counts[L:U] += 1

                if compartment_strand == "-":
                    temp_counts = temp_counts[::-1]
                relative_df.append(list(temp_counts))
        relative_df = pd.DataFrame(relative_df, 
                                   columns=range(-window_size, window_size+1))
        total_sum = float(relative_df.values.sum())
        assert total_overlap == total_sum, f"Overlap: {total_overlap} vs. Calculated overlap {total_sum}."
        return relative_df

    def bootstrap(self, relative_positions: pd.DataFrame, 
                  N: int = 1_000,
                  lower_quantile: float = 0.025,
                  upper_quantile: float = 0.975) -> tuple[pd.Series, pd.Series, pd.Series]:
        bootstrap_df = []
        total_samples = relative_positions.shape[0]
        for _ in tqdm(range(N), leave=True, position=0):
            density = relative_positions.sample(total_samples, replace=True)\
                                        .sum(axis=0)
            density /= np.mean(density)
            bootstrap_df.append(density)
        bootstrap_df = pd.DataFrame(bootstrap_df)
        average = bootstrap_df.mean()
        lower_bound = bootstrap_df.quantile(lower_quantile)
        upper_bound = bootstrap_df.quantile(upper_quantile)
        return average, lower_bound, upper_bound
            
    def plot_density(self, 
                     density, 
                     lower_quantile: Optional[pd.Series] = None,
                     upper_quantile: Optional[pd.Series] = None,
                     xlabel: str = '',
                     ylabel: str = 'Enrichment',
                     title: str = '',
                     xlabel_size: int = 14,
                     ylabel_size: int = 14,
                     title_size: int = 16,
                     alpha: float = 1.0,
                     lw: float = 1.5,
                     color: str = 'black',
                     linestyle: str = '',
                     style: Optional[dict] = None):
        default_style = {
                        'horizontal': {'linestyle': '-',
                                        'lw': 1.0,
                                        'color': 'crimson',
                                        'alpha': 1.0,
                                        },
                        'vertical': {'linestyle': '--',
                                    'lw': 1.0,
                                    'color': 'black',
                                    'alpha': 1.0
                                    },
                        'confidence': 
                                {'color': color_palette("Set3")[3],
                                 'alpha': 0.4
                                 }
                         }
        # set styling attributes; replace with default
        if isinstance(style, dict):
            for key in default_style:
                if key not in style:
                    style[key] = {}
                for _property in default_style[key]:
                    if _property not in style[key]:
                        style[key][_property] = default_style[_property]
        else:
            style = default_style
                
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
        window_size = (len(density) - 1) // 2
        ax.plot(range(-window_size, window_size+1),
                density,
                lw=lw,
                color=color,
                linestyle=linestyle,
                zorder=3,
                alpha=alpha,
                )
        if isinstance(lower_quantile, pd.Series) and isinstance(upper_quantile, pd.Series):
            ax.fill_between(range(-window_size, window_size+1),
                        y1=lower_quantile,
                        y2=upper_quantile,
                        color=style['confidence']['color'],
                        alpha=style['confidence']['alpha'],
                        zorder=2
                        )
        ax.axhline(1.0, 
                   linestyle=style['horizontal']['linestyle'],
                   lw=style['horizontal']['lw'],
                   color=style['horizontal']['color'],
                   alpha=style['horizontal']['alpha']
                   )
        ax.axvline(0.0,
                   linestyle=style['vertical']['linestyle'],
                   lw=style['vertical']['lw'],
                   color=style['vertical']['color'],
                   alpha=style['vertical']['alpha']
                   )
        ax.grid(lw=0.4, alpha=0.6, zorder=0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.xaxis.label.set_size(xlabel_size)
        ax.yaxis.label.set_size(ylabel_size)
        ax.legend(loc=0, prop={'size': 12}, title='')
        ax.title.set_size(title_size)
        ax.set_title(title)
        ax.tick_params(axis="both", labelsize=13)
        return fig, ax
        
    def extract_density(self, intersect_df: pd.DataFrame, 
                        window_size: int,
                        return_array: bool = True,
                        return_frame: bool = False,
                        enrichment: bool = False,
                        ) -> list[int] | np.ndarray:
        total_counts = np.zeros(2*window_size+1)
        total_overlap = 0
        for _, row in intersect_df.iterrows():
            compartment_strand = row['strand']
            if compartment_strand == "?":
                continue
            start = int(row['start'])
            end = int(row['end'])
            # 10 -- 21 -- 31 ] 32 = 2 * w + 10 + 2
            # if end != window_size * 2 + 2 + start:
            #    print('woops?')
            #    end = start + window_size * 2 + 2
            
            motif_start = int(row['motif_start'])
            motif_end = int(row['motif_end'])

            temp_counts = np.zeros(2*window_size+1)
            overlap = int(row['overlap'])
            total_overlap += overlap
            origin = end - window_size - 1
            L = max(0, window_size - (origin - motif_start))
            U = min(2 * window_size + 1, window_size - (origin - motif_end))

            assert L <= U
            overlap_start = max(motif_start, start)
            overlap_end = min(motif_end, end)
            overlap_length = overlap_end - overlap_start
            assert overlap == overlap_length, f"{row}-{U}-{L}-{origin}"
            assert overlap == U-L, f"{row}-{U}-{L}-{origin}"
            temp_counts[L:U] += 1

            if compartment_strand == "-":
                temp_counts = temp_counts[::-1]
            total_counts += temp_counts
        total_sum = int(np.sum(total_counts))
        assert total_overlap == total_sum, f"Overlap: {total_overlap} vs. Calculated overlap {total_sum}."

        if enrichment:
            total_counts = total_counts / np.mean(total_counts)
            name = "Enrichment"
        else:
            name = "Occurrences"

        if return_frame:
            total_counts = pd.Series(total_counts)\
                                .to_frame(name=name)
            total_counts.index = total_counts.index - window_size 
            total_counts = total_counts.reset_index()
        elif not return_array:
            total_counts = list(total_counts)
        return total_counts

    def extract_PWM(self, 
                    intersect_df: pd.DataFrame, 
                    window_size: int, 
                    return_frame: bool = True
                    ) -> dict[str, list[int]]:
        total_counts = {n: [0 for _ in range(2*window_size+1)] for n in self.nucleotides}
        total_overlap = 0
        for _, row in intersect_df.iterrows():
                compartment_strand = row['strand']
                if compartment_strand == "?":
                    continue
                start = int(row['start'])
                end = int(row['end'])
                motif_start = int(row['motif_start'])
                motif_end = int(row['motif_end'])
                sequence = row['sequence']
                if compartment_strand == "-":
                    sequence = "".join(PWMExtractor.invert(n) for n in sequence)[::-1]
                overlap = int(row['overlap'])
                total_overlap += overlap
                origin = end - window_size - 1
                L = max(0, window_size - (origin - motif_start))
                U = min(2 * window_size + 1, window_size - (origin - motif_end))

                assert L <= U
                overlap_start = max(motif_start, start)
                overlap_end = min(motif_end, end)
                overlap_length = overlap_end - overlap_start
                assert overlap == overlap_length

                overlapping_sequence = sequence[max(0, start-motif_start): min(end-motif_start, len(sequence))]
                assert len(overlapping_sequence) == overlap == U-L

                for idx, pos in enumerate(range(L, U)):
                    if compartment_strand == "-":
                        index = 2 * window_size - pos
                    else:
                        index = pos
                    nucl = overlapping_sequence[idx]
                    total_counts[nucl][index] += 1
        assert total_overlap == sum(sum(v) for v in total_counts.values())
        return total_counts
