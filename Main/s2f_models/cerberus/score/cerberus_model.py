"""Loading a Cerberus fold and predicting coverage for a reference/alternate pair.

`hound_snp` covers bulk variant scoring; this is for the handful of places that
need the predicted coverage *profile* itself — the example-locus panels, and the
credible-set scoring where the statistic is restricted to one gene.

The package is `baskerville`, the PyTorch port of Calico's baskerville. Note
that a separate TensorFlow package of the same name provides `bw_w5.py` used in
../../data; they are different codebases.
"""
import json

import numpy as np
import pandas as pd
import torch
from baskerville import dataset, dna, seqnn


class CerberusFold:
    """One trained fold, ready to predict coverage over a genomic window."""

    def __init__(self, params_path, model_path, targets_path, head=0):
        with open(params_path) as handle:
            params = json.load(handle)
        model_params = params["model"]

        targets = pd.read_csv(targets_path, sep="\t", index_col=0)

        # Collapsing strand pairs halves the track count: the two strands of an
        # RNA assay are averaged into one output. Every track index used
        # downstream refers to this collapsed ordering, not the targets file's
        # row numbers.
        if "strand_pair" in targets.columns:
            dataset.targets_prep_strand(targets)
            old_to_new = dict(zip(targets.index, np.arange(len(targets))))
            model_params["strand_pair"] = np.array(
                [old_to_new[t] for t in targets.strand_pair]
            )

        self.targets = targets
        self.head = head

        self.model = seqnn.SeqNN(model_params, output_slice=targets.index)
        self.model.restore(model_path)
        self.model.ensemble_rc = True
        self.model.ensemble_shifts = [0]
        self.model.mix_dtype = torch.bfloat16
        self.model.model.eval()

        self.seq_length = model_params["seq_length"]
        self.stride = self.model.output_stride()
        self.output_length = self.model.output_length()

    def describe(self, track_index):
        return self.targets.iloc[track_index]["description"]

    def check_track(self, track_index, expected):
        """Guard against a targets file whose ordering has changed."""
        description = self.describe(track_index)
        if expected not in description:
            raise ValueError(
                f"track {track_index} is {description!r}, expected {expected!r} - "
                "the targets file ordering does not match this track index"
            )

    def bin_centers(self, center):
        """Genomic coordinate of each output bin for a window centred on `center`."""
        out_bp = self.output_length * self.stride
        out_start = center - out_bp // 2
        return out_start + (np.arange(self.output_length) + 0.5) * self.stride

    def predict(self, sequence):
        one_hot = dna.dna_1hot(sequence)
        batch = torch.tensor(
            one_hot.T[np.newaxis], dtype=torch.float32, device=self.model.device
        )
        with torch.no_grad():
            output = self.model(batch, hi=self.head)
        return output.coverage[0].float().cpu().numpy().T  # (bins, tracks)

    def predict_alleles(self, fasta, chrom, pos, ref, alt):
        """Predict coverage for both alleles of a variant.

        `pos` is 1-based, as in a VCF. Returns (ref_coverage, alt_coverage,
        bin_centers), each coverage array shaped (bins, tracks).
        """
        center = pos - 1
        seq_start = center - self.seq_length // 2
        variant_offset = center - seq_start

        ref_seq = fasta.fetch(chrom, seq_start, seq_start + self.seq_length).upper()
        if ref_seq[variant_offset] != ref:
            print(
                f"  warning: {chrom}:{pos} genome has "
                f"{ref_seq[variant_offset]}, variant says ref={ref}"
            )
        alt_seq = ref_seq[:variant_offset] + alt + ref_seq[variant_offset + 1:]

        return self.predict(ref_seq), self.predict(alt_seq), self.bin_centers(center)


def log_ratio(ref_sum, alt_sum, pseudocount=1e-6):
    """Alternate-over-reference effect, log2, as reported throughout."""
    return np.log2(alt_sum + pseudocount) - np.log2(ref_sum + pseudocount)


def exon_intervals(gtf_path, gene_name):
    """Merged exon intervals of one gene, as (start, end) 1-based inclusive."""
    import re

    pattern = re.compile(r'gene_name "([^"]+)"')
    exons = set()
    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            match = pattern.search(fields[8])
            if match and match.group(1) == gene_name:
                exons.add((int(fields[3]), int(fields[4])))
    if not exons:
        raise ValueError(f"no exons for {gene_name} in {gtf_path}")
    return sorted(exons)


def mask_intervals(bin_centers, intervals):
    """Boolean mask over output bins whose centre falls in any interval."""
    mask = np.zeros(len(bin_centers), dtype=bool)
    for start, end in intervals:
        mask |= (bin_centers >= start) & (bin_centers <= end)
    return mask
