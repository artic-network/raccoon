import csv

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from raccoon.commands import alignment


def _write_genbank(path, record):
    record.annotations["molecule_type"] = "DNA"
    with open(path, "w") as handle:
        from Bio import SeqIO
        SeqIO.write(record, handle, "genbank")


def _build_alignment(tmp_path):
    # reference sequence and CDS covering full length
    ref_seq = "ATG" * 10
    ref_record = SeqRecord(Seq(ref_seq), id="ref", name="ref")
    ref_record.features.append(SeqFeature(FeatureLocation(0, len(ref_seq)), type="CDS"))

    # baseline identical sequence
    s1 = SeqRecord(Seq(ref_seq), id="s1")

    # clustered SNPs within 10bp window (positions 4, 6, 8 are unique)
    s_cluster_list = list(ref_seq)
    for pos, nt in [(4, "C"), (6, "G"), (8, "T")]:
        s_cluster_list[pos] = nt
    s_cluster = SeqRecord(Seq("".join(s_cluster_list)), id="s_cluster")

    # SNP next to N plus some Ns to trigger high-N sequence
    s_n_list = list(ref_seq)
    s_n_list[12] = "G"
    s_n_list[13] = "N"
    s_n_list[14] = "N"
    s_n = SeqRecord(Seq("".join(s_n_list)), id="s_n")

    # frame-breaking indel (single gap within CDS)
    s_frame_list = list(ref_seq)
    s_frame_list[5] = "-"
    s_frame = SeqRecord(Seq("".join(s_frame_list)), id="s_frame")

    aln = MultipleSeqAlignment([ref_record, s1, s_cluster, s_n, s_frame])

    aln_path = tmp_path / "test_alignment.fasta"
    AlignIO.write(aln, str(aln_path), "fasta")

    gb_path = tmp_path / "ref.gb"
    _write_genbank(gb_path, ref_record)

    return aln_path, gb_path


class _Args:
    def __init__(self, alignment_path, outdir, genbank_path, **overrides):
        self.alignment = str(alignment_path)
        self.output_dir = str(outdir)
        self.genbank = str(genbank_path)
        self.reference_id = "ref"
        self.max_n_content = 0.05
        self.cluster_window = 10
        self.cluster_count = 3
        self.mask_clustered = True
        self.mask_n_adjacent = True
        self.mask_gap_adjacent = True
        self.mask_frame_break = True
        for key, value in overrides.items():
            setattr(self, key, value)


def test_alignment_qc_end_to_end(tmp_path):
    aln_path, gb_path = _build_alignment(tmp_path)
    outdir = tmp_path / "out"

    args = _Args(aln_path, outdir, gb_path)
    result = alignment.main(args)
    assert result == 0

    summary_path = outdir / "alignment_qc_summary.txt"
    mask_path = outdir / "mask_sites.csv"
    assert summary_path.exists()
    assert mask_path.exists()

    rows = list(csv.DictReader(open(mask_path)))
    assert rows, "Expected mask_sites.csv to contain at least one row"

    site_rows = [row for row in rows if row.get("type") == "site"]
    notes_by_site = [row["note"] for row in site_rows]
    present_by_site = [row["present_in"] for row in site_rows]

    assert any("frame_break" in note for note in notes_by_site)
    assert any("s_frame" in present for present in present_by_site)

    assert any("clustered_snps" in note for note in notes_by_site)
    assert any("s_cluster" in present for present in present_by_site)


def test_alignment_qc_disable_n_adjacent(tmp_path):
    aln_path, gb_path = _build_alignment(tmp_path)
    outdir = tmp_path / "out"

    args = _Args(aln_path, outdir, gb_path, mask_n_adjacent=False)
    result = alignment.main(args)
    assert result == 0

    mask_path = outdir / "mask_sites.csv"
    rows = list(csv.DictReader(open(mask_path)))
    notes_by_site = [row["note"] for row in rows if row.get("type") == "site"]

    assert not any("N_adjacent" in note for note in notes_by_site)


