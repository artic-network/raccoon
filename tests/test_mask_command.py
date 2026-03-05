import csv

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from raccoon.commands import mask


def _write_alignment(tmp_path):
    records = [
        SeqRecord(Seq("ACGTA"), id="s1"),
        SeqRecord(Seq("ACGTA"), id="s2"),
    ]
    aln = MultipleSeqAlignment(records)
    aln_path = tmp_path / "alignment.fasta"
    AlignIO.write(aln, str(aln_path), "fasta")
    return aln_path


def _write_mask(tmp_path):
    mask_path = tmp_path / "mask_sites.csv"
    with mask_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["flagged", "type", "Minimum", "Maximum", "Length", "present_in", "note"], lineterminator="\n")
        writer.writeheader()
        writer.writerow({
            "flagged": "2-3",
            "type": "site",
            "Minimum": "2",
            "Maximum": "3",
            "Length": "2",
            "present_in": "s1",
            "note": "masked",
        })
    return mask_path


class _Args:
    def __init__(self, alignment_path, mask_file, outdir, **overrides):
        self.alignment = str(alignment_path)
        self.mask_file = str(mask_file)
        self.outdir = str(outdir)
        self.outfile = None
        self.sequence_type = "nt"
        for key, value in overrides.items():
            setattr(self, key, value)


def test_mask_command_writes_masked_alignment(tmp_path):
    aln_path = _write_alignment(tmp_path)
    mask_path = _write_mask(tmp_path)
    outdir = tmp_path / "out"

    args = _Args(aln_path, mask_path, outdir)
    result = mask.main(args)
    assert result == 0

    masked_path = outdir / "alignment.masked.fasta"
    assert masked_path.exists()

    masked_aln = AlignIO.read(str(masked_path), "fasta")
    for rec in masked_aln:
        assert rec.seq[1] == "N"
        assert rec.seq[2] == "N"


def test_mask_command_preserves_unmasked_sites(tmp_path):
    records = [
        SeqRecord(Seq("ACGTAC"), id="s1"),
        SeqRecord(Seq("ACGTAC"), id="s2"),
    ]
    aln = MultipleSeqAlignment(records)
    aln_path = tmp_path / "alignment.fasta"
    AlignIO.write(aln, str(aln_path), "fasta")

    mask_path = tmp_path / "mask_sites.csv"
    with mask_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["flagged", "type", "Minimum", "Maximum", "Length", "present_in", "note"], lineterminator="\n")
        writer.writeheader()
        writer.writerow({
            "flagged": "5",
            "type": "site",
            "Minimum": "",
            "Maximum": "",
            "Length": "1",
            "present_in": "s2",
            "note": "masked",
        })

    outdir = tmp_path / "out"
    args = _Args(aln_path, mask_path, outdir)
    result = mask.main(args)
    assert result == 0

    masked_path = outdir / "alignment.masked.fasta"
    masked_aln = AlignIO.read(str(masked_path), "fasta")
    for rec in masked_aln:
        assert str(rec.seq)[:4] == "ACGT"
        assert rec.seq[4] == "N"
        assert rec.seq[5] == "C"


def test_mask_command_removes_sequence_records(tmp_path):
    records = [
        SeqRecord(Seq("ACGTAC"), id="s1"),
        SeqRecord(Seq("ACGTAC"), id="s2"),
    ]
    aln = MultipleSeqAlignment(records)
    aln_path = tmp_path / "alignment.fasta"
    AlignIO.write(aln, str(aln_path), "fasta")

    mask_path = tmp_path / "mask_sites.csv"
    with mask_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["flagged", "type", "Minimum", "Maximum", "Length", "present_in", "note"], lineterminator="\n")
        writer.writeheader()
        writer.writerow({
            "flagged": "s2",
            "type": "sequence_record",
            "Minimum": "1",
            "Maximum": "6",
            "Length": "6",
            "present_in": "s2",
            "note": "2,3,4",
        })

    outdir = tmp_path / "out"
    args = _Args(aln_path, mask_path, outdir)
    result = mask.main(args)
    assert result == 0

    masked_path = outdir / "alignment.masked.fasta"
    masked_aln = AlignIO.read(str(masked_path), "fasta")
    ids = [rec.id for rec in masked_aln]
    assert ids == ["s1"]
