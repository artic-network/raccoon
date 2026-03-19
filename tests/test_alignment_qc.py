import textwrap
import os

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from raccoon.utils import alignment_functions as af


def test_find_high_n_sequences(tmp_path):
    r1 = SeqRecord(Seq('ATGATGATGATG'), id='seq1')
    r2 = SeqRecord(Seq('ATGNNNNNNATG'), id='seq2')
    aln = MultipleSeqAlignment([r1, r2])
    path = tmp_path / 'test.fasta'
    AlignIO.write(aln, str(path), 'fasta')

    aln2 = AlignIO.read(str(path), 'fasta')
    flagged = af.find_high_N_sequences(aln2, threshold=0.3)
    assert len(flagged) == 1
    assert flagged[0][0] == 'seq2'


def test_analyze_alignment_unique_snps_and_flags():
    # alignment where s2 has three unique SNPs, one near N and one near a gap
    ref_seq = 'AAAAAAAAAAAA'
    ref = SeqRecord(Seq(ref_seq), id='ref')
    s1 = SeqRecord(Seq(ref_seq), id='s1')
    s2 = SeqRecord(Seq('AATAGNA-CAAA'), id='s2')

    aln = MultipleSeqAlignment([ref, s1, s2])
    unique_mutations, snps_near_n, snps_near_gap, clustered_snps = af.analyze_alignment(
        aln, n_window=2, gap_window=1, snp_window=10, snp_count=3
    )

    assert unique_mutations['s2'] == {2, 4, 8}
    assert snps_near_n['s2'] == {4}
    assert snps_near_gap['s2'] == {8}
    assert clustered_snps['s2'] == {2, 4, 8}


def test_sliding_window():
    assert list(af.sliding_window([1, 2, 3, 4], 2)) == [[1, 2], [2, 3], [3, 4]]
    assert list(af.sliding_window([1, 2], 3)) == [[1, 2]]


def test_find_frame_breaking_indels(tmp_path):
    # reference sequence with a CDS covering the full length
    ref_record = SeqRecord(Seq('ATGAAATTT'), id='ref', name='ref')
    ref_record.annotations['molecule_type'] = 'DNA'
    ref_record.features.append(SeqFeature(FeatureLocation(0, 9), type='CDS'))

    genbank_path = tmp_path / 'ref.gb'
    with open(genbank_path, 'w') as handle:
        from Bio import SeqIO
        SeqIO.write(ref_record, handle, 'genbank')

    # alignment where seq2 has a single gap in the CDS region -> frame break
    ref_aln = SeqRecord(Seq('ATGAAATTT'), id='ref')
    seq2 = SeqRecord(Seq('ATG-AATTT'), id='seq2')
    aln = MultipleSeqAlignment([ref_aln, seq2])

    sites = af.find_frame_breaking_indels(aln, str(genbank_path), reference_id='ref')
    assert sites, "Expected frame-breaking sites to be reported"
    # ensure at least one site includes seq2 in present_in
    assert any('seq2' in v['present_in'] for v in sites.values())


def test_run_alignment_qc_with_reference_fasta_matching_length(tmp_path):
    """Test that run_alignment_qc validates reference FASTA with matching length."""
    # Create alignment (12 bp)
    ref = SeqRecord(Seq('ATGATGATGATG'), id='seq1')
    seq2 = SeqRecord(Seq('ATGATGATGATG'), id='seq2')
    aln = MultipleSeqAlignment([ref, seq2])
    
    aln_path = tmp_path / 'test.fasta'
    AlignIO.write(aln, str(aln_path), 'fasta')
    
    # Create reference FASTA with same length (12 bp)
    ref_rec = SeqRecord(Seq('ATGATGATGATG'), id='reference')
    ref_path = tmp_path / 'reference.fasta'
    with open(ref_path, 'w') as handle:
        from Bio import SeqIO
        SeqIO.write(ref_rec, handle, 'fasta')
    
    # Run with reference FASTA
    outdir = tmp_path / 'output'
    outdir.mkdir()
    
    result = af.run_alignment_qc(
        str(aln_path),
        outdir=str(outdir),
        reference_fasta_path=str(ref_path)
    )
    
    assert result, "Expected result dict to be non-empty"
    assert 'high_n_sequences' in result
    assert 'mask_sites.csv' in result.get('mask_file', '')


def test_run_alignment_qc_with_reference_fasta_mismatched_length(tmp_path):
    """Test that run_alignment_qc warns when reference FASTA has different length."""
    # Create alignment (12 bp)
    ref = SeqRecord(Seq('ATGATGATGATG'), id='seq1')
    seq2 = SeqRecord(Seq('ATGATGATGATG'), id='seq2')
    aln = MultipleSeqAlignment([ref, seq2])
    
    aln_path = tmp_path / 'test.fasta'
    AlignIO.write(aln, str(aln_path), 'fasta')
    
    # Create reference FASTA with different length (10 bp instead of 12)
    ref_rec = SeqRecord(Seq('ATGATGATAG'), id='reference')
    ref_path = tmp_path / 'reference.fasta'
    with open(ref_path, 'w') as handle:
        from Bio import SeqIO
        SeqIO.write(ref_rec, handle, 'fasta')
    
    # Run with reference FASTA - should still work but warn about length mismatch
    outdir = tmp_path / 'output'
    outdir.mkdir()
    
    result = af.run_alignment_qc(
        str(aln_path),
        outdir=str(outdir),
        reference_fasta_path=str(ref_path)
    )
    
    # Should still return results even with length mismatch
    assert result, "Expected result dict to be non-empty even with length mismatch"
    assert 'high_n_sequences' in result


def test_aln_qc_command_rejects_both_genbank_and_reference_fasta(tmp_path):
    """Test that aln-qc command rejects both --genbank and --reference-fasta."""
    from raccoon.commands import alignment as aln_cmd
    
    # Create test FASTA alignment
    aln = MultipleSeqAlignment([
        SeqRecord(Seq('ATGATGATGATG'), id='seq1'),
        SeqRecord(Seq('ATGATGATGATG'), id='seq2')
    ])
    aln_path = tmp_path / 'test.fasta'
    AlignIO.write(aln, str(aln_path), 'fasta')
    
    # Create dummy genbank and reference fasta
    ref_rec = SeqRecord(Seq('ATGATGATGATG'), id='reference')
    ref_fasta = tmp_path / 'reference.fasta'
    with open(ref_fasta, 'w') as handle:
        from Bio import SeqIO
        SeqIO.write(ref_rec, handle, 'fasta')
    
    ref_gb = SeqRecord(Seq('ATGATGATGATG'), id='ref', name='ref')
    ref_gb.annotations['molecule_type'] = 'DNA'
    ref_gb.features.append(SeqFeature(FeatureLocation(0, 12), type='CDS'))
    
    genbank_path = tmp_path / 'ref.gb'
    with open(genbank_path, 'w') as handle:
        from Bio import SeqIO
        SeqIO.write(ref_gb, handle, 'genbank')
    
    # Create mock args with both genbank and reference_fasta
    class MockArgs:
        def __init__(self):
            self.alignment = str(aln_path)
            self.outdir = str(tmp_path / 'output')
            self.genbank = str(genbank_path)
            self.reference_fasta = str(ref_fasta)
            self.reference_id = 'ref'
            self.sequence_type = 'nt'
            self.max_n_content = 0.5
            self.cluster_window = 10
            self.cluster_count = 3
            self.no_flag_clustered = False
            self.no_flag_n_adjacent = False
            self.no_flag_gap_adjacent = False
            self.no_flag_frame_break = False
            self.flag_removal_threshold = 20
    
    args = MockArgs()
    result = aln_cmd.main(args)
    
    # Should return error code 1 due to conflicting arguments
    assert result == 1, "Expected error when both --genbank and --reference-fasta are provided"
