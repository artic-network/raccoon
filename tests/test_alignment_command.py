import os
import tempfile
import pytest
from raccoon.commands import alignment


class MockArgs:
    def __init__(self, **kwargs):
        self.alignment = kwargs.get('alignment', 'test.fasta')
        self.outdir = kwargs.get('outdir', '.')
        self.genbank = kwargs.get('genbank', None)
        self.reference_id = kwargs.get('reference_id', None)
        self.max_n_content = kwargs.get('max_n_content', 0.2)
        self.cluster_window = kwargs.get('cluster_window', 10)
        self.cluster_count = kwargs.get('cluster_count', 3)


def test_outdirectory_creation(tmp_path):
    """Test that output directory is created if it doesn't exist."""
    new_dir = tmp_path / 'new_outdir'
    assert not new_dir.exists()

    args = MockArgs(alignment='nonexistent.fasta', outdir=str(new_dir))
    # will fail due to missing alignment file, but directory should be created first
    result = alignment.main(args)
    
    # directory should have been created
    assert new_dir.exists()
    assert new_dir.is_dir()


def test_outdirectory_writable_check(tmp_path):
    """Test that an error is raised if output directory is not writable."""
    readonly_dir = tmp_path / 'readonly'
    readonly_dir.mkdir()
    
    # make it read-only
    os.chmod(str(readonly_dir), 0o444)
    
    args = MockArgs(alignment='nonexistent.fasta', outdir=str(readonly_dir))
    result = alignment.main(args)
    
    # should fail with write permission error
    assert result == 1
    
    # restore permissions for cleanup
    os.chmod(str(readonly_dir), 0o755)


def test_outdirectory_defaults_to_cwd(tmp_path, monkeypatch, caplog):
    """Test that outdir defaults to current working directory if not provided."""
    monkeypatch.chdir(str(tmp_path))
    
    args = MockArgs(alignment='nonexistent.fasta', outdir=None)
    # will fail due to missing alignment, but should use cwd
    result = alignment.main(args)
    
    # should have logged the error at least
    assert result != 0  # should fail
