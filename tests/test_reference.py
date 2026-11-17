from array import array

import numpy as np
import pysam
import pytest

from ont_cpg_filter.filter import (
    Stats,
    filter_read,
    parse_mm,
    query_to_ref,
    query_to_ref_at,
    ref_blocks,
)
from ont_cpg_filter.index import build_index
from ont_cpg_filter.reference import Reference


@pytest.fixture
def fasta(tmp_path):
    """chr1 = ACGTCAGA: one CpG at 1-2; C4 is followed by A, so not a CpG."""
    path = tmp_path / "ref.fa"
    path.write_text(">chr1\nACGTCAGA\n")
    pysam.faidx(str(path))
    return str(path)


def test_scalar_cpg_check(fasta):
    with Reference(fasta, motif="CG") as ref:
        # 0-based: A0 C1 G2 T3 C4 A5 G6 A7
        assert ref.keeps("chr1", 1, "C")
        assert ref.keeps("chr1", 2, "G")
        assert not ref.keeps("chr1", 4, "C")
        assert not ref.keeps("chr1", 0, "C")
        assert not ref.keeps("chr1", 1, "G")
        assert not ref.keeps("chr1", 7, "A")


def test_bitmap_site_mask(fasta):
    with Reference(fasta, motif="CG", use_index=False) as ref:
        rpos = np.array([1, 2, 4, -1, 99], dtype=np.int64)
        fwd = ref.site_mask("chr1", rpos, forward=True)
        rev = ref.site_mask("chr1", rpos, forward=False)
        assert fwd.tolist() == [True, False, False, False, False]
        assert rev.tolist() == [False, True, False, False, False]


def test_index_roundtrip_matches_built_bitmaps(fasta):
    build_index(fasta, "CG")
    rpos = np.array([1, 2, 4, -1, 99], dtype=np.int64)
    with Reference(fasta, motif="CG", use_index=True) as ref:
        fwd = ref.site_mask("chr1", rpos, forward=True)
        rev = ref.site_mask("chr1", rpos, forward=False)
    assert fwd.tolist() == [True, False, False, False, False]
    assert rev.tolist() == [False, True, False, False, False]


def test_query_to_ref_handles_clips_indels():
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 100}]})
    read = pysam.AlignedSegment(header)
    read.query_sequence = "A" * 10
    read.cigarstring = "2S3M1I2M1D2M"
    read.reference_start = 10
    # S: -1 -1 | M3: 10 11 12 | I: -1 | M2: 13 14 | D skips 15 | M2: 16 17
    assert query_to_ref(read, 10).tolist() == [-1, -1, 10, 11, 12, -1, 13, 14, 16, 17]
    blocks = ref_blocks(read)
    idx = np.array([0, 2, 5, 9, 1], dtype=np.int64)
    assert query_to_ref_at(blocks, idx).tolist() == [-1, 10, -1, 17, -1]


def test_invalid_motif_raises(fasta):
    with pytest.raises(ValueError):
        Reference(fasta, motif="CGG")


def test_parse_mm_roundtrip():
    entries = parse_mm("C+mh?,0,3;G-m,2;")
    assert entries[0].codes == "mh" and entries[0].flag == "?"
    assert entries[0].skips.tolist() == [0, 3]
    assert entries[1].base == "G" and entries[1].strand == "-"
    assert entries[1].skips.tolist() == [2]
    assert entries[0].to_mm() == "C+mh?,0,3;"


def test_filter_read_keeps_only_cpg_calls(fasta):
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 8}]})
    read = pysam.AlignedSegment(header)
    read.query_name = "r1"
    read.query_sequence = "ACGTCAGA"
    read.cigarstring = "8M"
    read.reference_id = 0
    read.reference_start = 0
    read.flag = 0
    # Candidate C's in SEQ: index 1 (CpG) and index 4 (non-CpG); two calls.
    read.set_tag("MM", "C+m,0,0;")
    read.set_tag("ML", array("B", [200, 150]))

    stats = Stats()
    with Reference(fasta, motif="CG", use_index=False) as ref:
        filter_read(read, ref, stats)

    assert read.get_tag("MM") == "C+m?,0;"
    assert list(read.get_tag("ML")) == [200]
    assert stats.kept == 1 and stats.removed == 1
