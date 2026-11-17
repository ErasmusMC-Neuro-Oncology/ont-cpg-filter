from array import array

import pysam
import pytest

from ont_cpg_filter.filter import Stats, filter_read, parse_mm
from ont_cpg_filter.reference import Reference


@pytest.fixture
def fasta(tmp_path):
    """chr1 = ACGTCAGA: CpG at 1-2 and 4-... only C1 and G2 form a CpG; C4 does not."""
    path = tmp_path / "ref.fa"
    path.write_text(">chr1\nACGTCAGA\n")
    pysam.faidx(str(path))
    return str(path)


def test_forward_and_reverse_cpg(fasta):
    with Reference(fasta, motif="CG") as ref:
        # 0-based: A0 C1 G2 T3 C4 A5 G6 A7
        assert ref.keeps("chr1", 1, "C")  # C of CpG on + strand
        assert ref.keeps("chr1", 2, "G")  # G of CpG (reverse-strand call)
        assert not ref.keeps("chr1", 4, "C")  # C4 is followed by A
        assert not ref.keeps("chr1", 0, "C")
        assert not ref.keeps("chr1", 1, "G")  # C1 is not preceded by C
        assert not ref.keeps("chr1", 7, "A")


def test_invalid_motif_raises(fasta):
    with pytest.raises(ValueError):
        Reference(fasta, motif="CGG")


def test_parse_mm_roundtrip():
    entries = parse_mm("C+mh?,0,3;G-m,2;")
    assert [e.seq_base for e in entries] == ["C", "C"]
    assert entries[0].codes == "mh" and entries[0].flag == "?"
    assert entries[0].skips == [0, 3]
    assert entries[1].seq_base == "C" and entries[1].skips == [2]
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
    with Reference(fasta, motif="CG") as ref:
        filter_read(read, ref, stats)

    assert read.get_tag("MM") == "C+m?,0;"
    assert list(read.get_tag("ML")) == [200]
    assert stats.kept == 1 and stats.removed == 1
