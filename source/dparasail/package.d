module dparasail;

public import dparasail.alignment;
public import dparasail.result;


unittest
{
    import std.stdio;

    auto p = Parasail("dnafull", 3, 2);
    auto q = p.sw_striped("GATTA", "GACTA");

    assert(q.beg_query == 0);
    assert(q.beg_ref == 0);
    assert(q.result.score == 16);
    assert(q.cigar.alignedLength == 5);
    assert(q.cigar.toString == "2=1X2=");
    q.writeAlignment;

    q = p.sw_striped("GGCTTCTGATCAGGCTTCT", "GACTA");

    assert(q.beg_query == 7);
    assert(q.beg_ref == 0);
    assert(q.result.score == 14);
    assert(q.cigar.alignedLength == 5);
    assert(q.cigar.toString == "7S2=1D1=1I1=7S");
    q.writeAlignment;

    q = p.sw_striped("GACTAA", "GGCTTCTGATCAGGCTTCT");

    assert(q.beg_query == 0);
    assert(q.beg_ref == 7);
    assert(q.result.score == 14);
    assert(q.cigar.alignedLength == 5);
    assert(q.cigar.toString == "2=1D1=1I1=1S");
    q.writeAlignment;

    q = p.nw_scan("GATTA", "GACTA");

    assert(q.beg_query == 0);
    assert(q.beg_ref == 0);
    assert(q.result.score == 16);
    assert(q.cigar.alignedLength == 5);
    assert(q.cigar.toString == "2=1X2=");
    q.writeAlignment;

    q = p.nw_scan("GGCTTCTGATCAGGCTTCT", "GACTA");

    assert(q.beg_query == 0);
    assert(q.beg_ref == 0);
    assert(q.result.score == -14);
    assert(q.cigar.alignedLength == 5);
    assert(q.cigar.toString == "1=1X2=4I1=10S");
    q.writeAlignment;

    q = p.nw_scan("GACTAA", "GGCTTCTGATCAGGCTTCT");
    assert(q.beg_query == 0);
    assert(q.beg_ref == 0);
    assert(q.result.score == -8);
    assert(q.cigar.alignedLength == 12);
    assert(q.cigar.toString == "1=1X2=4D1=2D1=");
    q.writeAlignment;
    // auto cigar=parasail_result_get_cigar(q.result,q.seq1,q.seq1Len,q.seq2,q.seq2Len,q.score_matrix);
    // writeln(parasail_cigar_decode(cigar));
}
