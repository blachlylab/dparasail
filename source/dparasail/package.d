module dparasail;

public import dparasail.alignment;
public import dparasail.result;


unittest
{
    import std.stdio;
    import std.math : approxEqual;

    auto p = Parasail("dnafull", 3, 2);
    auto q = p.sw_striped("GATTA", "GACTA");

    assert(q.queryPosition == 0);
    assert(q.position == 0);
    assert(q.score == 16);
    assert(q.cigar.alignedLength == 5);
    assert(q.cigar.toString == "2=1X2=");
    q.writeAlignment;
    assert(approxEqual(q.identity, 0.8));
    assert(approxEqual(q.similarity, 1.0));
    assert(approxEqual(q.gapRate, 0.0));
    

    auto q2 = p.sw_striped("GGCTTCTGATCAGGCTTCT", "GACTA");

    assert(q2.queryPosition == 7);
    assert(q2.position == 0);
    assert(q2.score == 14);
    assert(q2.cigar.alignedLength == 5);
    assert(q2.cigar.toString == "7S2=1D1=1I1=7S");
    q2.writeAlignment;
    assert(approxEqual(q2.identity, 0.667));
    assert(approxEqual(q2.similarity, 0.667));
    assert(approxEqual(q2.gapRate, 0.333));

    auto q3 = p.sw_striped("GACTAA", "GGCTTCTGATCAGGCTTCT");

    assert(q3.queryPosition == 0);
    assert(q3.position == 7);
    assert(q3.score == 14);
    assert(q3.cigar.alignedLength == 5);
    assert(q3.cigar.toString == "2=1D1=1I1=1S");
    q3.writeAlignment;
    assert(approxEqual(q3.identity, 0.667));
    assert(approxEqual(q3.similarity, 0.667));
    assert(approxEqual(q3.gapRate, 0.333));

    auto q4 = p.nw_scan("GATTA", "GACTA");

    assert(q4.queryPosition == 0);
    assert(q4.position == 0);
    assert(q4.score == 16);
    assert(q4.cigar.alignedLength == 5);
    assert(q4.cigar.toString == "2=1X2=");
    q4.writeAlignment;
    assert(approxEqual(q4.identity, 0.8));
    assert(approxEqual(q4.similarity, 1.0));
    assert(approxEqual(q4.gapRate, 0.0));

    auto q5 = p.nw_scan("GGCTTCTGATCAGGCTTCT", "GACTA");

    assert(q5.queryPosition == 0);
    assert(q5.position == 0);
    assert(q5.score == -14);
    assert(q5.cigar.alignedLength == 5);
    assert(q5.cigar.toString == "1=1X2=4I1=10S");
    q5.writeAlignment;
    writeln(q5.identity);
    writeln(q5.similarity);
    writeln(q5.gapRate);
    assert(approxEqual(q5.identity, 0.444));
    assert(approxEqual(q5.similarity, 0.556));
    assert(approxEqual(q5.gapRate, 0.444));

    // TODO: reconcile issues where cigar manipulation for 
    // use with bams changes output statistics

    // auto q52 = p.aligner!("nw","stats")("GGCTTCTGATCAGGCTTCT", "GACTA");

    // writeln(q52.queryPosition);
    // writeln(q52.position);
    // writeln(q52.score);
    // assert(q52.queryPosition == 0);
    // assert(q52.position == 0);
    // assert(q52.score == -14);
    // q5.writeAlignment;
    // writeln(q52.identity);
    // writeln(q52.similarity);
    // writeln(q52.gapRate);
    // assert(approxEqual(q5.identity, 0.444));
    // assert(approxEqual(q5.similarity, 0.556));
    // assert(approxEqual(q5.gapRate, 0.444));
    // assert(approxEqual(q5.identity, 0.211));
    // assert(approxEqual(q5.similarity, 0.263));
    // assert(approxEqual(q5.gapRate, 0.737));

    auto q6 = p.nw_scan("GACTAA", "GGCTTCTGATCAGGCTTCT");
    assert(q6.queryPosition == 0);
    assert(q6.position == 0);
    assert(q6.score == -8);
    assert(q6.cigar.alignedLength == 12);
    assert(q6.cigar.toString == "1=1X2=4D1=2D1=");
    q6.writeAlignment;
    writeln(q6.identity);
    writeln(q6.similarity);
    writeln(q6.gapRate);
    assert(approxEqual(q6.identity, 0.4167));
    assert(approxEqual(q6.similarity, 0.5));
    assert(approxEqual(q6.gapRate, 0.5));

    // TODO: report potiential bug as 
    // parasail similarity == identity 
    // when it shoudl not be

    // auto q62 = p.aligner!("nw","stats")("GACTAA", "GGCTTCTGATCAGGCTTCT");
    // writeln(q62.queryPosition);
    // writeln(q62.position);
    // writeln(q62.score);
    // writeln(q62.identity);
    // writeln(q62.similarity);
    // writeln(q62.gapRate);
    // assert(q6.queryPosition == 0);
    // assert(q6.position == 0);
    // assert(q6.score == -8);
    // assert(q6.cigar.alignedLength == 12);
    // assert(q6.cigar.toString == "1=1X2=4D1=2D1=");
    // q6.writeAlignment;
    // writeln(q6.identity);
    // writeln(q6.similarity);
    // writeln(q6.gapRate);
    // assert(approxEqual(q6.identity, 0.4167));
    // assert(approxEqual(q6.similarity, 0.5));
    // assert(approxEqual(q6.gapRate, 0.5));
    // assert(approxEqual(q6.identity, 0.263));
    // assert(approxEqual(q6.similarity, 0.316));
    // assert(approxEqual(q6.gapRate, 0.684));
    
}
