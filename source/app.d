import std.stdio;
import std.array:array;
import dparasail;
void main()
{
	auto p=Parasail("dnafull",3,2);
	auto q=p.sw_striped("GATTA","GACTA");
	q=p.aligner!("sw")("GATTA","GACTA");
	writeln(q.seq1[0..q.seq1Len]);
	writeln(q.get_cigar(p.score_matrix));
	writeln(q.beg_query);
	writeln(q.beg_ref);
	//assert(q.cigar=="2=1X2=");
	q=p.sw_striped("GGCTTCTGATCAGGCTTCT","GACTA");
	writeln(q.seq1[0..q.seq1Len]);
	writeln(q.get_cigar(p.score_matrix));
	writeln(q.beg_query);
	writeln(q.beg_ref);
	//assert(q.cigar=="7S2=1D1=1I1=7S");
	q=p.sw_striped("GACTAA","GGCTTCTGATCAGGCTTCT");
	writeln(q.seq1[0..q.seq1Len]);
	writeln(q.get_cigar(p.score_matrix));
	writeln(q.beg_query);
	writeln(q.beg_ref);
	q=p.nw_scan("GATTA","GACTA");
	writeln(q.seq1[0..q.seq1Len]);
	writeln(q.get_cigar(p.score_matrix));
	writeln(q.beg_query);
	writeln(q.beg_ref);
	//assert(q.cigar=="2=1X2=");
	q=p.nw_scan("GGCTTCTGATCAGGCTTCT","GACTA");
	writeln(q.seq1[0..q.seq1Len]);
	writeln(q.get_cigar(p.score_matrix));
	writeln(q.beg_query);
	writeln(q.beg_ref);
	//assert(q.cigar=="1=1X2=4I1=10S");
	q=p.nw_scan("GACTAA","GGCTTCTGATCAGGCTTCT");
	writeln(q.seq1[0..q.seq1Len]);
	writeln(q.get_cigar(p.score_matrix));
	writeln(q.beg_query);
	writeln(q.beg_ref);
	auto cigar=parasail_result_get_cigar(q.result,q.seq1,q.seq1Len,q.seq2,q.seq2Len,p.score_matrix);
	writeln(parasail_cigar_decode(cigar));
	p = Parasail("dnafull","GATTA",3,2);
	writeln(p.databaseAligner!"nw"("GACTA").result.score);
	writeln(p.databaseAligner!"nw"("GGCTTCTGATCAGGCTTCT").result.score);
	writeln(p.databaseAligner!"nw"("GGCTTCTGATCAGGCTTCTGGCTTCTGATCAGGCTTCT").result.score);
	p.close;
}
