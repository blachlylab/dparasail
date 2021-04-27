module dparasail.result;

import std.stdio;
import std.conv;
import std.utf;

// import bio.std.hts.bam.cigar;
import dhtslib.sam.cigar;
import std.algorithm : map, filter;
import std.algorithm.iteration : sum;
import dparasail.alignment;
import parasail;

/*
* Stores alignment results (scores, cigar, alignment starts).
* Contains pointers to underlying parasail data. Owns parasail_result_t.
*/
struct ParasailResult
{
    Parasail* alnSettings;
    char* seq1;
    char* seq2;
    int seq1Len;
    int seq2Len;
    //0-based
    int beg_query;
    int beg_ref;
    Cigar cigar;
    parasail_result_t* result;
    Cigar get_cigar()
    {
        parasail_cigar_t* cigar;
        Cigar cigar_string;
        cigar = parasail_result_get_cigar(result, seq1, seq1Len, seq2, seq2Len,
                alnSettings.score_matrix);
        beg_query = cigar.beg_query;
        beg_ref = cigar.beg_ref;
        cigar_string = Cigar(((cast(CigarOp*) cigar.seq)[0 .. cigar.len]).dup);
        parasail_cigar_free(cigar);
        //if *
        if (cigar_string.length == 0)
        {
            return cigar_string;
        }
        //if 30I8M make 30S8M
        if (cigar_string[0].op == Ops.INS)
        {
            //move beg_query as well as it is accounted for
            cigar_string[0] = CigarOp(cigar_string[0].length + this.beg_query, Ops.SOFT_CLIP);
            this.beg_query += cigar_string[0].length;
        }
        //else if 30D8M make 8M and move ref start
        else if (cigar_string[0].op == Ops.DEL)
        {
            this.beg_ref = this.beg_ref + cigar_string[0].length;
            cigar_string = Cigar(cigar_string[1 .. $]);
            //cigar_string[0]=CigarOperation(cigar_string[0].length+this.beg_query,'S');
        }
        //else if begin query not 0 add softclip
        else if (this.beg_query != 0)
        {
            assert(cigar_string[0].op != Ops.SOFT_CLIP);
            cigar_string = Cigar(CigarOp(this.beg_query, Ops.SOFT_CLIP) ~ cigar_string[]);
            this.beg_query = cigar_string[0].length;
        }
        ///////////////////////////////////////////////////////////
        int q_bases_covered = cast(int) cigar_string[].filter!(x => x.is_query_consuming())
            .map!(x => x.length)
            .sum;
        if (cigar_string[$ - 1].op == Ops.INS)
        {
            cigar_string[$ - 1] = CigarOp(cigar_string[$ - 1].length + this.seq1Len - q_bases_covered,
                    Ops.SOFT_CLIP);
            q_bases_covered = cigar_string[].filter!(x => x.is_query_consuming())
                .map!(x => x.length)
                .sum;
        }
        else if (cigar_string[$ - 1].op == Ops.DEL)
        {
            cigar_string = Cigar(cigar_string[0 .. ($ - 1)]);
        }
        else if (q_bases_covered != this.seq1Len)
        {
            cigar_string = Cigar(cigar_string[] ~ CigarOp(this.seq1Len - q_bases_covered,
                    Ops.SOFT_CLIP));
            q_bases_covered = cigar_string[].filter!(x => x.is_query_consuming())
                .map!(x => x.length)
                .sum;
        }
        assert(q_bases_covered == seq1Len);
        // parasail_cigar_free(cigar);  
        // Bug caused bc Cigar ctor takes reference not copy
        // Cigar implicit dtor destroys array ptr which cleans up our problem
        return cigar_string;
    }

    void writeAlignment()
    {
        parasail_traceback_generic(seq1, seq1Len, seq2, seq2Len, toUTFz!(char*)("Query:"),
                toUTFz!(char*)("Target:"), alnSettings.score_matrix, result,
                '|', '*', '*', 60, 7, 1);
    }

    void close()
    {
        parasail_result_free(result);
    }
}