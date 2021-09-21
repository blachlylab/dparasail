module dparasail.result;

import std.stdio;
import std.conv;
import std.utf;

// import bio.std.hts.bam.cigar;
import dhtslib.sam.cigar;
import std.algorithm : map, filter;
import std.algorithm.iteration : sum;
import std.array : array;
import dparasail.alignment;
import dparasail.memory;
import parasail;

/**
* Stores alignment results (scores, cigar, alignment starts).
* Contains pointers to underlying parasail data. Owns parasail_result_t.
*/
struct ParasailResult
{
    private Parasail* alnSettings;
    private string s1;
    private string s2;
    private const char* seq1;
    private const char* seq2;
    private const int seq1Len;
    private const int seq2Len;
    //0-based
    private int queryPos;
    private int refPos;
    private Cigar p_cigar;
    private bool hasCigar;
    private ParasailResultPtr result;

    private int num_gaps;
    private int num_mis;
    private int num_matches;
    private int total_length;

    /// need to use a constructor
    @disable this();

    /// ctor for strings
    this(Parasail * alnSettings, parasail_result_t * result, string s1, string s2)
    {
        this.result = ParasailResultPtr(result);
        this.alnSettings = alnSettings;
        this.s1 = s1;
        this.s2 = s2;
        auto seq1 = toUTFz!(const char*)(s1);
        auto seq2 = toUTFz!(const char*)(s2);
        auto seq1Len = cast(const int) s1.length;
        auto seq2Len = cast(const int) s2.length;
        this.seq1 = seq1;
        this.seq1Len = seq1Len;
        this.seq2 = seq2;
        this.seq2Len = seq2Len;
        assert(alnSettings.scoreMatrix);
        if(parasail_result_is_trace(this.result))
        {
            this.hasCigar = true; 
            this.p_cigar = this.get_cigar();  
        }
    }

    /// starting position of the alignment 
    /// on the reference/target sequence
    int position()
    {
        return this.refPos;
    }

    /// starting position of the alignment 
    /// on the query sequence
    int queryPosition()
    {
        return this.queryPos;
    }

    /// starting position of the alignment 
    /// on the reference/target sequence
    int endQuery()
    {
        return this.result.end_query;
    }

    /// starting position of the alignment 
    /// on the query sequence
    int endRef()
    {
        return this.result.end_ref;
    }

    /// alignment score
    int score()
    {
        return this.result.score;    
    }

    /// get identity
    float identity()
    {
        if(parasail_result_is_trace(this.result)){
            if(total_length == 0)
                cigarStats;
            return float(num_matches) / float(total_length);
        }
        else if(parasail_result_is_stats(this.result))
            return float(parasail_result_get_matches(this.result)) / float(parasail_result_get_length(this.result));
        throw new Exception("Alignment type is not trace or stats");
    }
    
    /// get similarity
    float similarity()
    {
        if(parasail_result_is_trace(this.result)){
            if(total_length == 0)
                cigarStats;
            return float(num_matches + num_mis) / float(total_length);
        }
        else if(parasail_result_is_stats(this.result))
            return float(parasail_result_get_similar(this.result)) / float(parasail_result_get_length(this.result));
        throw new Exception("Alignment type is not trace or stats");
    }

    /// get identity
    float gapRate()
    {
        if(parasail_result_is_trace(this.result)){
            if(total_length == 0)
                cigarStats;
            return float(num_gaps) / float(total_length);
        }
        else if(parasail_result_is_stats(this.result))
            return 1.0 - (float(parasail_result_get_similar(this.result)) / float(parasail_result_get_length(this.result)));
        throw new Exception("Alignment type is not trace or stats");
    }

    private void cigarStats(){
        assert(total_length == 0);
        assert(num_matches == 0);
        assert(num_mis == 0);
        assert(num_gaps == 0);
        auto coords = AlignedCoordinatesItr(cigar).filter!(x => x.cigar_op != Ops.SOFT_CLIP).array;
        foreach (i, op; coords)
        {
            switch(op.cigar_op){
                // should no longer happen
                // case Ops.SOFT_CLIP:
                //     continue;
                case Ops.DIFF:
                    num_mis++;
                    break;
                case Ops.DEL:
                case Ops.INS:
                    num_gaps++;
                    break;
                default:
                    num_matches++; 
                    break;
            }
            total_length++;
        }
    }

    /// get cigar string of alignment
    Cigar cigar(){
        if(hasCigar)    
            return this.p_cigar;
        else
            throw new Exception("This alignment type cannot produce cigar.");
    }

    /// get cigar from result and parse it
    /// default cigar from parasail is not compatible with 
    /// a bam per say. It needs to be padded with soft-clips,
    /// and have leading and lagging INDELs removed.
    /// Alignment positions also have to modified.
    private Cigar get_cigar()
    {
        // Get cigar from result
        parasail_cigar_t* cigar = parasail_result_get_cigar(this.result, seq1, seq1Len, seq2, seq2Len,
                alnSettings.scoreMatrix);
        this.queryPos = cigar.beg_query;
        this.refPos = cigar.beg_ref;
        Cigar cigar_string = Cigar(((cast(CigarOp*) cigar.seq)[0 .. cigar.len]).dup);
        parasail_cigar_free(cigar);
        
        //if *
        if (cigar_string.length == 0)
        {
            return cigar_string;
        }
        //if 30I8M make 30S8M
        if (cigar_string[0].op == Ops.INS)
        {
            //move queryPos as well as it is accounted for
            cigar_string[0] = CigarOp(cigar_string[0].length + this.queryPos, Ops.SOFT_CLIP);
            this.queryPos += cigar_string[0].length;
        }
        //else if 30D8M make 8M and move ref start
        else if (cigar_string[0].op == Ops.DEL)
        {
            this.refPos = this.refPos + cigar_string[0].length;
            cigar_string = Cigar(cigar_string[1 .. $]);
            //cigar_string[0]=CigarOperation(cigar_string[0].length+this.queryPos,'S');
        }
        //else if begin query not 0 add softclip
        else if (this.queryPos != 0)
        {
            assert(cigar_string[0].op != Ops.SOFT_CLIP);
            cigar_string = Cigar(CigarOp(this.queryPos, Ops.SOFT_CLIP) ~ cigar_string[]);
            this.queryPos = cigar_string[0].length;
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

    /// display parasail representation of alignment
    void writeAlignment()
    {
        parasail_traceback_generic(seq1, seq1Len, seq2, seq2Len, toUTFz!(char*)("Query:"),
                toUTFz!(char*)("Target:"), alnSettings.scoreMatrix, result,
                '|', '*', '*', 60, 7, 1);
    }

    /// clean up 
    private void close()
    {
        parasail_result_free(result);
    }
}