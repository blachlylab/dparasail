module dparasail.alignment;
import dparasail.result;
import parasail;
import std.utf;


/*
* Acts as a profile of settings for reusing 
* alignment settings
*/
struct Parasail
{
    const(parasail_matrix_t)* score_matrix = null;
    int gap;
    int open;
    parasail_profile_t* profile;
    string s1;

    @disable this();
    this(string matrix, int open, int gap)
    {
        this.score_matrix = parasail_matrix_lookup(toUTFz!(char*)(matrix));
        this.open = open;
        this.gap = gap;
    }

    this(string alphabet, int match, int mismatch, int open, int gap)
    {
        this.score_matrix = parasail_matrix_create(toUTFz!(char*)(alphabet), match, mismatch);
        this.open = open;
        this.gap = gap;
    }

    this(string matrix, string s1, int open, int gap)
    {
        this.score_matrix = parasail_matrix_lookup(toUTFz!(char*)(matrix));
        this.profile = parasail_profile_create_sat(toUTFz!(char*)(s1),
                cast(int) s1.length, this.score_matrix);
        this.s1 = s1;
        this.open = open;
        this.gap = gap;
    }

    this(string alphabet, string s1, int match, int mismatch, int open, int gap)
    {
        this.score_matrix = parasail_matrix_create(toUTFz!(char*)(alphabet), match, mismatch);
        this.profile = parasail_profile_create_sat(toUTFz!(char*)(s1),
                cast(int) s1.length, this.score_matrix);
        this.s1 = s1;
        this.open = open;
        this.gap = gap;
    }

    ~this()
    {
        this.close;
    }

    ParasailResult sw_striped(string s1, string s2)
    {
        return sw_striped(toUTFz!(char*)(s1), toUTFz!(char*)(s2),
                cast(int) s1.length, cast(int) s2.length);
    }

    ParasailResult sw_striped(char* s1, char* s2, int s1Len, int s2Len)
    {
        ParasailResult p;
        p.seq1 = s1;
        p.seq2 = s2;
        p.seq1Len = s1Len;
        p.seq2Len = s2Len;
        p.result = sw_striped(p);
        p.alnSettings = &this;
        p.cigar = p.get_cigar;
        return p;
    }

    parasail_result_t* sw_striped(ParasailResult p)
    {
        return parasail_sw_trace_striped_16(p.seq1, p.seq1Len, p.seq2,
                p.seq2Len, open, gap, score_matrix);
    }

    ParasailResult nw_scan(string s1, string s2)
    {
        return nw_scan(toUTFz!(char*)(s1), toUTFz!(char*)(s2),
                cast(int) s1.length, cast(int) s2.length);
    }

    ParasailResult nw_scan(char* s1, char* s2, int s1Len, int s2Len)
    {
        ParasailResult p;
        p.seq1 = s1;
        p.seq2 = s2;
        p.seq1Len = s1Len;
        p.seq2Len = s2Len;
        p.result = nw_scan(p);
        p.alnSettings = &this;
        p.cigar = p.get_cigar;
        return p;
    }

    parasail_result_t* nw_scan(ParasailResult p)
    {
        return parasail_nw_trace_scan_16(p.seq1, p.seq1Len, p.seq2, p.seq2Len,
                open, gap, score_matrix);
    }

    ParasailResult aligner(string alg, string output_option = "trace",
            string impl_option = "striped", string sol_width = "sat")(string s1, string s2)
    {
        return aligner!(alg, output_option, impl_option, sol_width)(toUTFz!(char*)(s1),
                toUTFz!(char*)(s2), cast(int) s1.length, cast(int) s2.length);
    }

    ParasailResult aligner(string alg, string output_option, string impl_option, string sol_width)(
            char* s1, char* s2, int s1Len, int s2Len)
    {
        ParasailResult p;
        p.seq1 = s1;
        p.seq2 = s2;
        p.seq1Len = s1Len;
        p.seq2Len = s2Len;
        p.result = aligner!(alg, output_option, impl_option, sol_width)(p);
        p.alnSettings = &this;
        p.cigar = p.get_cigar;
        return p;
    }

    parasail_result_t* aligner(string alg, string output_option,
            string impl_option, string sol_width)(ParasailResult p)
    {
        mixin("return parasail_" ~ alg ~ "_" ~ output_option ~ "_" ~ impl_option ~ "_"
                ~ sol_width ~ "(p.seq1,p.seq1Len,p.seq2,p.seq2Len,open,gap,score_matrix);");
    }

    ParasailResult databaseAligner(string alg, string impl_option = "striped")(string s2)
    {
        assert(!(this.profile is null));
        return databaseAligner!(alg, impl_option)(toUTFz!(char*)(s2), cast(int) s2.length);
    }

    ParasailResult databaseAligner(string alg, string impl_option)(char* s2, int s2Len)
    {
        ParasailResult p;
        p.seq2 = s2;
        p.seq2Len = s2Len;
        p.result = databaseAligner!(alg, impl_option)(p);
        p.alnSettings = &this;
        return p;
    }

    parasail_result_t* databaseAligner(string alg, string impl_option)(ParasailResult p)
    {
        mixin(
                "return parasail_" ~ alg ~ "_" ~ impl_option
                ~ "_profile_sat(this.profile,p.seq2,p.seq2Len,open,gap);");
    }

    void close()
    {
        parasail_matrix_free(cast(parasail_matrix*) score_matrix);
        if (!(this.profile is null))
        {
            parasail_profile_free(this.profile);
        }
    }
}