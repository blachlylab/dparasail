module dparasail.alignment;
import dparasail.result;
import parasail;
import std.utf;
import std.exception : enforce; 


/*
* Acts as a profile of settings for reusing 
* alignment settings
*/
struct Parasail
{
    /// alignment scoring matrix reference
    private const(parasail_matrix_t)* score_matrix = null;

    /// using a prebuilt parasail matrix?
    private bool prebuiltMatrix;

    /// using database query?
    private bool databaseQuery;

    /// gap open penalty (must be positive)
    private int gapOpen;
    
    /// gap extension penalty (must be positive)
    private int gapExt;

    /// parasail profile reference if doing
    /// if doing database alignment
    private parasail_profile_t* profile;

    /// query sequence if doing database alignment
    private string s1;

    private int refct = 1;

    invariant(){
        assert(this.refct >=0);
    }

    /// need to use a constructor
    @disable this();

    /// use a prebuilt parasail scoring matrix
    /// and a gap and ext penalty
    this(string matrix, int open, int ext)
    {
        auto m = parasail_matrix_lookup(toUTFz!(char*)(matrix));
        this.prebuiltMatrix = true;
        this(m, open, ext);
    }

    /// build your own parasail scoring matrix
    /// specify an alphabet (GATCN), specify a match
    /// score and mismatch penalty (must be a negative number),
    /// and a gap and ext penalty
    this(string alphabet, int match, int mismatch, int open, int ext)
    {
        enforce(match > 0, "match score must be > 0");
        enforce(mismatch < 0, "mismatch penalty/score must be < 0");
        auto m = parasail_matrix_create(toUTFz!(char*)(alphabet), match, mismatch);
        this(m, open, ext);
    }

    /// use a single query sequence to database alignment
    /// against many other sequences
    /// also
    /// use a prebuilt parasail scoring matrix
    /// and a gap and ext penalty
    this(string matrix, string s1, int open, int ext)
    {
        this.s1 = s1;
        this.databaseQuery = true;
        this.prebuiltMatrix = true;
        auto m = parasail_matrix_lookup(toUTFz!(char*)(matrix));
        this.profile = parasail_profile_create_sat(toUTFz!(char*)(s1),
                cast(int) s1.length, m);
        this(m, open, ext);
    }

    /// use a single query sequence to database alignment
    /// against many other sequences
    /// also
    /// build your own parasail scoring matrix
    /// specify an alphabet (GATCN), specify a match
    /// score and mismatch penalty (must be a negative number),
    /// and a gap and ext penalty
    this(string alphabet, string s1, int match, int mismatch, int open, int ext)
    {
        this.s1 = s1;
        this.databaseQuery = true;
        auto m = parasail_matrix_create(toUTFz!(char*)(alphabet), match, mismatch);
        this.profile = parasail_profile_create_sat(toUTFz!(char*)(s1),
                cast(int) s1.length, this.score_matrix);
        this(m, open, ext);
    }

    /// use a parasail_matrix_t directly
    /// and a gap and ext penalty
    /// not recommended, for internal use
    this(const(parasail_matrix_t) * matrix, int open, int ext)
    {
        this.score_matrix = matrix;
        enforce(open > 0, "gap open penalty must be greater than 0");
        enforce(ext > 0, "gap extension penalty must be greater than 0");
        this.gapOpen = open;
        this.gapExt = ext;
    }

    this(this){
        this.refct++;
    }

    ~this()
    {
        if(--refct == 0)
            this.close;
    }

    const(parasail_matrix_t) * scoreMatrix(){
        return this.score_matrix;
    }

    ParasailResult sw_striped(string s1, string s2)
    {
        auto seq1 = toUTFz!(const char*)(s1);
        auto seq2 = toUTFz!(const char*)(s2);
        auto seq1Len = cast(const int) s1.length;
        auto seq2Len = cast(const int) s2.length;
        auto res = parasail_sw_trace_striped_16(seq1, seq1Len, seq2,
                seq2Len, this.gapOpen, this.gapExt, this.score_matrix);
        return ParasailResult(&this, res, s1, s2);
    }

    ParasailResult nw_scan(string s1, string s2)
    {
        auto seq1 = toUTFz!(const char*)(s1);
        auto seq2 = toUTFz!(const char*)(s2);
        auto seq1Len = cast(const int) s1.length;
        auto seq2Len = cast(const int) s2.length;
        auto res = parasail_nw_trace_scan_16(seq1, seq1Len, seq2,
                seq2Len, this.gapOpen, this.gapExt, this.score_matrix);
        return ParasailResult(&this, res, s1, s2);
    }


    ParasailResult aligner(string alg, string output_option = "trace",
            string impl_option = "striped", string sol_width = "sat")(string s1, string s2)
    {
         auto seq1 = toUTFz!(const char*)(s1);
        auto seq2 = toUTFz!(const char*)(s2);
        auto seq1Len = cast(const int) s1.length;
        auto seq2Len = cast(const int) s2.length;
        mixin("auto res = parasail_" ~ alg ~ "_" ~ output_option ~ "_" ~ impl_option ~ "_"
                ~ sol_width ~ "(seq1, seq1Len, seq2, seq2Len, this.gapOpen, this.gapExt, this.score_matrix);");
        return ParasailResult(&this, res, s1, s2);
    }

    ParasailResult databaseAligner(string alg, string impl_option = "striped")(string s2)
    {
        auto seq2 = toUTFz!(const char*)(s2);
        auto seq2Len = cast(const int) s2.length;
        assert(!(this.profile is null));
        mixin("auto res = parasail_" ~ alg ~ "_" ~ impl_option
                ~ "_profile_sat(this.profile, seq2, seq2Len, gapOpen, gapExt);");
        return ParasailResult(&this, res, this.s1, s2);
    }

    void close()
    {
        if(!prebuiltMatrix)
            parasail_matrix_free(cast(parasail_matrix*) score_matrix);

        if (!(this.profile is null))
            parasail_profile_free(this.profile);
    }
}