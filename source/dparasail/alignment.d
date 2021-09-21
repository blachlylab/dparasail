module dparasail.alignment;

import std.utf;
import core.stdc.stdlib : exit;

import dparasail.result;
import dparasail.memory;
import parasail;

import htslib.hts_log;

/*
* Acts as a profile of settings for reusing 
* alignment settings
*/
struct Parasail
{
    /// alignment scoring matrix reference
    private ParasailMatrix score_matrix;

    /// using database query?
    private bool databaseQuery;

    /// gap open penalty (must be positive)
    private int gapOpen;
    
    /// gap extension penalty (must be positive)
    private int gapExt;

    /// parasail profile reference if doing
    /// if doing database alignment
    private ParasailProfile profile;

    /// query sequence if doing database alignment
    private string s1;

    @disable this();

    /// use a prebuilt parasail scoring matrix
    /// or an alphabet matrix
    /// and a gap and ext penalty
    this(string matrix, int open, int ext, int match = -1, int mismatch = 1)
    {
        this(matrix, "", open, ext, match, mismatch);
    }

    /// use a prebuilt parasail scoring matrix
    /// or an alphabet matrix
    /// with a database sequence
    /// and a gap and ext penalty
    this(string matrix, string databaseSequence, int open, int ext, int match = -1, int mismatch = 1){
        auto m = parasail_matrix_lookup(toUTFz!(char*)(matrix));

        if(m != null){
            auto mCopy = parasail_matrix_copy(m);
            this(mCopy, open, ext);
        }else{
            hts_log_warning(__FUNCTION__, "Couldn't load matrix named " ~matrix ~", assuming alphabet");
            if(m == null && (match == -1 || mismatch == 1)){
                hts_log_error(__FUNCTION__, "Couldn't load matrix and match and mismatch penalties are unset (-1, 1)");
                exit(-1);
            }
            auto alphabet = matrix;
            hts_log_warning(__FUNCTION__, "Attempting to create matrix from alphabet " ~alphabet);

            assert(match > 0, "match score must be > 0");
            assert(mismatch < 0, "mismatch penalty/score must be < 0");

            auto mCreate = parasail_matrix_create(toUTFz!(char*)(alphabet), match, mismatch);
            if(mCreate == null)
            {
                hts_log_error(__FUNCTION__, "Couldn't create matrix from alphabet " ~alphabet);
                exit(-1);
            }
            this(mCreate, databaseSequence, open, ext);
        }
    }

    /// use a parasail_matrix_t directly
    /// with a database sequence
    /// and a gap and ext penalty
    this(parasail_matrix_t * matrix, string databaseSequence, int open, int ext)
    {
        if(databaseSequence != "")
        {
            auto prof = parasail_profile_create_sat(toUTFz!(char*)(s1), cast(int) s1.length, matrix);
            if(prof == null)
            {
                hts_log_error(__FUNCTION__, "Couldn't create database profile from " ~databaseSequence);
            }
            this.profile = ParasailProfile(prof);
        }
        this(matrix, open, ext);
    }

    /// use a parasail_matrix_t directly
    /// and a gap and ext penalty
    this(parasail_matrix_t * matrix, int open, int ext)
    {
        this.profile.rcPtr.refCountedStore.ensureInitialized;
        this.score_matrix = ParasailMatrix(matrix);
        assert(open > 0, "gap open penalty must be greater than 0");
        assert(ext > 0, "gap extension penalty must be greater than 0");
        this.gapOpen = open;
        this.gapExt = ext;
    }

    /// get scoring matrix
    auto scoreMatrix()
    {
        return this.score_matrix;
    }

    /// run sw alignment
    ParasailResult sw_striped(string s1, string s2)
    {
        if(!(this.profile is null))
            hts_log_warning(__FUNCTION__, "You are using sw align but a database profile is set");

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
        if(!(this.profile is null))
            hts_log_warning(__FUNCTION__, "You are using nw align but a database profile is set");

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
        if(!(this.profile is null))
            hts_log_warning(__FUNCTION__, "You are using aligner but a database profile is set");

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
        if(this.profile is null){
            hts_log_error(__FUNCTION__, "Cannot use databaseAligner when profile is not set");
            exit(-1);
        }
        mixin("auto res = parasail_" ~ alg ~ "_" ~ impl_option
                ~ "_profile_sat(this.profile, seq2, seq2Len, gapOpen, gapExt);");
        return ParasailResult(&this, res, this.s1, s2);
    }
}