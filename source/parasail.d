module parasail;
import std.meta:AliasSeq;
import std.stdio:FILE;

extern(C):

//Struct declarations
struct parasail_result_t {
    int score;
    int end_query;
    int end_ref;
    int flag;
    void *extra;
}
struct parasail_matrix_t;
struct parasail_cigar_t{
    uint *seq;
    int len;
    int beg_query;
    int beg_ref;
}
//Matrix functions
parasail_matrix_t* parasail_matrix_lookup(char * matrix);
parasail_matrix_t* parasail_matrix_create(const char *alphabet, const int match, const int mismatch);
void parasail_matrix_set_value(parasail_matrix_t *matrix, int row, int col, int value);
void parasail_matrix_free(parasail_matrix_t* matrix);


alias ALIGNMENT_TYPES = AliasSeq!("nw","sg","sg_qb","sg_qe","sg_qx","sg_db","sg_de","sg_dx","sg_qb_de","sg_qe_db","sw");
alias OPTIONS = AliasSeq!("_stats","_table","_rowcol","_scan");
alias VEC_TYPES = AliasSeq!("_striped","_scan","_diag");
alias VEC_BITS = AliasSeq!("_8","_16","_32","_64","_sat");
alias VEC_OPTIONS = AliasSeq!("_sse2_128","_sse41_128","_avx2_256","_altivec_128","_neon_128");

enum FUNC_SIG="(const char * s1,const int s1Len,"~
                "const char *s2,const int s2Len,"~
                "const int open,const int gap,"~
                "const parasail_matrix_t* matrix);\n";

enum FUNC_PRE="parasail_result_t* parasail_";
//Meta-programming to generate all alignment functions
//Non-vectorized, reference implementations.
//parasail_ {nw,sg,sg_qb,sg_qe,sg_qx,sg_db,sg_de,sg_dx,sg_qb_de,sg_qe_db,sw} [_stats] [{_table,_rowcol}] [_scan]
string generateNonVecFuns(){
    string funcs;
    foreach (aln; ALIGNMENT_TYPES)
    {
        funcs~=FUNC_PRE~aln~FUNC_SIG;
    }
    foreach (aln; ALIGNMENT_TYPES)
    {
        funcs~=FUNC_PRE~aln~OPTIONS[0]~FUNC_SIG;
        funcs~=FUNC_PRE~aln~OPTIONS[1]~FUNC_SIG;
        funcs~=FUNC_PRE~aln~OPTIONS[2]~FUNC_SIG;
        funcs~=FUNC_PRE~aln~OPTIONS[3]~FUNC_SIG;

        funcs~=FUNC_PRE~aln~OPTIONS[0]~OPTIONS[1]~FUNC_SIG;
        funcs~=FUNC_PRE~aln~OPTIONS[0]~OPTIONS[2]~FUNC_SIG;
        funcs~=FUNC_PRE~aln~OPTIONS[0]~OPTIONS[3]~FUNC_SIG;

        funcs~=FUNC_PRE~aln~OPTIONS[0]~OPTIONS[1]~OPTIONS[3]~FUNC_SIG;
        funcs~=FUNC_PRE~aln~OPTIONS[0]~OPTIONS[2]~OPTIONS[3]~FUNC_SIG;

        funcs~=FUNC_PRE~aln~OPTIONS[1]~OPTIONS[3]~FUNC_SIG;
        funcs~=FUNC_PRE~aln~OPTIONS[2]~OPTIONS[3]~FUNC_SIG;
    }
    return funcs;
}
mixin(generateNonVecFuns);

//Non-vectorized, traceback-capable reference implementations.
//parasail_ {nw,sg,sg_qb,sg_qe,sg_qx,sg_db,sg_de,sg_dx,sg_qb_de,sg_qe_db,sw} _trace [_scan]
string generateNonVecTraceFuns(){
    string funcs;
    foreach (aln; ALIGNMENT_TYPES)
    {
        funcs~=FUNC_PRE~aln~"_trace"~FUNC_SIG;
    }
    foreach (aln; ALIGNMENT_TYPES)
    {
        funcs~=FUNC_PRE~aln~"_trace"~"_scan"~FUNC_SIG;
    }
    return funcs;
}
mixin(generateNonVecTraceFuns);

//Vectorized.
//parasail_ {nw,sg,sg_qb,sg_qe,sg_qx,sg_db,sg_de,sg_dx,sg_qb_de,sg_qe_db,sw} [_stats] [{_table,_rowcol}] {_striped,_scan,_diag} [{_sse2_128,_sse41_128,_avx2_256,_altivec_128,_neon_128}] {_8,_16,_32,_64,_sat}
string generateVecFuns(){
    string funcs;
    foreach (aln; ALIGNMENT_TYPES)
    {
        foreach (type; VEC_TYPES)
        {
            foreach (bit; VEC_BITS)
            {
                funcs~=FUNC_PRE~aln~type~bit~FUNC_SIG;
                funcs~=FUNC_PRE~aln~OPTIONS[0]~type~bit~FUNC_SIG;
                funcs~=FUNC_PRE~aln~OPTIONS[1]~type~bit~FUNC_SIG;
                funcs~=FUNC_PRE~aln~OPTIONS[2]~type~bit~FUNC_SIG;

                funcs~=FUNC_PRE~aln~OPTIONS[0]~OPTIONS[1]~type~bit~FUNC_SIG;
                funcs~=FUNC_PRE~aln~OPTIONS[0]~OPTIONS[2]~type~bit~FUNC_SIG;   
            }
        }
    }
    foreach (aln; ALIGNMENT_TYPES)
    {
        foreach (type; VEC_TYPES)
        {
            foreach (bit; VEC_BITS)
            {
                foreach (vec; VEC_OPTIONS)
                {
                    funcs~=FUNC_PRE~aln~type~vec~bit~FUNC_SIG;    
                    funcs~=FUNC_PRE~aln~OPTIONS[0]~type~vec~bit~FUNC_SIG;
                    funcs~=FUNC_PRE~aln~OPTIONS[1]~type~vec~bit~FUNC_SIG;
                    funcs~=FUNC_PRE~aln~OPTIONS[2]~type~vec~bit~FUNC_SIG;

                    funcs~=FUNC_PRE~aln~OPTIONS[0]~OPTIONS[1]~type~vec~bit~FUNC_SIG;
                    funcs~=FUNC_PRE~aln~OPTIONS[0]~OPTIONS[2]~type~vec~bit~FUNC_SIG;
                }
            }
        }
    }
    return funcs;
}
mixin(generateVecFuns);

//Vectorized, traceback-capable.
//parasail_ {nw,sg,sg_qb,sg_qe,sg_qx,sg_db,sg_de,sg_dx,sg_qb_de,sg_qe_db,sw} _trace {_striped,_scan,_diag} [{_sse2_128,_sse41_128,_avx2_256,_altivec_128,_neon_128}] {_8,_16,_32,_64,_sat}
string generateVecTraceFuns(){
    string funcs;
    foreach (aln; ALIGNMENT_TYPES)
    {
        foreach (type; VEC_TYPES)
        {
            foreach (bit; VEC_BITS)
            {
                funcs~=FUNC_PRE~aln~"_trace"~type~bit~FUNC_SIG;   
            }
        }
    }
    foreach (aln; ALIGNMENT_TYPES)
    {
        foreach (type; VEC_TYPES)
        {
            foreach (bit; VEC_BITS)
            {
                foreach (vec; VEC_OPTIONS)
                {
                    funcs~=FUNC_PRE~aln~"_trace"~type~vec~bit~FUNC_SIG;
                }
            }
        }
    }
    return funcs;
}
mixin(generateVecTraceFuns);

mixin(FUNC_PRE~"_nw_banded"~FUNC_SIG);

void parasail_result_free(parasail_result_t *result);

//Cigar Functions
void parasail_cigar_free(parasail_cigar_t *cigar);
char parasail_cigar_decode_op(uint cigar_int);
uint parasail_cigar_decode_len(uint cigar_int);
char* parasail_cigar_decode(parasail_cigar_t *cigar);
parasail_cigar_t* parasail_result_get_cigar(
        parasail_result_t *result,
        const char *seqA, int lena,
        const char *seqB, int lenb,
        const parasail_matrix_t *matrix);

int parasail_result_get_score(parasail_result_t *result);
int parasail_result_get_end_query(parasail_result_t *result);
int parasail_result_get_end_ref(parasail_result_t *result);
int parasail_result_get_matches(parasail_result_t *result);
int parasail_result_get_similar(parasail_result_t *result);
int parasail_result_get_length(parasail_result_t *result);
int* parasail_result_get_score_table(parasail_result_t *result);
int* parasail_result_get_matches_table(parasail_result_t *result);
int* parasail_result_get_similar_table(parasail_result_t *result);
int* parasail_result_get_length_table(parasail_result_t *result);
int* parasail_result_get_score_row(parasail_result_t *result);
int* parasail_result_get_matches_row(parasail_result_t *result);
int* parasail_result_get_similar_row(parasail_result_t *result);
int* parasail_result_get_length_row(parasail_result_t *result);
int* parasail_result_get_score_col(parasail_result_t *result);
int* parasail_result_get_matches_col(parasail_result_t *result);
int* parasail_result_get_similar_col(parasail_result_t *result);
int* parasail_result_get_length_col(parasail_result_t *result);
int* parasail_result_get_trace_table(parasail_result_t *result);
int* parasail_result_get_trace_ins_table(parasail_result_t *result);
int* parasail_result_get_trace_del_table(parasail_result_t *result);

int parasail_result_is_nw(parasail_result_t *result);
int parasail_result_is_sg(parasail_result_t *result);
int parasail_result_is_sw(parasail_result_t *result);
int parasail_result_is_saturated(parasail_result_t *result);
int parasail_result_is_banded(parasail_result_t *result);
int parasail_result_is_scan(parasail_result_t *result);
int parasail_result_is_striped(parasail_result_t *result);
int parasail_result_is_diag(parasail_result_t *result);
int parasail_result_is_blocked(parasail_result_t *result);
int parasail_result_is_stats(parasail_result_t *result);
int parasail_result_is_stats_table(parasail_result_t *result);
int parasail_result_is_stats_rowcol(parasail_result_t *result);
int parasail_result_is_table(parasail_result_t *result);
int parasail_result_is_rowcol(parasail_result_t *result);
int parasail_result_is_trace(parasail_result_t *result);

void parasail_traceback_generic(
        const char *seqA, int lena,
        const char *seqB, int lenb,
        const char *nameA,
        const char *nameB,
        const parasail_matrix_t *matrix,
        parasail_result_t *result,
        char match, /* character to use for a match */
        char pos,   /* character to use for a positive-value mismatch */
        char neg,   /* character to use for a negative-value mismatch */
        int width,  /* width of traceback to display before wrapping */
        int name_width,
        int use_stats); /* if 0, don't display stats, if non-zero, summary stats displayed */

void parasail_traceback_generic_extra(
        const char *seqA, int lena,
        const char *seqB, int lenb,
        const char *nameA,
        const char *nameB,
        const parasail_matrix_t *matrix,
        parasail_result_t *result,
        char match, /* character to use for a match */
        char pos,   /* character to use for a positive-value mismatch */
        char neg,   /* character to use for a negative-value mismatch */
        int width,  /* width of traceback to display before wrapping */
        int name_width,
        int use_stats, /* if 0, don't display stats, if non-zero, summary stats displayed */
        int int_width, /* width used for reference and query indexes */
        FILE *stream); /* to print to custom file stream */