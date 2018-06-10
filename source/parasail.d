module dparasail;
import std.conv;
import std.utf;
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
parasail_matrix_t* parasail_matrix_create(
        const char *alphabet, const int match, const int mismatch);
void parasail_matrix_free(parasail_matrix_t* matrix);

//SW functions
parasail_cigar_t* parasail_result_get_cigar(
        parasail_result_t *result,
        const char *seqA, int lena,
        const char *seqB, int lenb,
        const parasail_matrix_t *matrix);


parasail_result_t* parasail_sw_trace_striped_sat(
        const char * s1,const int s1Len,
        const char *s2,const int s2Len,
        const int open,const int gap,
        const parasail_matrix_t* matrix
        );
void parasail_result_free(parasail_result_t *result);

//Cigar Functions
void parasail_cigar_free(parasail_cigar_t *cigar);
char parasail_cigar_decode_op(uint cigar_int);
uint parasail_cigar_decode_len(uint cigar_int);
char* parasail_cigar_decode(parasail_cigar_t *cigar);


//D wrapping
extern(D):

struct Parasail{
    parasail_matrix_t* score_matrix = null;
    parasail_result_t* result =null;
    int gap;
    int open;
    char * seq1;
    char * seq2;
    int seq1Len;
    int seq2Len;
    this(string matrix, int open, int gap){
        this.score_matrix=parasail_matrix_lookup(toUTFz!(char *)(matrix));
        this.open=open;
        this.gap=gap;
    }
    void sw_striped(string s1,string s2){
        seq1=toUTFz!(char *)(s1);
        seq2=toUTFz!(char *)(s2);
        seq1Len=cast(int) s1.length;
        seq2Len=cast(int) s2.length;
        sw_striped();
    }
    void sw_striped(char *s1,char * s2, int s1Len,int s2Len){
        seq1=s1;
        seq2=s2;
        seq1Len=s1Len;
        seq2Len=s2Len;
        sw_striped();
    }
    void sw_striped(){
        if (result!=null){
            parasail_result_free(result);
        }
        result=parasail_sw_trace_striped_sat(seq1,seq1Len,seq2,seq2Len,open,gap,score_matrix);
    }
    string get_cigar(){
        parasail_cigar_t* cigar;
        string cigar_string;
        int i=0;
        cigar=parasail_result_get_cigar(result,seq1,seq1Len,seq2,seq2Len,score_matrix);
        for (i=0; i<cigar.len; ++i) {
            char letter=parasail_cigar_decode_op(cigar.seq[i]);
            uint length = parasail_cigar_decode_len(cigar.seq[i]);
            cigar_string=cigar_string~to!string(length)~letter;
        }
        parasail_cigar_free(cigar);
        return cigar_string;
    }
    void close(){
        parasail_result_free(result);
        parasail_matrix_free(score_matrix);
    }
}

unittest{
    auto p=Parasail("dnafull",3,2);
    p.sw_striped("GATCA","GACTA");
    assert(p.get_cigar()=="2=1D1=1I1=");
}

