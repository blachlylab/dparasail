module dparasail;
import std.stdio;
import std.conv;
import std.utf;
import bio.bam.cigar;
import std.algorithm:map,filter;
import std.algorithm.iteration:sum;
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

parasail_result_t* parasail_sw_trace_striped_16(
        const char * s1,const int s1Len,
        const char *s2,const int s2Len,
        const int open,const int gap,
        const parasail_matrix_t* matrix
        );

parasail_result_t* parasail_nw_trace_scan_16(
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
parasail_cigar_t* parasail_result_get_cigar(
        parasail_result_t *result,
        const char *seqA, int lena,
        const char *seqB, int lenb,
        const parasail_matrix_t *matrix);



//D wrapping
extern(D):

/*
* As a note nothing about this struct is thread-safe
* so just don't
*
*/
struct parasail_query{
    char * seq1;
    char * seq2;
    int seq1Len;
    int seq2Len;
    int beg_query;
    int beg_ref;
    parasail_result_t* result;
    //string get_cigar(parasail_matrix_t* score_matrix){
    //    parasail_cigar_t* cigar;
    //    string cigar_string;
    //    int i=0;
    //    cigar=parasail_result_get_cigar(result,seq1,seq1Len,seq2,seq2Len,score_matrix);
    //    beg_query=cigar.beg_query;
    //    beg_ref=cigar.beg_ref;
    //    char letter;
    //    uint length;
    //    int count=0;
    //    if(this.beg_query!=0){
    //        length=this.beg_query;
    //        letter=parasail_cigar_decode_op(cigar.seq[0]);
    //        if(letter=='I'){
    //            length+= parasail_cigar_decode_len(cigar.seq[0]);
    //            i=1;
    //        }else if(letter=='*'){
    //            this.cigar="*";
    //            return;
    //        }
    //        count+=length;
    //        letter='S';
    //        cigar_string=cigar_string~to!string(length)~letter;
    //    }
    //    for (;i<cigar.len; ++i) {
    //        letter=parasail_cigar_decode_op(cigar.seq[i]);
    //        length= parasail_cigar_decode_len(cigar.seq[i]);
    //        if(i==0 && letter=='I'){
    //            letter='S';
    //        }
    //        if(i==0 && letter=='D'){
    //            letter='S';
    //            this.beg_ref+=length;
    //            continue;
    //        }
    //        if(i==(cigar.len-1) && letter=='I'){
    //            letter='S';
    //        }
    //        if(letter!='D'){
    //            count+=length;
    //        }
    //        //if((i==0 && letter=='D')|(i==(cigar.len-1) && letter=='D')){
    //        //    continue;
    //        //}
    //        cigar_string=cigar_string~to!string(length)~letter;
    //    }
    //    letter='S';
    //    length=seq1Len-result.end_query-1;
    //    if (length>0){
    //        count+=length;
    //        cigar_string=cigar_string~to!string(length)~letter;
    //    }
    //    assert(count==seq1Len);
    //    parasail_cigar_free(cigar);
    //    return cigar_string;
    //}
    CigarOperations get_cigar(parasail_matrix_t* score_matrix){
        parasail_cigar_t* cigar;
        CigarOperations cigar_string;
        cigar=parasail_result_get_cigar(result,seq1,seq1Len,seq2,seq2Len,score_matrix);
        beg_query=cigar.beg_query;
        beg_ref=cigar.beg_ref;
        for (auto i=0;i<cigar.len; ++i) {
            cigar_string~=CigarOperation(parasail_cigar_decode_len(cigar.seq[i]),parasail_cigar_decode_op(cigar.seq[i]));
        }
        if(cigar_string.is_unavailable){
            return cigar_string;
        }
        if(cigar_string[0].type=='I'){
            cigar_string[0]=CigarOperation(cigar_string[0].length+this.beg_query,'S');
            this.beg_query=0;
        }
        else if(cigar_string[0].type=='D'){
            this.beg_ref=this.beg_ref+cigar_string[0].length;
            //cigar_string[0]=CigarOperation(cigar_string[0].length+this.beg_query,'S');
        }
        else if(this.beg_query!=0){
            cigar_string=CigarOperation(this.beg_query,'S')~cigar_string;
            this.beg_query=0;
        }
        ///////////////////////////////////////////////////////////
        int q_bases_covered=cast(int) cigar_string.filter!(x=>x.is_query_consuming()).map!(x=>x.length).sum;
        if(cigar_string[$-1].type=='I'){
            cigar_string[$-1]=CigarOperation(cigar_string[$-1].length+this.seq1Len-q_bases_covered,'S');
            q_bases_covered=cigar_string.filter!(x=>x.is_query_consuming()).map!(x=>x.length).sum;
        }
        else if(cigar_string[$-1].type=='D'){
            cigar_string=cigar_string[0..($-1)];
        }
        else if(q_bases_covered!=this.seq1Len){
            cigar_string~=CigarOperation(this.seq1Len-q_bases_covered,'S');
            q_bases_covered=cigar_string.filter!(x=>x.is_query_consuming()).map!(x=>x.length).sum;
        }
        assert(q_bases_covered==seq1Len);
        parasail_cigar_free(cigar);
        return cigar_string;
    }
    void close(){
        parasail_result_free(result);
    }
}
struct Parasail{
    parasail_matrix_t* score_matrix = null;
    int gap;
    int open;
    @disable this();
    this(string matrix, int open, int gap){
        this.score_matrix=parasail_matrix_lookup(toUTFz!(char *)(matrix));
        this.open=open;
        this.gap=gap;
    }
    this(string alphabet,int match,int mismatch, int open, int gap){
        this.score_matrix=parasail_matrix_create(toUTFz!(char *)(alphabet),match,mismatch);
        this.open=open;
        this.gap=gap;
    }
    parasail_query sw_striped(string s1,string s2){
        return sw_striped(toUTFz!(char *)(s1),toUTFz!(char *)(s2),cast(int) s1.length,cast(int) s2.length);
    }
    parasail_query sw_striped(char *s1,char * s2, int s1Len,int s2Len){
        parasail_query p;
        p.seq1=s1;
        p.seq2=s2;
        p.seq1Len=s1Len;
        p.seq2Len=s2Len;
        p.result=sw_striped(p);
        //p.cigar=p.get_cigar(this.score_matrix);
        return p;
    }
    parasail_result_t* sw_striped(parasail_query p){
        return parasail_sw_trace_striped_16(p.seq1,p.seq1Len,p.seq2,p.seq2Len,open,gap,score_matrix);
    }
    parasail_query nw_scan(string s1,string s2){
        return nw_scan(toUTFz!(char *)(s1),toUTFz!(char *)(s2),cast(int) s1.length,cast(int) s2.length);
    }
    parasail_query nw_scan(char *s1,char * s2, int s1Len,int s2Len){
        parasail_query p;
        p.seq1=s1;
        p.seq2=s2;
        p.seq1Len=s1Len;
        p.seq2Len=s2Len;
        p.result=nw_scan(p);
        return p;
    }
    parasail_result_t* nw_scan(parasail_query p){
        return parasail_nw_trace_scan_16(p.seq1,p.seq1Len,p.seq2,p.seq2Len,open,gap,score_matrix);
    }
    void close(){
        parasail_matrix_free(score_matrix);
    }
}

unittest{
    import std.stdio;
    auto p=Parasail("dnafull",3,2);
    auto q=p.sw_striped("GATTA","GACTA");
    writeln(q.seq1[0..q.seq1Len]);
    writeln(q.cigar);
    writeln(q.beg_query);
    writeln(q.beg_ref);
    assert(q.cigar=="2=1X2=");
    q=p.sw_striped("GGCTTCTGATCAGGCTTCT","GACTA");
    writeln(q.seq1[0..q.seq1Len]);
    writeln(q.cigar);
    writeln(q.beg_query);
    writeln(q.beg_ref);
    assert(q.cigar=="7S2=1D1=1I1=7S");
    q=p.sw_striped("GACTAA","GGCTTCTGATCAGGCTTCT");
    writeln(q.seq1[0..q.seq1Len]);
    writeln(q.cigar);
    writeln(q.beg_query);
    writeln(q.beg_ref);
    q=p.nw_scan("GATTA","GACTA");
    writeln(q.seq1[0..q.seq1Len]);
    writeln(q.cigar);
    writeln(q.beg_query);
    writeln(q.beg_ref);
    assert(q.cigar=="2=1X2=");
    q=p.nw_scan("GGCTTCTGATCAGGCTTCT","GACTA");
    writeln(q.seq1[0..q.seq1Len]);
    writeln(q.cigar);
    writeln(q.beg_query);
    writeln(q.beg_ref);
    assert(q.cigar=="1=1X2=4I1=10S");
    q=p.nw_scan("GACTAA","GGCTTCTGATCAGGCTTCT");
    writeln(q.seq1[0..q.seq1Len]);
    writeln(q.cigar);
    writeln(q.beg_query);
    writeln(q.beg_ref);
    auto cigar=parasail_result_get_cigar(q.result,q.seq1,q.seq1Len,q.seq2,q.seq2Len,q.score_matrix);
    writeln(parasail_cigar_decode(cigar));
}

