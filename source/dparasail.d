module dparasail;
import std.stdio;
import std.conv;
import std.utf;
// import bio.std.hts.bam.cigar;
import dhtslib.sam.cigar;
import std.algorithm:map,filter;
import std.algorithm.iteration:sum;
public import parasail;

/*
* As a note nothing about this struct is thread-safe
* so just don't
*
*/
struct parasail_query
{
    Parasail * alnSettings;
    char * seq1;
    char * seq2;
    int seq1Len;
    int seq2Len;
    //0-based
    int beg_query;
    int beg_ref;
    Cigar cigar;
    parasail_result_t* result;
    Cigar get_cigar(){
        parasail_cigar_t* cigar;
        Cigar cigar_string;
        cigar=parasail_result_get_cigar(result,seq1,seq1Len,seq2,seq2Len,alnSettings.score_matrix);
        beg_query=cigar.beg_query;
        beg_ref=cigar.beg_ref;
        cigar_string = Cigar(((cast(CigarOp*) cigar.seq)[0..cigar.len]).dup);
        parasail_cigar_free(cigar);
        //if *
        if(cigar_string.length==0){
            return cigar_string;
        }
        //if 30I8M make 30S8M
        if(cigar_string[0].op==Ops.INS){
            //move beg_query as well as it is accounted for
            cigar_string[0]=CigarOp(cigar_string[0].length+this.beg_query,Ops.SOFT_CLIP);
            this.beg_query+=cigar_string[0].length;
        }
        //else if 30D8M make 8M and move ref start
        else if(cigar_string[0].op==Ops.DEL){
            this.beg_ref=this.beg_ref+cigar_string[0].length;
            cigar_string=Cigar(cigar_string[1..$]);
            //cigar_string[0]=CigarOperation(cigar_string[0].length+this.beg_query,'S');
        }
        //else if begin query not 0 add softclip
        else if(this.beg_query!=0){
            assert(cigar_string[0].op!=Ops.SOFT_CLIP);
            cigar_string=Cigar(CigarOp(this.beg_query,Ops.SOFT_CLIP)~cigar_string[]);
            this.beg_query=cigar_string[0].length;
        }
        ///////////////////////////////////////////////////////////
        int q_bases_covered=cast(int) cigar_string[].filter!(x=>x.is_query_consuming()).map!(x=>x.length).sum;
        if(cigar_string[$-1].op==Ops.INS){
            cigar_string[$-1]=CigarOp(cigar_string[$-1].length+this.seq1Len-q_bases_covered,Ops.SOFT_CLIP);
            q_bases_covered=cigar_string[].filter!(x=>x.is_query_consuming()).map!(x=>x.length).sum;
        }
        else if(cigar_string[$-1].op==Ops.DEL){
            cigar_string=Cigar(cigar_string[0..($-1)]);
        }
        else if(q_bases_covered!=this.seq1Len){
            cigar_string=Cigar(cigar_string[]~CigarOp(this.seq1Len-q_bases_covered,Ops.SOFT_CLIP));
            q_bases_covered=cigar_string[].filter!(x=>x.is_query_consuming()).map!(x=>x.length).sum;
        }
        assert(q_bases_covered==seq1Len);
        // parasail_cigar_free(cigar);  
        // Bug caused bc Cigar ctor takes reference not copy
        // Cigar implicit dtor destroys array ptr which cleans up our problem
        return cigar_string;
    }
    void writeAlignment(){
        parasail_traceback_generic(
            seq1, seq1Len, 
            seq2, seq2Len, 
            toUTFz!(char *)("Query:"), toUTFz!(char *)("Target:"), 
            alnSettings.score_matrix, result, 
            '|', '*', '*', 60, 7, 1);
    }

    void close(){
        parasail_result_free(result);
    }
}
struct Parasail{
    const(parasail_matrix_t)* score_matrix = null;
    int gap;
    int open;
    parasail_profile_t * profile;
    string s1;

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

    this(string matrix, string s1, int open, int gap){
        this.score_matrix=parasail_matrix_lookup(toUTFz!(char *)(matrix));
        this.profile = parasail_profile_create_sat(toUTFz!(char *)(s1),cast(int) s1.length,this.score_matrix);
        this.s1 = s1;
        this.open=open;
        this.gap=gap;
    }
    this(string alphabet,string s1, int match,int mismatch, int open, int gap){
        this.score_matrix=parasail_matrix_create(toUTFz!(char *)(alphabet),match,mismatch);
        this.profile = parasail_profile_create_sat(toUTFz!(char *)(s1),cast(int) s1.length,this.score_matrix);
        this.s1 = s1;
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
        p.alnSettings = &this;
        p.cigar=p.get_cigar;
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
        p.alnSettings = &this;
        p.cigar=p.get_cigar;
        return p;
    }
    parasail_result_t* nw_scan(parasail_query p){
        return parasail_nw_trace_scan_16(p.seq1,p.seq1Len,p.seq2,p.seq2Len,open,gap,score_matrix);
    }

    parasail_query aligner(string alg, string output_option="trace", string impl_option="striped", string sol_width="sat")(string s1,string s2){
        return aligner!(alg, output_option, impl_option, sol_width)(toUTFz!(char *)(s1),toUTFz!(char *)(s2),cast(int) s1.length,cast(int) s2.length);
    }
    parasail_query aligner(string alg, string output_option, string impl_option, string sol_width)(char *s1,char * s2, int s1Len,int s2Len){
        parasail_query p;
        p.seq1=s1;
        p.seq2=s2;
        p.seq1Len=s1Len;
        p.seq2Len=s2Len;
        p.result=aligner!(alg, output_option, impl_option, sol_width)(p);
        p.alnSettings = &this;
        p.cigar=p.get_cigar;
        return p;
    }
    parasail_result_t* aligner(string alg, string output_option, string impl_option, string sol_width)(parasail_query p){
        mixin("return parasail_"~alg~"_"~output_option~"_"~impl_option~"_"~sol_width~"(p.seq1,p.seq1Len,p.seq2,p.seq2Len,open,gap,score_matrix);");
    }

    parasail_query databaseAligner(string alg, string impl_option="striped")(string s2){
        assert(!(this.profile is null));
        return databaseAligner!(alg, impl_option)(toUTFz!(char *)(s2),cast(int) s2.length);
    }
    parasail_query databaseAligner(string alg, string impl_option)(char * s2,int s2Len){
        parasail_query p;
        p.seq2=s2;
        p.seq2Len=s2Len;
        p.result=databaseAligner!(alg, impl_option)(p);
        p.alnSettings = &this;
        return p;
    }
    parasail_result_t* databaseAligner(string alg, string impl_option)(parasail_query p){
        mixin("return parasail_"~alg~"_"~impl_option~"_profile_sat(this.profile,p.seq2,p.seq2Len,open,gap);");
    }

    

    void close(){
        parasail_matrix_free(cast(parasail_matrix *)score_matrix);
        if(!(this.profile is null)){
            parasail_profile_free(this.profile);
        }
    }
}

unittest{
    import std.stdio;
    auto p=Parasail("dnafull",3,2);
    auto q=p.sw_striped("GATTA","GACTA");
    
    assert(q.beg_query == 0);
    assert(q.beg_ref == 0);
    assert(q.result.score == 16);
    assert(q.cigar.alignedLength == 5);
    assert(q.cigar.toString=="2=1X2=");
    q.writeAlignment;
    

    q=p.sw_striped("GGCTTCTGATCAGGCTTCT","GACTA");

    assert(q.beg_query == 7);
    assert(q.beg_ref == 0);
    assert(q.result.score == 14);
    assert(q.cigar.alignedLength == 5);
    assert(q.cigar.toString=="7S2=1D1=1I1=7S");
    q.writeAlignment;

    q=p.sw_striped("GACTAA","GGCTTCTGATCAGGCTTCT");

    assert(q.beg_query == 0);
    assert(q.beg_ref == 7);
    assert(q.result.score == 14);
    assert(q.cigar.alignedLength == 5);
    assert(q.cigar.toString == "2=1D1=1I1=1S");
    q.writeAlignment;

    q=p.nw_scan("GATTA","GACTA");

    assert(q.beg_query == 0);
    assert(q.beg_ref == 0);
    assert(q.result.score == 16);
    assert(q.cigar.alignedLength == 5);
    assert(q.cigar.toString=="2=1X2=");
    q.writeAlignment;

    q=p.nw_scan("GGCTTCTGATCAGGCTTCT","GACTA");

    assert(q.beg_query == 0);
    assert(q.beg_ref == 0);
    assert(q.result.score == -14);
    assert(q.cigar.alignedLength == 5);
    assert(q.cigar.toString=="1=1X2=4I1=10S");
    q.writeAlignment;
    

    q=p.nw_scan("GACTAA","GGCTTCTGATCAGGCTTCT");
    assert(q.beg_query == 0);
    assert(q.beg_ref == 0);
    assert(q.result.score == -8);
    assert(q.cigar.alignedLength == 12);
    assert(q.cigar.toString == "1=1X2=4D1=2D1=");
    q.writeAlignment;
    // auto cigar=parasail_result_get_cigar(q.result,q.seq1,q.seq1Len,q.seq2,q.seq2Len,q.score_matrix);
    // writeln(parasail_cigar_decode(cigar));
}

