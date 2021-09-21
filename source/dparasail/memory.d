module dparasail.memory;

import std.traits : isPointer, isSomeFunction, ReturnType;
import core.lifetime : move;
import std.typecons : RefCounted, RefCountedAutoInitialize;
import std.traits : Unqual;
import parasail;
import htslib.hts_log;

/// Template struct that performs reference
/// counting on htslib pointers and destroys with specified function
template ParasailMemory(T)
{
    enum freeMix = T.stringof ~ "_free";
    mixin("alias destroy = "~freeMix~";");
    static assert(isSomeFunction!destroy && is(ReturnType!destroy == void));
    struct ParasailMemoryImpl(T)
    if(!isPointer!T)
    {
        @safe:
        /// Pointer Wrapper
        static struct ParasailPtr
        {
            /// data pointer
            T * ptr;

            /// no copying this as that could result
            /// in premature destruction
            @disable this(this);

            /// destroy 
            ~this() @trusted
            {
                /// if destroy function return is void 
                /// just destroy
                /// else if int
                /// destroy then check return value 
                /// else don't compile
                if(this.ptr){
                    destroy(this.ptr);
                }
            }
        }

        /// reference counted HtslibPtr
        RefCounted!(ParasailPtr, RefCountedAutoInitialize.yes) rcPtr;

        /// get underlying data pointer
        @property nothrow pure @nogc
        ref inout(T*) getPtr() inout return
        {
            return rcPtr.refCountedPayload.ptr;
        }

        /// allow ParasailMemory to be used as 
        /// underlying ptr type
        alias getPtr this;

        /// ctor from raw pointer
        this(T * rawPtr) @trusted
        {
            auto wrapped = ParasailPtr(rawPtr);
            move(wrapped,this.rcPtr.refCountedPayload);
        }
    }
    alias ParasailMemory = ParasailMemoryImpl!T;
}
/// reference counted bam1_t wrapper
/// can be used directly as a bam1_t *
alias ParasailMatrix = ParasailMemory!(parasail_matrix_t);

/// reference counted bam_hdr_t wrapper
/// can be used directly as a bam_hdr_t *
alias ParasailProfile = ParasailMemory!(parasail_profile_t);

/// reference counted bam_hdr_t wrapper
/// can be used directly as a bam_hdr_t *
alias ParasailResultPtr = ParasailMemory!(parasail_result_t);

unittest
{
    import std.range : chunks;
    auto dna = parasail_matrix_lookup("dnafull");
    auto rc1 = ParasailMatrix(parasail_matrix_copy(dna));
    parasail_matrix_set_value(rc1, 2, 3, 0);
    assert(rc1 != null);
    // No more allocation, add just one extra reference count
    auto rc2 = rc1;
    // Reference semantics
    parasail_matrix_set_value(rc1, 2, 3, 0);

    import std.stdio;
    writeln(rc1.matrix[0..rc1.size * rc1.size].chunks(16));
    assert(rc1.size == 16);
    assert(rc1.matrix[2*rc1.size + 3] == 0);
}