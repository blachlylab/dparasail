module dparasail.memory;

import std.traits : isPointer, isSomeFunction, ReturnType, isSafe;
import core.lifetime : move;
import std.typecons : RefCounted, RefCountedAutoInitialize;
import std.traits : Unqual;
import parasail;
import core.stdc.stdlib : calloc, free;

/// can we use @live for scope checking? 
enum dip1000Enabled = isSafe!((int x) => *&x);

static if(dip1000Enabled)
    pragma(msg, "Using -dip1000 for scope checking and safety");

/// Template struct that performs reference
/// counting on htslib pointers and destroys with specified function
template ParasailMemory(T)
{
    enum freeMix = T.stringof ~ "_free";
    mixin("alias destroy = "~freeMix~";");
    static assert(isSomeFunction!destroy && is(ReturnType!destroy == void));
    
    struct SafeParasailPtr(T)
    if(!isPointer!T)
    {
        @safe @nogc nothrow:

        /// data pointer
        T * ptr;
        /// reference counting
        int* refct;
        /// initialized?
        bool initialized;

        /// ctor that respects scope
        this(T * rawPtr) @trusted return scope
        {
            this.ptr = rawPtr;
            this.refct = cast(int *) calloc(int.sizeof,1);
            (*this.refct) = 1;
            this.initialized = true;
        }
        
        /// postblit that respects scope
        this(this) @trusted return scope
        {
            if(initialized)(*this.refct)++;
        }

        /// allow SafeHtslibPtr to be used as 
        /// underlying ptr type
        alias getRef this;

        /// get underlying data pointer
        @property nothrow pure @nogc
        ref inout(T*) getRef() inout return
        {
            return ptr;
        }

        /// take ownership of underlying data pointer
        @property nothrow pure @nogc
        T* moveRef()
        {
            T * ptr;
            move(this.getRef, ptr);
            return ptr;
        }
        ~this() @trusted return scope
        {
            /// if destroy function return is void 
            /// just destroy
            /// else if int
            /// destroy then check return value 
            /// else don't compile
            if(!this.initialized) return;
            if(--(*this.refct)) return;
            if(this.ptr){
                free(this.refct);
                destroy(this.ptr);
            }
        }
    }
    alias ParasailMemory = SafeParasailPtr!T;
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