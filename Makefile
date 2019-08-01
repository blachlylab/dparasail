all: parasail/build/libparasail.a

parasail/src/cigar.c:
	git submodule init
	git submodule update

parasail/build/libparasail.a: parasail/src/cigar.c
	cd parasail;mkdir build;cd build;cmake -DBUILD_SHARED_LIBS=OFF ..;make -j8;
