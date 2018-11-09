import std.stdio;
import dparasail;
void main()
{
	writeln("Edit source/app.d to start your project.");
	auto p=Parasail("dnafull",3,2);
	p.sw_striped("GGGATCA","GACTAAA");
	writeln(p.result.score);
	writeln(p.get_cigar());
}
