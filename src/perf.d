module perf;

@safe:

import core.time : MonoTime;
import std.stdio : writeln;

Perf newPerf() {
    return Perf(MonoTime.currTime);
}

struct Perf {
	MonoTime before;

	void mark(scope immutable string name) scope {
		immutable MonoTime now = MonoTime.currTime;
		writeln(name, " took ", now - before);
		before = now;
	}
}
