module util;

@safe pure nothrow:

import std.algorithm : max, min;
import std.algorithm.iteration : sum;
import std.conv : to;

immutable(Out[]) map(Out, In)(
	scope immutable In[] values,
	scope immutable(Out) delegate(scope ref immutable In) @safe pure nothrow cb,
) {
	Out[] res;
	foreach (ref immutable In value; values)
		res ~= cb(value);
	return castImmutable(res);
}

immutable(T[]) fillArr(T)(
	immutable size_t size,
	scope immutable(T) delegate(immutable size_t) @safe pure nothrow cb,
) {
	T[] res;
	foreach (immutable size_t i; 0..size)
		res ~= cb(i);
	return res;
}

@trusted immutable(T[size]) fillArrN(T, size_t size)(scope immutable(T) delegate(immutable size_t) @safe pure nothrow cb) {
	T[size] res;
	foreach (immutable size_t i; 0..size)
		res[i] = cast(T) cb(i);
	return cast(immutable) res;
}

immutable(double) mean(scope ref immutable double[] values) {
	return sum(values) / double(values.length);
}

immutable(double) minOver(
	immutable size_t n,
	scope immutable(double) delegate(immutable size_t) @safe pure nothrow cb,
) {
	double res = double.infinity;
	foreach (immutable size_t i; 0..n)
		res = min(res, cb(i));
	return res;
}

immutable(double) maxOver(
	immutable size_t n,
	scope immutable(double) delegate(immutable size_t) @safe pure nothrow cb,
) {
	double res = -double.infinity;
	foreach (immutable size_t i; 0..n)
		res = max(res, cb(i));
	return res;
}

immutable(string) toOrdinal(immutable uint i) {
	return i == 1 ? "1st" : (i.to!string ~ "th");
}

@trusted immutable(T[]) castImmutable(T)(T[] values) {
	return cast(immutable) values;
}

@trusted T[] castMutable(T)(immutable T[] values) {
	return cast(T[]) values;
}
