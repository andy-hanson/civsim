module test;

@safe pure:

import sim : distance, distanceSquared, isPrecluded, isPrecludedBy, Vec, timeToReach, timeToSpread, wrappingDistance;
import std.math : abs;
import std.stdio : writeln;

void test() {
	testWrappingDistance();
	testDistance();
	testTimeToSpread();
	testTimeToReach();
	testIsPrecludedBy();
	testIsPrecluded();
}

private:

immutable Vec!3 v000 = [0.0, 0.0, 0.0];
immutable Vec!3 v111 = [1.0, 1.0, 1.0];
immutable Vec!3 vMid00 = [0.5, 0.0, 0.0];
immutable Vec!3 vMid = [0.5, 0.5, 0.5];

void testWrappingDistance() {
	assertEq(0.1, wrappingDistance(0.05, 0.95));
	assertEq(0.1, wrappingDistance(0.95, 0.05));
	assertEq(0.4, wrappingDistance(0.55, 0.95));
	assertEq(0.4, wrappingDistance(0.95, 0.55));
}

void testDistance() {
	assertEq(0.0, distance!3(v000, v111));
	assertEq(0.75, distanceSquared!3(v000, vMid));
	assertEq(0.5, distance!3(v000, vMid00));
}

void testTimeToSpread() {
	assertEq(1.0, timeToSpread!3(v000, vMid00, 0.5));
	assertEq(0.5, timeToSpread!3(v000, vMid00, 1.0));
	assertEq(0.25, timeToSpread!3(v000, vMid00, 2.0));
}

void testTimeToReach() {
	assertEq(2.23, timeToReach!3(v000, 1.23, 0.5, vMid00));
}

void testIsPrecludedBy() {
	immutable Vec!1 a = [0.0];
	immutable Vec!1 b = [0.5];
	assert(!isPrecludedBy!1(a, 0.0, 0.99, b, 0.5));
	assert(isPrecludedBy!1(a, 0.0, 1.01, b, 0.5));
	assert(!isPrecludedBy!1(b, 0.0, 0.99, a, 0.5));
	assert(isPrecludedBy!1(b, 0.0, 1.01, a, 0.5));

	assert(!isPrecludedBy!3(v000, 0.1, 0.5, vMid00, 1.09));
	assert(isPrecludedBy!3(v000, 0.1, 0.5, vMid00, 1.1));
	assert(isPrecludedBy!3(v000, 0.1, 0.5, vMid00, 1.11));
}

void testIsPrecluded() {
	scope immutable Vec!3[] prevLocations = [v000, v111];
	scope immutable double[] prevOriginTimes = [0.1, 0.2];
	immutable double speedOfSpread = 0.5;
	
	assert(!isPrecluded!3(prevLocations, prevOriginTimes, speedOfSpread, vMid00, 1.09));
	assert(isPrecluded!3(prevLocations, prevOriginTimes, speedOfSpread, vMid00, 1.1));
}

void assertEq(immutable double expected, immutable double actual) {
	if (abs(actual - expected) > 0.001) {
		debug {
			writeln("Expected ", expected, " == ", actual);
		}
		assert(false);
	}
}
