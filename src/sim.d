module sim;

@safe:

import std.algorithm : max, min;
import std.algorithm.sorting : sort;
import std.array : array;
import std.math : pow, round, sqrt;
import std.parallelism : TaskPool;
import std.random : Mt19937, uniform01;
import std.stdio : writeln;
import perf : Perf;
import util : castImmutable, fillArr, fillArrN, map;

struct SimOptions {
	immutable uint nPoints;
	immutable double speedOfSpread;
	immutable double civOriginPower;
}

struct SimResult(uint D) {
	@safe @nogc pure nothrow:

	immutable double[][D] coordValues;
	immutable double[] originTimes;
	immutable double[] waitTimes;
	immutable double[] distanceToClosestAtOrigin;
	immutable double[] biggestAngle;
	immutable double[] visibleCount;

	immutable(size_t) nCivs() immutable {
		return originTimes.length;
	}
}

immutable(SimResult!D) runSim(uint D)(scope ref TaskPool taskPool, scope ref immutable SimOptions options) {
	immutable double speedOfLight = 1.0;
	immutable uint randomSeed = 1337;
	immutable Setup!D setup =
		generateSetup!D(taskPool, randomSeed, options.nPoints, options.speedOfSpread, options.civOriginPower);

	immutable Calculated calculated = calculate(taskPool, setup, options.speedOfSpread, speedOfLight);
	
	immutable double[][D] coordValues = fillArrN!(double[], D)((immutable size_t d) =>
		map!(double, Vec!D)(setup.civLocations, (scope ref immutable Vec!D it) =>
			it[d]));

	return immutable SimResult!D(
		coordValues,
		setup.civOriginTimes,
		calculated.waitTimes,
		calculated.distanceToClosestAtOrigin,
		calculated.biggestAngle,
		calculated.visibleCount);
}

struct Setup(uint D) {
	@safe @nogc pure nothrow:

	Vec!D[] civLocations;
	double[] civOriginTimes;

	immutable(size_t) nCivs() scope immutable {
		return civLocations.length;
	}
}

alias Vec(uint D) = double[D];

pure nothrow immutable(double) distance(uint D)(scope ref immutable Vec!D a, scope ref immutable Vec!D b) {
	return sqrt(distanceSquared!D(a, b));
}

pure nothrow immutable(double) distanceSquared(uint D)(scope ref immutable Vec!D a, scope ref immutable Vec!D b) {
	double res = 0.0;
	assert(a.length == b.length);
	foreach (immutable uint i; 0..a.length) {
		immutable double d = wrappingDistance(a[i], b[i]);
		res += d * d;
	}
	return res;
}

private:

immutable size_t workUnitSize = 64;

@trusted immutable(Setup!D) generateSetup(uint D)(
	scope ref TaskPool taskPool,
	immutable uint randomSeed,
	immutable uint nPoints,
	immutable double speedOfSpread,
	immutable double civOriginPower,
) {
	Random random;
	random.seed(randomSeed);
	immutable Vec!D[] allLocations = fillArr!(Vec!D)(nPoints, (immutable size_t) =>
		randomVec!D(random));
	immutable double[] allOriginTimes = map!(double, double)(randomSortedArray(random, nPoints), (scope ref immutable double r) =>
		pow(r, 1.0 / (1.0 + civOriginPower)));
	//writeln("allOriginTimes: ", allOriginTimes);

	bool[] precluded = new bool[nPoints];
	foreach (size_t i, ref immutable Vec!D location; taskPool.parallel(allLocations, workUnitSize))
		precluded[i] = isPrecluded!D(allLocations[0..i], allOriginTimes[0..i], speedOfSpread, location, allOriginTimes[i]);
	//writeln("precluded: ", precluded);

	Vec!D[] civLocations;
	double[] civOriginTimes;
	foreach (immutable size_t i; 0..nPoints)
		if (!precluded[i]) {
			civLocations ~= allLocations[i];
			civOriginTimes ~= allOriginTimes[i];
		}
	return immutable Setup!D(castImmutable(civLocations), castImmutable(civOriginTimes));
}

immutable(double[]) randomSortedArray(ref Random random, immutable size_t size) {
	return array(sort(fillArr!double(size, (immutable size_t) =>
		random.uniform01!double()).dup));

}

struct Calculated {
	immutable double[] waitTimes;
	immutable double[] distanceToClosestAtOrigin;
	immutable double[] biggestAngle;
	immutable double[] visibleCount;
}

@trusted immutable(Calculated) calculate(uint D)(
	scope ref TaskPool taskPool,
	scope ref immutable Setup!D setup,
	immutable double speedOfSpread,
	immutable double speedOfLight,
) {
	double[] waitTimes = new double[setup.nCivs];
	double[] distanceToClosestAtOrigin = new double[setup.nCivs];
	double[] biggestAngles = new double[setup.nCivs];
	double[] visibleCounts = new double[setup.nCivs];

	assert(speedOfSpread < speedOfLight);

	foreach (size_t i, immutable double ti; taskPool.parallel(setup.civOriginTimes, workUnitSize)) {
		immutable Vec!D locationI = setup.civLocations[i];
		double waitTime = double.infinity;
		double distanceToClosest = double.infinity;
		double biggestAngle = -double.infinity;
		uint countVisible = 0;
		foreach (immutable size_t j; 0..setup.nCivs) {
			immutable double tj = setup.civOriginTimes[j];
			immutable double dt = ti - tj;
			if (dt > 0) {
				immutable double dist = distance!D(locationI, setup.civLocations[j]);
				immutable double thisWaitTime = dist / speedOfSpread - dt;
				assert(thisWaitTime > 0);
				waitTime = min(waitTime, thisWaitTime);
				distanceToClosest = min(distanceToClosest, dist);

				immutable double firstVisibleTime = tj + dist / speedOfLight;
				immutable double visibleDt = ti - firstVisibleTime;
				if (visibleDt > 0) {
					countVisible++;
					biggestAngle = max(biggestAngle, speedOfSpread * visibleDt / dist);
				}
			}
		}
		waitTimes[i] = waitTime;
		distanceToClosestAtOrigin[i] = distanceToClosest;
		biggestAngles[i] = biggestAngle;
		visibleCounts[i] = countVisible;
	}

	return immutable Calculated(
		castImmutable(waitTimes),
		castImmutable(distanceToClosestAtOrigin),
		castImmutable(biggestAngles),
		castImmutable(visibleCounts));
}

pure nothrow:

alias Random = Mt19937;

public immutable(double) wrappingDistance(immutable double a, immutable double b) {
	immutable double low = min(a, b);
	immutable double high = max(a, b);
	return min(high - low, low + 1.0 - high);
}

public double timeToSpread(uint D)(scope ref immutable Vec!D a, scope ref immutable Vec!D b, immutable double speedOfSpread) {
	return distance!D(a, b) / speedOfSpread;
}

public immutable(bool) isPrecluded(uint D)(
	scope immutable Vec!D[] prevLocations,
	scope immutable double[] prevOriginTimes,
	immutable double speedOfSpread,
	scope ref immutable Vec!D location,
	immutable double originTime,
) {
	assert(prevLocations.length == prevOriginTimes.length);
	foreach (immutable size_t j; 0..prevOriginTimes.length)
		if (isPrecludedBy!D(prevLocations[j], prevOriginTimes[j], speedOfSpread, location, originTime))
			return true;
	return false;
}

public immutable(bool) isPrecludedBy(uint D)(
	scope ref immutable Vec!D prevLocation,
	immutable double prevOriginTime,
	immutable double speedOfSpread,
	scope ref immutable Vec!D location,
	immutable double originTime,
) {
	assert(prevOriginTime <= originTime);
	immutable double timeToReach = timeToReach!D(prevLocation, prevOriginTime, speedOfSpread, location);
	return timeToReach <= originTime;
}

public immutable(double) timeToReach(uint D)(
	scope ref immutable Vec!D prevLocation,
	immutable double prevOriginTime,
	immutable double speedOfSpread,
	scope ref immutable Vec!D location,
) {
	return prevOriginTime + timeToSpread!D(prevLocation, location, speedOfSpread);
}

immutable(Vec!D) randomVec(uint D)(ref Random random) {
	Vec!D res;
	foreach (ref double cur; res)
		cur = random.uniform01!double();
	return res;
}
