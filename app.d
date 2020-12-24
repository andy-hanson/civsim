@safe:

import std.algorithm : max, min;
import std.algorithm.iteration : filter, sum;
import std.algorithm.sorting : sort;
import std.array : array;
import std.conv : to;
import std.file : write;
import std.math : isInfinity, pow, round, sqrt;
import std.random : Mt19937, uniform01;
import std.stdio : writeln;

immutable uint D = 3;
void main(scope immutable string[] args) {
	immutable Options options = parseOptions(args[1..$]); // Arg 0 is path to exe
	immutable double speedOfLight = 1.0;
	immutable uint randomSeed = 1337;
	immutable Setup setup = generateSetup(randomSeed, options.nPoints, options.speedOfSpread, options.civOriginPower);

	immutable double[] waitTimes = calcWaitTimes(setup, options.speedOfSpread);
	immutable double[] distanceToClosestAtOrigin = calcDistanceToClosestAtOrigin(setup);
	immutable VisibleCivs visibleCivs = calcVisibleCivs(setup, options.speedOfSpread, speedOfLight);

	immutable double[] xValues = map!(double, Vec)(setup.civLocations, (scope ref immutable Vec it) =>
		it[0]);
	immutable double[] yValues = map!(double, Vec)(setup.civLocations, (scope ref immutable Vec it) =>
		it[1]);
	immutable double[] zValues = map!(double, Vec)(setup.civLocations, (scope ref immutable Vec it) =>
		it[2]);

	writeln(showOptions(options));
	writeln("x values: ", summarize(xValues, [1, 10]));
	writeln("origin times: ", summarize(setup.civOriginTimes, [1, 10]));
	writeln("wait times: ", summarize(waitTimes, [1, 10]));
	writeln("distance to closest: ", summarize(distanceToClosestAtOrigin, [1, 10]));
	writeln("count visible: ", summarize(visibleCivs.visibleCount, [90, 99]));
	writeln("biggest angle: ", summarize(visibleCivs.biggestAngle, [90, 99]));

	if (options.writeCsv) {
		immutable string csv = toCsv(
			["x", "y", "z", "origin time", "wait time", "closest", "n visible", "biggest angle"],
			[xValues, yValues, zValues, setup.civOriginTimes, waitTimes, distanceToClosestAtOrigin, visibleCivs.visibleCount, visibleCivs.biggestAngle]);
		immutable string csvName = "results_" ~ showOptions(options) ~ ".csv";
		write(csvName, csv);
		writeln("Wrote to ", csvName);
	}
}

immutable(string) toCsv(
	scope immutable string[] columnNames,
	scope immutable double[][] columns,
) {
	assert(columnNames.length == columns.length);
	string res;
	writeWithCommas!string(res, columnNames, (scope ref immutable string it) {
		res ~= it;
	});
	res ~= '\n';
	foreach (immutable size_t row; 0..columns[0].length) {
		writeWithCommas!(double[])(res, columns, (scope ref immutable double[] column) {
			res ~= column[row].to!string;
		});
		res ~= '\n';
	}
	return res;
}

void writeWithCommas(T)(
	scope ref string res,
	scope ref immutable T[] values,
	scope void delegate(scope ref immutable T) @safe cb,
) {
	bool first = true;
	foreach (ref immutable T value; values) {
		if (first)
			first = false;
		else
			res ~= ',';
		cb(value);
	}
}

immutable(string) summarize(scope ref immutable double[] values, scope immutable uint[] percentiles) {
	immutable double[] filteredValues = array(filter!((immutable double s) => !isInfinity(s))(values));
	scope immutable double[] sortedValues = array(sort(filteredValues.dup));
	string res = sortedValues.length.to!string ~ " values"
		~ ", mean=" ~ mean(sortedValues).to!string
		~ ", median=" ~ getPercentile(sortedValues, 50).to!string
		~ ", min=" ~ sortedValues[0].to!string
		~ ", max=" ~ sortedValues[$ - 1].to!string;
	foreach (immutable uint percentile; percentiles)
		res ~= ", " ~ toOrdinal(percentile) ~ " percentile=" ~ getPercentile(sortedValues, percentile).to!string;
	return res;
}

immutable(double) getPercentile(scope ref immutable double[] sortedValues, uint pct) {
	assert(sortedValues.length > 1);
	return sortedValues[cast(size_t) round((sortedValues.length - 1) * pct / 100.0)];
}

struct Options {
	immutable uint nPoints;
	immutable double speedOfSpread;
	immutable double civOriginPower;
	immutable bool writeCsv;
}

immutable(string) showOptions(scope ref immutable Options options) {
	return "nPoints=" ~ options.nPoints.to!string
		~ "~speedOfSpread=" ~ options.speedOfSpread.to!string
		~ "~civOriginPower=" ~ options.civOriginPower.to!string;
}

immutable(Options) parseOptions(scope immutable string[] args) {
	uint nPoints = 1000;
	double speedOfSpread = 0.5;
	double civOriginPower = 6.0;
	bool writeCsv;
	for (uint i = 0; i < args.length;) {
		immutable string arg = args[i];
		i++;
		if (arg == "--nPoints") {
			nPoints = args[i].to!uint;
			i++;
		} else if (arg == "--speedOfSpread") {
			speedOfSpread = arg[i].to!double;
			i++;
		} else if (arg == "civOriginPower") {
			civOriginPower = arg[i].to!double;
			i++;
		} else if (arg == "--writeCsv") {
			writeCsv = true;
		} else {
			writeln("Unexpected argument ", arg);
		}
	}
	return immutable Options(nPoints, speedOfSpread, civOriginPower, writeCsv);
}

pure nothrow:

immutable(double[]) calcWaitTimes(scope ref immutable Setup setup, immutable double speedOfSpread) {
	return fillArr!double(setup.nCivs, (immutable size_t i) =>
		minOver(setup.nCivs, (immutable size_t j) {
			immutable double ti = setup.civOriginTimes[i];
			immutable double tj = setup.civOriginTimes[j];
			immutable double dt = ti - tj;
			if (dt > 0) {
				immutable double dist = distance(setup.civLocations[i], setup.civLocations[j]);
				immutable double thisWaitTime = dist / speedOfSpread - dt;
				assert(thisWaitTime > 0);
				return thisWaitTime;
			} else
				return double.infinity;
		}));
}

immutable(double[]) calcDistanceToClosestAtOrigin(scope ref immutable Setup setup) {
	return fillArr!double(setup.nCivs, (immutable size_t i) =>
		minOver(setup.nCivs, (immutable size_t j) =>
			setup.civOriginTimes[j] < setup.civOriginTimes[i]
				? distance(setup.civLocations[i], setup.civLocations[j])
				: double.infinity));
}

struct VisibleCivs {
	immutable double[] biggestAngle;
	immutable double[] visibleCount;
}
immutable(VisibleCivs) calcVisibleCivs(
	scope ref immutable Setup setup,
	immutable double speedOfSpread,
	immutable double speedOfLight,
) {
	double[] biggestAngles;
	double[] visibleCounts;

	assert(speedOfSpread < speedOfLight);
	foreach (immutable size_t i; 0..setup.nCivs) {
		uint count = 0;
		immutable double biggest = maxOver(setup.nCivs, (immutable size_t j) {
			if (setup.civOriginTimes[j] < setup.civOriginTimes[i]) {
				immutable double dist = distance(setup.civLocations[i], setup.civLocations[j]);
				immutable double firstVisibleTime = setup.civOriginTimes[j] + dist / speedOfLight;
				immutable double dt = setup.civOriginTimes[i] - firstVisibleTime;
				if (dt > 0) {
					count++;
					return speedOfSpread * dt / dist;
				}
			}
			return -double.infinity;
		});
		biggestAngles ~= biggest;
		visibleCounts ~= count;
	}
	return immutable VisibleCivs(castImmutable(biggestAngles), castImmutable(visibleCounts));
}

struct Setup {
	@safe @nogc pure nothrow:

	Vec[] civLocations;
	double[] civOriginTimes;

	immutable(size_t) nCivs() immutable {
		return civLocations.length;
	}
}

immutable(Setup) generateSetup(
	immutable uint randomSeed,
	immutable uint nPoints,
	immutable double speedOfSpread,
	immutable double civOriginPower,
) {
	Random random;
	random.seed(randomSeed);
	immutable Vec[] allLocations = fillArr!Vec(nPoints, (immutable size_t) => randomVec(random));
	immutable double[] allOriginTimeR = fillArr!double(nPoints, (immutable size_t) => random.uniform01!double());
	immutable double[] allOriginTimesUnsorted = map!(double, double)(allOriginTimeR, (scope ref immutable double r) =>
		pow(r, 1.0 / (1.0 + civOriginPower)));
	immutable double[] allOriginTimes = array(sort(allOriginTimesUnsorted.dup));

	Vec[] civLocations;
	double[] civOriginTimes;
	foreach (immutable uint i; 0..nPoints)
		if (!isPrecluded(allLocations[0..i], allOriginTimes[0..i], speedOfSpread, allLocations[i], allOriginTimes[i])) {
			civLocations ~= allLocations[i];
			civOriginTimes ~= allOriginTimes[i];
		}
	return immutable Setup(castImmutable(civLocations), castImmutable(civOriginTimes));
}

immutable(bool) isPrecluded(
	scope immutable Vec[] locations,
	scope immutable double[] originTimes,
	immutable double speedOfSpread,
	scope ref immutable Vec location,
	immutable double time,
) {
	foreach (immutable size_t j; 0..originTimes.length)
		if (originTimes[j] + timeToSpread(locations[j], location, speedOfSpread) <= time)
			return true;
	return false;
}

immutable(Vec) randomVec(ref Random random) {
	Vec res;
	foreach (immutable uint i; 0..D)
		res[i] = random.uniform01!double();
	return res;
}

immutable(double) wrappingDistance(immutable double a, immutable double b) {
	immutable double low = min(a, b);
	immutable double high = max(a, b);
	return min(high - low, low + 1.0 - high);
}

immutable(double) distanceSquared(ref const Vec a, ref const Vec b) {
	double res = 0.0;
	foreach (immutable uint i; 0..D) {
		immutable double d = wrappingDistance(a[i], b[i]);
		res += d * d;
	}
	return res;
}

immutable(double) distance(ref const Vec a, ref const Vec b) {
	return sqrt(distanceSquared(a, b));
}

double timeToSpread(ref const Vec a, ref const Vec b, immutable double speedOfSpread) {
	return distance(a, b) / speedOfSpread;
}

alias Random = Mt19937;
alias Vec = double[D];

@trusted immutable(T[]) castImmutable(T)(T[] values) {
	return cast(immutable) values;
}

immutable(Out[]) map(Out, In)(
	scope ref immutable In[] values,
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
