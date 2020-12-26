@safe:

import std.algorithm : max, min;
import std.algorithm.iteration : filter;
import std.algorithm.sorting : sort;
import std.array : array;
import std.conv : to;
import std.file : write;
import std.math : isInfinity, round;
import std.parallelism : taskPool, TaskPool;
import std.stdio : writeln;

import perf : newPerf, Perf;
import sim : distance, runSim, Setup, SimOptions, SimResult, Vec;
import test : test;
import util : castImmutable, castMutable, mean, toOrdinal;

void main(scope immutable string[] args) {
	//TODO: clear and remake 'results' directory on each run

	// immutable Options options = parseOptions(args[1..$]); // Arg 0 is path to exe
	TaskPool pool = taskPool();
	test();
	loop(pool);
}

void loop(ref TaskPool taskPool) {
	immutable uint D = 1;
	writeln("D=", D);

	writeln("nPoints,speedOfSpread,nCivs,medianCivOrigin,medianWaitTime");
	for (uint nPoints = 16_384; nPoints <= 2_097_152; nPoints *= 2) {
		for (double speedOfSpread = 0.125; speedOfSpread <= 0.5; speedOfSpread *= 2) {
			immutable SimOptions options = immutable SimOptions(nPoints, speedOfSpread, 6.0);
			immutable SimResult!D result = runSim!D(taskPool, options);
			immutable size_t nCivs = result.nCivs;
			immutable double medianCivOrigin = getMedian(result.originTimes);
			immutable double medianWaitTime = getMedian(result.waitTimes);
			string temp;
			writeWithCommas(temp, [double(nPoints), speedOfSpread, double(nCivs), medianCivOrigin, medianWaitTime]);
			writeln(temp);
			writeCsv("results/" ~ showOptions!D(options) ~ ".csv", result);
		}
	}
}

void printSimResult(uint D)(scope ref immutable SimResult!D result) {
	writeln("coord[0] values: ", summarize(result.coordValues[0], [1, 10]));
	writeln("origin times: ", summarize(result.originTimes, [1, 10]));
	writeln("wait times: ", summarize(result.waitTimes, [1, 10]));
	writeln("distance to closest: ", summarize(result.distanceToClosestAtOrigin, [1, 10]));
	writeln("count visible: ", summarize(result.visibleCount, [90, 99]));
	writeln("biggest angle: ", summarize(result.biggestAngle, [90, 99]));
}

void writeCsv(uint D)(scope immutable string path, scope ref immutable SimResult!D result) {
	scope immutable string[] columnNames = ["x", /*"y", "z",*/ "origin time", "wait time", "closest", "n visible", "biggest angle"];
	scope double[][D + 5] columns;
	foreach (immutable uint d; 0..D)
		columns[d] = castMutable(result.coordValues[d]);
	columns[D + 0] = castMutable(result.originTimes);
	columns[D + 1] = castMutable(result.waitTimes);
	columns[D + 2] = castMutable(result.distanceToClosestAtOrigin);
	columns[D + 3] = castMutable(result.visibleCount);
	columns[D + 4] = castMutable(result.biggestAngle);
	write(path, toCsv(columnNames, castImmutable(columns)));
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

void writeWithCommas(scope ref string res, scope immutable double[] values) {
	writeWithCommas!double(res, values, (scope ref immutable double it) {
		res ~= it.to!string;
	});
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

immutable(double) getMedian(scope ref immutable double[] values) {
	return getPercentile(array(sort(values.dup)), 50);
}

immutable(double) getPercentile(scope immutable double[] sortedValues, uint pct) {
	assert(sortedValues.length > 1);
	return sortedValues[cast(size_t) round((sortedValues.length - 1) * pct / 100.0)];
}

immutable(string) showOptions(uint D)(scope ref immutable SimOptions options) {
	return "D=" ~ D.to!string 
		~ "~nPoints=" ~ options.nPoints.to!string
		~ "~speedOfSpread=" ~ options.speedOfSpread.to!string
		~ "~civOriginPower=" ~ options.civOriginPower.to!string;
}

immutable(SimOptions) parseOptions(scope immutable string[] args) {
	uint nPoints = 1000;
	double speedOfSpread = 0.5;
	double civOriginPower = 6.0;
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
		} else {
			writeln("Unexpected argument ", arg);
		}
	}
	return immutable SimOptions(nPoints, speedOfSpread, civOriginPower);
}
