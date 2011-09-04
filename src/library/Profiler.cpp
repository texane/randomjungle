#include "Profiler.h"

Profiler::Profiler() {
	isTicking = false;
	timeStamp = 0;
	totalTime = 0;
}

Profiler::~Profiler() {
}

void Profiler::start() {
	// Watch is not running yet
	assert(!isTicking);
	isTicking = true;

	// get time stamp
	timeStamp = clock();
}

void Profiler::stop() {
	// watch is ticking
	assert(isTicking);
	isTicking = false;

	// stop the watch and add time
	totalTime += difftime(clock(), timeStamp)/1000;
}

