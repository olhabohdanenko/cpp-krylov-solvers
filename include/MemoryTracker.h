// MemoryTracker.h
#include <sys/resource.h>
#include <unistd.h>
#include <iostream>

class MemoryTracker
{
public:
	static size_t getPeakMemory();

	static size_t getCurrentMemory();
};
