#include <sys/resource.h>
#include <unistd.h>
#include <iostream>

#include "MemoryTracker.h"

size_t MemoryTracker::getPeakMemory()
{
	struct rusage usage;
	if (getrusage(RUSAGE_SELF, &usage) == 0)
		return static_cast<size_t>(usage.ru_maxrss) * 1024;
		
	return 0;
}

size_t MemoryTracker::getCurrentMemory()
{
	long pages = 0;
	FILE* file = fopen("/proc/self/statm", "r");
	if (file)
	{
		if (fscanf(file, "%*s %ld", &pages) != 1) pages = 0;
		fclose(file);
	}

	return static_cast<size_t>(pages) * sysconf(_SC_PAGESIZE);
}
