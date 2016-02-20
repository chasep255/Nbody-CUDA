#include "Universe.cuh"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include "timer.hpp"

int main()
{
	cudaSetDevice(0);
	srand(time(0));
	UniverseConfiguration cfg;
	cfg.graviational_constant = 1.0f;
	cfg.time_step = 0.1f;
	cfg.threshold_angle = 15.0f;
	
	Universe u(cfg, 1000000);
	for(int i = 0; i < 1000000; i++)
	{
		Object o;
		o.p.x = 1000.0 * rand() / RAND_MAX;
		o.p.y = 1000.0 * rand() / RAND_MAX;
		o.p.z = 1000.0 * rand() / RAND_MAX;
		o.m = 100.0 * rand() / RAND_MAX;
		u.addObject(o);
	}
	double start = timer_now();
	u.timeStep(0.01);
	double end = timer_now();
	
	std::cout << (end - start) << std::endl;
	return 0;
}
