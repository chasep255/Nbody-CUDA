#include <chrono>

#ifndef _TIMER_HPP_
#define _TIMER_HPP_

inline double timer_now()
{
	auto t = std::chrono::high_resolution_clock::now().time_since_epoch();
	return 1.0e-6 * std::chrono::duration_cast<std::chrono::microseconds>(t).count();
}

#endif
