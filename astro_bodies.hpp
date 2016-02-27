#ifndef _ASTRO_BODIES_HPP_
#define _ASTRO_BODIES_HPP_

#include <random>
#include <cmath>
#include "Universe.hpp"
#include <vector>

void galaxy_generate(Universe& u, int n, float max_r, float max_m, float min_m, float3 p0 = {0}, float3 v0 = {0})
{
	static std::default_random_engine gen;
	static bool seeded = false;
	if(!seeded)
	{
		gen.seed(time(0));
		seeded = true;
	}
	
	std::uniform_real_distribution<float> theta_dist(0, 2.0f * 3.14159256f);
	std::uniform_real_distribution<float> r_dist(0, max_r);
	std::uniform_real_distribution<float> z_dist(0, 0.01 * max_r);
	std::uniform_real_distribution<float> m_dist(min_m, max_m);
	std::uniform_real_distribution<float> voff_dist(0.9, 1.1);
	
	std::vector<Object> objects;
	
	double radial_masses[100] = {0};
	for(int i = 0; i < n; i++)
	{
		float r = r_dist(gen);
		float theta = theta_dist(gen);
		float z = z_dist(gen);
		float m = m_dist(gen);
		
		Object o;
		o.p.x = r;
		o.p.y = theta;
		o.p.z = z * voff_dist(gen);
		o.m = m;
		
		int r_bucket = 100.0f * r / max_r;
		radial_masses[r_bucket] += m;
		
		objects.push_back(o);
	}
	
	for(Object& o : objects)
	{
		double m_enc = 0.0;
		int r_bucket = 100.0f * o.p.x / max_r;
		for(int i = 0; i < r_bucket; i++)
			m_enc += radial_masses[i];
		
		float r = o.p.x;
		float theta = o.p.y;
		
		float v = sqrtf(u.getConfiguration().graviational_constant * m_enc / (r + 10));
		float s = sinf(theta);
		float c = cosf(theta);
		
		o.p.x = c * r + p0.x;
		o.p.y = s * r + p0.y;
		o.p.z += p0.z;
		o.v.x = -s * v * voff_dist(gen) + v0.x;
		o.v.y = c * v * voff_dist(gen) + v0.y;
		o.v.z  = v0.z * voff_dist(gen);
		u.addObject(o);
	}
}

#endif
