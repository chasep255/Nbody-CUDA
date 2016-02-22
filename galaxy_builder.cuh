#ifndef _GALAXY_BUILDER_HPP_
#define _GALAXY_BUILDER_HPP_

#include <cmath>
#include <random>
#include <vector>
#include "Universe.cuh"


void build_galaxy(Universe& u, int n, float mass_max, float mass_min, float rad)
{
	static std::default_random_engine gen;
	std::uniform_real_distribution<float> r_dist(0, rad);
	std::uniform_real_distribution<float> theta_dist(0, M_PI * 2.0);
	std::uniform_real_distribution<float> m_dist(mass_min, mass_max);
	std::uniform_real_distribution<float> z_dist(-rad / 2000.0, rad / 2000.0);
	
	std::vector<Object> objects;
	double radial_mass_dis[100] = {0};
	for(int i = 0; i < n; i++)
	{
		float r = r_dist(gen);
		float theta = theta_dist(gen);
		float z = z_dist(gen);
		float m = m_dist(gen);
		
		r = fminf(r, rad);
		
		radial_mass_dis[(int)(100.0 * r / rad)] += m;
		
		Object o;
		o.m = m;
		o.p.x = r;
		o.p.y = theta;
		o.p.z = z;
		
		objects.push_back(o);
	}
	
	for(Object& o : objects)
	{
		int r = 100.0 * o.p.x / rad;
		double m_enc = 0.0;
		for(int i = 0; i < r; i++)
		{
			m_enc += radial_mass_dis[i];
		}
		
		
		float v = sqrt(u.getConfiguration().graviational_constant * m_enc / (o.p.x + 1));
		v *= (rad - o.p.x) / rad;
		float x = o.p.x * cos(o.p.y);
		float y = o.p.x * sin(o.p.y);
		
		o.v.x = v * -sin(o.p.y);
		o.v.y = v * cos(o.p.y);
		
		o.p.x = x;
		o.p.y = y;
		
		u.addObject(o);
	}
}

#endif
