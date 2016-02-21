#ifndef _UNIVERSE_CUH_
#define _UNIVERSE_CUH_

#include <vector>

struct Object
{
	float m;
	float3 p, v;
	
	Object() :
		m(0.0f), p({0}), v({0}) 
	{ }
};

struct Container
{
	float total_mass;
	float radius;
	float3 logical_center;
	float3 mass_center;
	int quads[8];
	
	Container() :
		total_mass(-1.0f), logical_center{0}, mass_center{0}, quads{0}
	{ }
};

struct UniverseConfiguration
{
	float time_step;
	float threshold_angle;
	float graviational_constant;
	float softening_factor;
	float max_velocity;
};

class Universe
{
	public:
	
	Universe(UniverseConfiguration c, size_t expected_objects = 0);
	Universe(size_t expected_objects = 0);
	
	void timeStep();
	
	template<typename Z>
	void addObject(Z&& o)
	{
		objects.push_back(o);
	}
	
	Object* getObjects()
	{
		return objects.data();
	}
	
	const Object* getObjects() const
	{
		return objects.data();
	}
	
	size_t size() const
	{
		return objects.size();
	}
	
	UniverseConfiguration& getConfiguration()
	{
		return config;
	}
	
	const UniverseConfiguration& getConfiguration() const
	{
		return config;
	}
	
	float3 getCenterOfMass() const
	{
		if(containers.empty())
			return {0};
		else
			return containers[0].mass_center;
	}
	
	private:
	
	void makeTree();
	void makeRoot();
	void computeComs();
	void reorderObjects();
	void computeTimeStep();
	void reorderObjectsRecursive(int c, int current_depth, std::vector<Object>& objs, std::vector<Container>& ctrs);
	static int quadrant(float3 origin, float3 position);
	
	void checkTree();
	
	UniverseConfiguration config;
	std::vector<Object> objects;
	std::vector<Container> containers;
	int tree_depth;
};

#endif
