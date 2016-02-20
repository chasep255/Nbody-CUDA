#include "Universe.cuh"
#include <cfloat>
#include <cmath>
#include <iostream>
#include <cassert>
#include <utility>

inline void cucheck()
{
	cudaError_t err = cudaGetLastError();
	if(err != cudaSuccess)
	{
		throw std::runtime_error(cudaGetErrorString(err));
	}
}


Universe::Universe(UniverseConfiguration c, size_t expected_objects) :
	config(c) 
{
	objects.reserve(expected_objects);
	containers.reserve(expected_objects);
}

int Universe::quadrant(float3 origin, float3 position)
{
	float dx = origin.x - position.x;
	float dy = origin.y - position.y;
	float dz = origin.z - position.z;
	
	return ((dz < 0.0f) << 2) | ((dy < 0.0f) << 1) | (dx < 0.0f);
}

void Universe::makeTree()
{
	for(int o = 0; o < objects.size(); o++)
	{
		Object& obj = objects[o];
		
		int ctr = -1;
		while(true)
		{
			int q = quadrant(containers[-ctr - 1].logical_center, obj.p);
			if(containers[-ctr - 1].quads[q] == 0)
			{
				containers[-ctr - 1].quads[q] = o + 1;
				break;
			}
			else if(containers[-ctr - 1].quads[q] > 0)
			{
				int found_obj = containers[-ctr - 1].quads[q];
				
				float3 foundp = objects[found_obj - 1].p;
				if(foundp.x == obj.p.x && foundp.y == obj.p.y && foundp.z == obj.p.z)
					break;
				
				Container new_container;
				new_container.radius = containers[-ctr - 1].radius * 0.5;
				new_container.logical_center = containers[-ctr - 1].logical_center;
				new_container.logical_center.z += (q >> 2) ? new_container.radius : -new_container.radius;
				new_container.logical_center.y += ((q >> 1) & 1) ? new_container.radius : -new_container.radius;
				new_container.logical_center.x += (q & 1) ? new_container.radius : -new_container.radius;	
				new_container.quads[quadrant(new_container.logical_center, objects[found_obj - 1].p)] = found_obj;
				containers.push_back(new_container);
				containers[-ctr - 1].quads[q] = -containers.size();
				ctr = -containers.size();
			}
			else
			{
				ctr = containers[-ctr - 1].quads[q];
			}
		}
	}
}

void Universe::computeComs()
{
	for(int i = containers.size() - 1; i >= 0; i--)
	{
		Container& c = containers[i];
		
		double total_mass = 0.0;
		double3 com = {0};
		for(int q = 0; q < 8; q++)
		{
			if(c.quads[q] > 0)
			{
				Object& o = objects[c.quads[q] - 1];
				total_mass += o.m;
				com.x += o.p.x * o.m;
				com.y += o.p.y * o.m;
				com.z += o.p.z * o.m;
			}
			else if(c.quads[q] < 0)
			{
				Container& c2 = containers[-c.quads[q] - 1];
				total_mass += c2.total_mass;
				com.x += c2.mass_center.x * c2.total_mass;
				com.y += c2.mass_center.y * c2.total_mass;
				com.z += c2.mass_center.z * c2.total_mass;
			}
		}
		c.total_mass = total_mass;
		c.mass_center.x = com.x / total_mass;
		c.mass_center.y = com.y / total_mass;
		c.mass_center.z = com.z / total_mass;
	}
}

void Universe::makeRoot()
{
	containers.clear();
	float max_x = FLT_MIN, min_x = FLT_MAX;
	float max_y = FLT_MIN, min_y = FLT_MAX;
	float max_z = FLT_MIN, min_z = FLT_MAX;
	
	for(int o = 0; o < objects.size(); o++)
	{
		float3 p = objects[o].p;
		max_x = fmaxf(p.x, max_x);
		min_x = fminf(p.x, min_x);
		max_y = fmaxf(p.x, max_y);
		min_y = fminf(p.x, min_y);
		max_z = fmaxf(p.x, max_z);
		min_z = fminf(p.x, min_z);
	}
	
	float dx = max_x - min_x;
	float dy = max_y - min_y;
	float dz = max_z - min_z;
	
	float r = fmaxf(dx, fmaxf(dy, dz));
	float cx = 0.5f * (max_x + min_x);
	float cy = 0.5f * (max_y + min_y);
	float cz = 0.5f * (max_z + min_z);
	
	containers.push_back(Container());
	containers[0].logical_center.x = cx;
	containers[0].logical_center.y = cy;
	containers[0].logical_center.z = cz;
	containers[0].radius = r;
}

void Universe::checkTree()
{
	std::cout << "ROOT RADIUS: " << containers[0].radius << std::endl;
	std::cout << "ROOT CENTER: (" << containers[0].logical_center.x << ", " <<
			containers[0].logical_center.y << ", " <<
			containers[0].logical_center.z << ")" << std::endl;
	
	std::cout << "TREE DEPTH: " << tree_depth << std::endl;
	
	int objs_found = 0;
	double total_mass = 0.0;
	for(int i = 0; i < containers.size(); i++)
	{
		Container& c = containers[i];
		for(int q = 0; q < 8; q++)
		{
			float3 p;
			bool check = false;
			if(c.quads[q] < 0)
			{
				p = containers[-c.quads[q] - 1].logical_center;
				check = true;
			}
			else if(c.quads[q] > 0)
			{
				p = objects[c.quads[q] - 1].p;
				total_mass += objects[c.quads[q] - 1].m;
				check = true;
				objs_found++;
			}
			if(check)
			{
				float dx = c.logical_center.x - p.x;
				float dy = c.logical_center.y - p.y;
				float dz = c.logical_center.z - p.z;
				
				assert((dx < 0.0f) == (q & 1));
				assert((dy < 0.0f) == ((q >> 1) & 1));
				assert((dz < 0.0f) == ((q >> 2) & 1));
			}
		}
	}
	
	std::cout << "OBJECTS FOUND: " << objs_found << std::endl;
	std::cout << "EXPECTED TOTAL MASS: " << total_mass << std::endl;
	std::cout << "ROOT TOTAL MASS: " << containers[0].total_mass << std::endl;
}


void Universe::reorderObjects()
{
	std::vector<Object> new_objects;
	std::vector<Container> new_containers;
	
	new_objects.reserve(objects.size());
	new_containers.reserve(new_containers.size());
	
	tree_depth = 0;
	reorderObjectsRecursive(-1, 1, new_objects, new_containers);
	
	objects = std::move(new_objects);
	containers = std::move(new_containers);
	tree_depth += 1;
}

void Universe::reorderObjectsRecursive(int c, int current_depth, std::vector<Object>& objs, std::vector<Container>& ctrs)
{
	if(current_depth > tree_depth)
		tree_depth = current_depth;
	
	Container& ct = containers[-c - 1];
	int new_index = ctrs.size();
	ctrs.push_back(ct);
	
	for(int q = 0; q < 8; q++)
	{
		if(ct.quads[q] < 0)
		{
			ctrs[new_index].quads[q] = -ctrs.size() - 1;
			reorderObjectsRecursive(ct.quads[q], current_depth + 1, objs, ctrs);
		}
	}
	
	for(int q = 0; q < 8; q++)
	{
		if(ct.quads[q] > 0)
		{
			ctrs[new_index].quads[q] = objs.size() + 1;
			objs.push_back(objects[ct.quads[q] - 1]);
		}
	}
}

__global__
void acceleration_kernel(Object* objects, Container* containers, float3* acceleration_buffer, int objects_size, int tree_depth, float threshold_tangent, float g)
{
	extern __shared__ int visit_stack[];
	
	int oid = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if(oid > objects_size)
		return;
	
	int warp_num = threadIdx.x / warpSize;
	int wid = threadIdx.x % warpSize;
	
	int* stack = visit_stack + tree_depth * warp_num;
	int stack_ptr = 1;
	
	if(wid == 0)
		stack[0] = -1;
	
	float ax = 0.0f, ay = 0.0f, az = 0.0f;
	
	float3 op = objects[oid - 1].p;
	float x, y, z, m, r;
	int c, t;
	
	while(stack_ptr)
	{
		int stack_base_ptr;
		if(wid == 0)
		{
			t = stack[--stack_ptr];
			stack_base_ptr = stack_ptr;
			
			if(t < 0)
			{
				Container& ct = containers[-t - 1];
				x = ct.mass_center.x;
				y = ct.mass_center.y;
				z = ct.mass_center.z;
				m = ct.total_mass;
				r = ct.radius;
				c = 1;
				
				for(int q = 0; q < 8; q++)
				{
					if(ct.quads[q])
						stack[stack_ptr++] = ct.quads[q];
				}
			}
			else
			{
				Object& o = objects[t - 1];
				x = o.p.x;
				y = o.p.y;
				z = o.p.z;
				m = o.m;
				r = 0.0f;
				c = 0;
			}
		}
		
		x = __shfl(x, 0);
		y = __shfl(y, 0);
		z = __shfl(z, 0);
		m = __shfl(m, 0);
		r = __shfl(r, 0);
		c = __shfl(c, 0);
		t = __shfl(t, 0);
		
		float dx = x - op.x;
		float dy = y - op.y;
		float dz = z - op.z;
		float over_d = rnorm3df(dx, dy, dz);
		
		if(__all(!c || (2.0f * r * over_d) < threshold_tangent))
		{
			if(oid != t)
			{
				float a = g * m * over_d * over_d;	
				ax += a * dx * over_d;
				ay += a * dy * over_d;
				az += a * dz * over_d;
			}
			
			if(wid == 0)
				stack_ptr = stack_base_ptr;
		}
		
		stack_ptr = __shfl(stack_ptr, 0);
	}
	
	acceleration_buffer[oid - 1].x = ax;
	acceleration_buffer[oid - 1].y = ay;
	acceleration_buffer[oid - 1].z = az;
}

__global__
void time_step_kernel(Object* objects, float3* acceleration_buffer, int objects_size, float dt)
{
	int oid = blockDim.x * blockIdx.x + threadIdx.x;
	if(oid >= objects_size)
		return;
	
	objects[oid].v.x += acceleration_buffer[oid].x * dt;
	objects[oid].v.y += acceleration_buffer[oid].y * dt;
	objects[oid].v.z += acceleration_buffer[oid].z * dt;
	
	objects[oid].p.x += objects[oid].v.x * dt;
	objects[oid].p.y += objects[oid].v.y * dt;
	objects[oid].p.z += objects[oid].v.z * dt;
}

void Universe::computeTimeStep(float dt)
{
	Object* objects_d = nullptr;
	Container* containers_d = nullptr;
	float3* acceleration_buffer_d = nullptr;
	
	auto free_memory = [&]()
	{
		if(objects_d)
			cudaFree(objects_d);
		if(containers_d)
			cudaFree(containers_d);
		if(acceleration_buffer_d)
			cudaFree(acceleration_buffer_d);
	};
	
	try
	{	
		cudaMalloc(&objects_d, sizeof(Object) * objects.size());
		cucheck();
		
		cudaMalloc(&containers_d, sizeof(Container) * containers.size());
		cucheck();
		
		cudaMalloc(&acceleration_buffer_d, sizeof(float3) * objects.size());
		cucheck();
		
		cudaMemcpy(objects_d, objects.data(), sizeof(Object) * objects.size(), cudaMemcpyHostToDevice);
		cucheck();
		
		cudaMemcpy(containers_d, containers.data(), sizeof(Container) * containers.size(), cudaMemcpyHostToDevice);
		cucheck();
		
		int device;
		cudaGetDevice(&device);
		cucheck();
		
		cudaDeviceProp p;
		cudaGetDeviceProperties(&p, device);
		cucheck();
		
		int blocks, threads;
		cudaOccupancyMaxPotentialBlockSize(&blocks, &threads, acceleration_kernel, 0);
		
		blocks = (objects.size() + (threads - 1)) / threads;
		
		acceleration_kernel<<< blocks, threads, (threads / p.warpSize) * tree_depth * sizeof(int) * 8 >>>(objects_d, containers_d, acceleration_buffer_d, 
				objects.size(), tree_depth, tan(config.threshold_angle * 3.14159f / 180.0f), config.graviational_constant);
		cucheck();
		
		cudaOccupancyMaxPotentialBlockSize(&blocks, &threads, time_step_kernel, 0, 0);
		cucheck();
		
		blocks = (objects.size() + (threads - 1)) / threads;
		
		time_step_kernel<<< blocks, threads, 0 >>>(objects_d, acceleration_buffer_d, objects.size(), dt);
		cucheck();
		
		cudaMemcpy(objects.data(), objects_d, sizeof(Object) * objects.size(), cudaMemcpyDeviceToHost);
		cucheck();
		
		cudaDeviceSynchronize();
		cucheck();
	}
	catch(...)
	{
		free_memory();
		throw;
	}
	free_memory();
}

#include "timer.hpp"

void Universe::timeStep(float dt)
{
	double mkrt = timer_now();
	makeRoot();
	mkrt = timer_now() - mkrt;
	
	double mktr = timer_now();
	makeTree();
	mktr = timer_now() - mktr;
	
	double ccms = timer_now();
	computeComs();
	ccms = timer_now() - ccms;
	
	double reorder = timer_now();
	reorderObjects();
	reorder = timer_now() - reorder;
	
	double ts = timer_now();
	computeTimeStep(dt);
	ts = timer_now() - ts;
	
	double total = mkrt + mktr + ccms + reorder + ts;
	
	std::cout << "MAKE ROOT: %" << (100.0 * mkrt / total) << std::endl;
	std::cout << "MAKE TREE: %" << (100.0 * mktr / total) << std::endl;
	std::cout << "CCM: %" << (100.0 * ccms / total) << std::endl;
	std::cout << "REORDER : %" << (100.0 * reorder / total) << std::endl;
	std::cout << "TIME STEP: %" << (100.0 * ts / total) << std::endl;
	std::cout << "TOTAL: " << total << std::endl;
}
