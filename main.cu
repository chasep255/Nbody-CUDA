#include "Universe.hpp"
#include "timer.hpp"
#include <iostream>
#include <GL/glut.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <thread>
#include "astro_bodies.hpp"
#include <cassert>

#define N 30000

Universe u;
int2 window_dim;
float3 com;
double xangle = 0.0;
double yangle = 0.0;
double zoom = 1.0;
float* vertex_buffer = nullptr;
float* color_buffer = nullptr;
float alpha = 1.0f;
bool paused = true;

void display()
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (float)window_dim.x / window_dim.y, 50, 10000);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0, 0.0, 1000.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	glScaled(zoom, zoom, zoom);
	glRotated(xangle, 1.0, 0.0, 0.0);
	glRotated(yangle, 0.0, 1.0, 0.0);
	
	double alpha_norm = M_2_PI * atan(alpha);
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, vertex_buffer);
	glColorPointer(4, GL_FLOAT, 0, color_buffer);
	glDrawArrays(GL_POINTS, 0, u.size());
	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
	
	
//	glColor3f(1.0, 0.0, 0.0);
//	glutWireCube(2.0 * u.getConfiguration().max_distance_from_com);
	
	glutSwapBuffers();
}

void update_universe_display()
{
	if(vertex_buffer == nullptr)
		vertex_buffer = new float[3 * u.size()];
	
	if(color_buffer == nullptr)
		color_buffer = new float[4 * u.size()];
	
	Object* objects = u.getObjects();
	
	#pragma omp parallel for
	for(int i = 0; i < u.size(); i++)
	{
		reinterpret_cast<float3*>(vertex_buffer)[i] =  objects[i].p;
		
		color_buffer[4 * i + 0] = (1000.0 - objects[i].m) / 1000.0;
		color_buffer[4 * i + 1] = 0.55;
		color_buffer[4 * i + 2] = objects[i].m / 1000.0;
		color_buffer[4 * i + 3] = alpha;
	}
}

void create_universe()
{
	srand(time(0));
	
	UniverseConfiguration cfg;
	cfg.graviational_constant = 1000.0;
	cfg.max_velocity = 1e12;
	cfg.softening_factor = 1;
	cfg.threshold_angle = 15;
	cfg.time_step = 0.02;
	cfg.max_distance_from_com = 1e9;
	u.setConfiguration(cfg);
	
	auto rand3 = [](float max) 
	{
		float3 v;
		v.x = max * rand() / RAND_MAX; 
		v.y = max * rand() / RAND_MAX; 
		v.z = 0; 
		return v;
	};
	
	galaxy_generate(u, N / 2, 1500, 1000, 1, rand3(500000), rand3(5));
	galaxy_generate(u, N / 3, 1000, 1000, 1, rand3(500000), rand3(5));
	galaxy_generate(u, N, 3000, 1000, 1, rand3(500000), rand3(5));
	galaxy_generate(u, N, 3000, 1000, 1, rand3(500000), rand3(5));
	galaxy_generate(u, N, 3000, 1000, 1, rand3(500000), rand3(5));
	
	update_universe_display();
}
bool record = false;
void save()
{
	static int count = 0;
	static int pn = 0;
	
	if(!record)
		return;
	count++;
	
	if(count % 15 == 0)
	{
		count = 0;
		char fn[1000];
		sprintf(fn, "convert /dev/stdin /media/chase/3161D67803D8C5BE/galaxy/img%06d.jpg", pn++);
		std::cout << fn << std::endl;
		FILE* f = popen(fn, "w");
		assert(f);
		
		fprintf(f, "P6 1280 720 255 ");
		
		unsigned char* pix = new unsigned char[3 * 1000 * 1000];
		assert(pix);
		glReadPixels(0, 0, 1280, 720, GL_RGB, GL_UNSIGNED_BYTE, pix);
		
		fwrite(pix, 1000 * 1000, 3, f);
		delete[] pix;
		fclose(f);
		
		if(pn >= 30000)
			exit(0);
	}
}

void idle()
{
	static bool update = true;
	
	if(!paused && update)
	{
		update_universe_display();
		update = false;
		std::thread([&]()
		{
			u.timeStep();
			update = true;
		}).detach();
		
		save();
		glutPostRedisplay();
	}
	
	if(paused)
	{
		update = true;
	}
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitWindowSize(1280, 720);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutCreateWindow("N-BODY");
	
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_BLEND);
	glPointSize(1.0f);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	glutReshapeFunc([](int w, int h)
	{
		window_dim = {w, h};
		glViewport(0, 0, w, h);
	});
	
	glutKeyboardFunc([](unsigned char k, int x, int y)
	{
		bool redisplay = false;
		switch(k)
		{
			case 'w': xangle += 1.0; redisplay = true; break;
			case 's': xangle -= 1.0; redisplay = true; break;
			case 'a': yangle += 1.0; redisplay = true; break;
			case 'd': yangle -= 1.0; redisplay = true; break;
			case 'z': alpha *= 1.05; redisplay = true; break;
			case 'x': alpha *= 0.95; redisplay = true; break;
			case '=': zoom *= 1.2; redisplay = true; break;
			case '-': zoom /= 1.2; redisplay = true; break;
			case 'r': record = !record; break;
			case ' ': paused = !paused; break;
		}
		
		if(redisplay)
		{
			glutPostRedisplay();
		}
	});
	
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	
	create_universe();
	
	glutMainLoop();
	
}
