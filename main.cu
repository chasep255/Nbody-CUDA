#include "Universe.cuh"
#include "timer.hpp"
#include <iostream>
#include <GL/glut.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <thread>

#define N 200000

Universe u;
int2 window_dim;
float3 com;
double xangle = 0.0;
double yangle = 0.0;
double zoom = 1.0;
float* vertex_buffer = nullptr;
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
	
	glColor4f(1.0f, 1.0f, 1.0f, alpha_norm);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, vertex_buffer);
	glDrawArrays(GL_POINTS, 0, N);
	glDisableClientState(GL_VERTEX_ARRAY);
	
	
	glColor3f(1.0, 0.0, 0.0);
	glutWireCube(2.0 * u.getConfiguration().max_distance_from_com);
	
//	for(Container& c : u.ctrs())
//	{
//		glColor3f(1, 1, 1);
//		glPushMatrix();
//		glTranslated(c.logical_center.x,c.logical_center.y, c.logical_center.z);
//		glutWireCube(c.radius * 2);
//		glPopMatrix();
//	}
	
	glutSwapBuffers();
}

void update_universe_display()
{
	if(vertex_buffer == nullptr)
		vertex_buffer = new float[3 * N];
	
	Object* objects = u.getObjects();
	
	#pragma omp parallel for
	for(int i = 0; i < N; i++)
	{
		reinterpret_cast<float3*>(vertex_buffer)[i] =  objects[i].p;
	}
}

void create_universe()
{
	srand(time(0));
	
	UniverseConfiguration cfg;
	cfg.graviational_constant = 1.0;
	cfg.max_velocity = 1e9;
	cfg.softening_factor = 0.1;
	cfg.threshold_angle = 30;
	cfg.time_step = 0.01;
	cfg.max_distance_from_com = 1000;
	u.setConfiguration(cfg);
	
	for(int i = 0; i < N; i++)
	{
		Object obj;
		obj.m = 100;
		obj.p.x = 500.0 * ((double)rand() / RAND_MAX) - 250.0;
		obj.p.y = 500.0 * ((double)rand() / RAND_MAX) - 250.0;
		obj.p.z = 500.0 * ((double)rand() / RAND_MAX) - 250.0;
		u.addObject(obj);
	}
	
	update_universe_display();
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
	glutInitWindowSize(1000, 1000);
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
