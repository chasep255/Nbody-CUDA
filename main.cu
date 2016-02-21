#include "Universe.cuh"
#include "timer.hpp"
#include <iostream>
#include <GL/glut.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <thread>

#define N 1000000

Universe u;
int2 window_dim;
float3 com;
double angle = 0.0;
double zoom = 1.0;
float* vertex_buffer = nullptr;
bool paused = true;

void display()
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, window_dim.x / window_dim.y, 50, 10000);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0, 0.0, 1000.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	glScaled(zoom, zoom, zoom);
	glRotated(angle, 1.0, 0.0, 0.0);
	glTranslatef(-com.x, -com.y, -com.z);
	
	glBegin(GL_POINTS);
	{
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3fv((GLfloat*)&com);
		
		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 0.0f);
	}
	glEnd();
	
	glColor3f(1.0f, 1.0f, 1.0f);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, vertex_buffer);
	glDrawArrays(GL_POINTS, 0, N);
	glDisableClientState(GL_VERTEX_ARRAY);
	
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
	
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	double m = 0.0;
	
	Object* objects = u.getObjects();
	
	#pragma omp parallel for reduction(+: x, y, z, m)
	for(int i = 0; i < N; i++)
	{
		x += (double)objects[i].p.x * objects[i].m;
		y += (double)objects[i].p.y * objects[i].m;
		z += (double)objects[i].p.z * objects[i].m;
		m += objects[i].m;
		reinterpret_cast<float3*>(vertex_buffer)[i] =  objects[i].p;
	}
	com.x = x / m;
	com.y = y / m;
	com.z = z / m;
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
	if(!paused)
	{
		u.timeStep();
		update_universe_display();
		glutPostRedisplay();
	}
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitWindowSize(1000, 1000);
	glEnable(GL_DOUBLE);
	
	glutCreateWindow("N-BODY");
	
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
			case 'w': angle += 1.0; redisplay = true; break;
			case 's': angle -= 1.0; redisplay = true; break;
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
