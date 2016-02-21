#include "Universe.cuh"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include "timer.hpp"
#include <GL/glut.h>

Universe u;
float* vbuf;
double angle = 0.0, zoom = 1.0;
void display()
{	
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();
	float3 com = u.getCenterOfMass();
	
	glScalef(zoom, zoom, zoom);
	gluLookAt(com.x, com.y, com.z - 10000, com.x, com.y, com.z, 0, 1, 0);
	glRotated(angle, 1, 0, 0);
	
	
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, vbuf);
	glDrawArrays(GL_POINTS, 0, u.size());
	glDisableClientState(GL_VERTEX_ARRAY);
	
	
	glutSwapBuffers();
}

void reshape(int w, int h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (double)w / (double)h, 1.0, 100000);
	glMatrixMode(GL_MODELVIEW);
}

void copy_to_vbuf()
{
	for(int i = 0; i < u.size(); i++)
	{
		vbuf[3 * i + 0] = u.getObjects()[i].p.x;
		vbuf[3 * i + 1] = u.getObjects()[i].p.y;
		vbuf[3 * i + 2] = u.getObjects()[i].p.z;
	}
}

void keyboard(unsigned char c, int x, int y)
{
	if(c == 'w')
		angle += 1;
	else if(c == 's')
		angle -= 1;
	else if(c == '=')
		zoom *= 1.2;
	else if(c == '-')
		zoom /= 1.2;
	glutPostRedisplay();
}

void idle()
{
	u.timeStep();
	copy_to_vbuf();
	glutPostRedisplay();
}

int main(int argc, char** argv)
{
	cudaSetDevice(0);
	srand(time(0));
	
	u.getConfiguration().max_velocity = 1e9;
	u.getConfiguration().softening_factor = 0.01;
	u.getConfiguration().threshold_angle = 35;
	u.getConfiguration().time_step = 0.1;
	
	const int N = 200000;
	vbuf = new float[3 * N];
	for(int i = 0; i < N/2; i++)
	{
		Object o;
		o.m = rand() % 100 + 1;
		o.p.x = 500.0 * rand() / RAND_MAX - 250.0;
		o.p.y = 500.0 * rand() / RAND_MAX - 250.0;
		o.p.z = 500.0 * rand() / RAND_MAX - 250.0;
		u.addObject(o);
	}
	
	for(int i = 0; i < N/2; i++)
	{
		Object o;
		o.m = 1000 + rand() % 100 + 1;
		o.p.x = 500.0 * rand() / RAND_MAX + 1000;
		o.p.y = 500.0 * rand() / RAND_MAX + 1000;
		o.p.z = 500.0 * rand() / RAND_MAX + 1000;
		u.addObject(o);
	}
	copy_to_vbuf();
	
	glutInit(&argc, argv);
	glutInitWindowSize(1000, 1000);
	glutCreateWindow("N-Body");
	glutInitDisplayMode(GLUT_DOUBLE);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutIdleFunc(idle);
	
	glClearColor(0.0, 0.0, 0.0, 0.0);
	
	glutMainLoop();
	return 0;
}
