#include "common.h"
#include "spherical.h"
#include <cmath>
#include <GL/glut.h>

using namespace std;

double Time = 0.0, timestep;
int nstep, nprint, nenergy;
Vector3d G(0.0,0.0,0.0);
double lx, ly, lz, x_0, y_0, z_0, radii_max;
ParticleList particle(N);
unsigned int no_of_particles; // number of active particles
double dvx=0.1,vx,abs_vx,eps;
double sx,sy,sxy,sxq,sqx,a_val,b_val,GT;
int cnt;

ofstream fphase("phase.dat"), fenergy("energy.dat");

const double a1 = 0.92388, a2 = 0.382683, b = 0.707107;
const double cos0  = 1.0, sin0  = 0.0;
const double cos1  =  a1, sin1  =  a2;
const double cos2  =   b, sin2  =   b;
const double cos3  =  a2, sin3  =  a1;
const double cos4  = 0.0, sin4  = 1.0;
const double cos5  = -a2, sin5  =  a1;
const double cos6  =  -b, sin6  =   b;
const double cos7  = -a1, sin7  =  a2;
const double cos8  =-1.0, sin8  = 0.0;
const double cos9  = -a1, sin9  = -a2;
const double cos10 =  -b, sin10 =  -b;
const double cos11 = -a2, sin11 = -a1;
const double cos12 = 0.0, sin12 =-1.0;
const double cos13 =  a2, sin13 = -a1;
const double cos14 =   b, sin14 =  -b;
const double cos15 =  a1, sin15 = -a2;

void DrawCircle(GLfloat x, GLfloat y, GLfloat r)
{
        // Drawing a circle with OpenGL functions
        glBegin(GL_TRIANGLE_FAN);
                glVertex2f(x+r*cos0 ,y+r*sin0 );
                glVertex2f(x+r*cos1 ,y+r*sin1 );
                glVertex2f(x+r*cos2 ,y+r*sin2 );
                glVertex2f(x+r*cos3 ,y+r*sin3 );
                glVertex2f(x+r*cos4 ,y+r*sin4 );
                glVertex2f(x+r*cos5 ,y+r*sin5 );
                glVertex2f(x+r*cos6 ,y+r*sin6 );
                glVertex2f(x+r*cos7 ,y+r*sin7 );
                glVertex2f(x+r*cos8 ,y+r*sin8 );
                glVertex2f(x+r*cos9 ,y+r*sin9 );
                glVertex2f(x+r*cos10,y+r*sin10);
                glVertex2f(x+r*cos11,y+r*sin11);
                glVertex2f(x+r*cos12,y+r*sin12);
                glVertex2f(x+r*cos13,y+r*sin13);
                glVertex2f(x+r*cos14,y+r*sin14);
                glVertex2f(x+r*cos15,y+r*sin15);
        glEnd();

}

void SetupRC(void)
{
        glClearColor(0.0,0.0,0.0,1.0);
        glDisable(GL_DEPTH_TEST);
}

void RenderScene(void)
{
        GLfloat x,y,r;

        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f,0.0f,0.0f);
        for (unsigned int i=0; i<N_inactive; i++)
        {
                x = particle[i].x();
                y = particle[i].y();
                r = particle[i].r();
                DrawCircle(x,y,r);
        }
	glColor3f(1.0f,1.0f,1.0f);
        for (unsigned int i=N_inactive; i<particle.size(); i++)
        {
		x = particle[i].x();
               	y = particle[i].y();
	        r = particle[i].r();
                DrawCircle(x,y,r);
        }


        glutSwapBuffers();
}

void ChangeSize(GLsizei w, GLsizei h)
{
        GLfloat aspectRatio;
	double xmin, xmax, ymin, ymax;

	xmin = x_0; xmax = lx+x_0;
	ymin = y_0; ymax = ly+y_0;

        if (h==0)
                h=1;

        glViewport(0,0,w,h);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        aspectRatio=(GLfloat)w/(GLfloat)h;

        if (w<=h)
              gluOrtho2D(xmin,xmax,ymin/aspectRatio,ymax/aspectRatio);
        else
              gluOrtho2D(xmin*aspectRatio,xmax*aspectRatio,ymin,ymax);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

}

void TimerFunction(int value)
{
//	cout << Time << endl;
	step();
/*
	if ((particle[0].x()<100.0)&&(particle[0].vx()<0.0))
	{
		abs_vx = fabs(particle[0].vx());
		eps = abs_vx/vx;
		fphase << 2*vx <<" "<< eps <<" "<< 1.0-eps <<" "<< log(1.0-eps) <<" "<< log(2*vx) << endl;
		cnt++;
		sx += log(2*vx); sy += log(1.0-eps); sxy += log(2*vx)*log(1.0-eps); sxq = log(2*vx)*log(2*vx);
		sqx = sx*sx;
		a_val = (sx*sy-cnt*sxy)/(sqx-cnt*sxq);
		b_val = (sy-a_val*sx)/cnt;
		cout << 2*vx <<" "<< eps <<" "<<" "<< a_val <<" "<< b_val << endl;
		vx += dvx;
		particle[0].set_x(double(100));
		particle[0].set_vx(vx);
		particle[1].set_x(double(135));
		particle[1].set_vx(-vx);			
	}
*/
	GT = 0.0;
	for (unsigned int i=0; i<particle.size(); i++){
		GT += particle[i].kinetic_energy();
	}
	GT /= particle.size();

	fenergy << Time <<" "<< GT <<" "<< log(Time) <<" "<< log(GT) << endl;
			
        glutPostRedisplay();
        glutTimerFunc(1,TimerFunction,1);
}


int main(int argc, char* argv[])
{
	glutInit(&argc,argv);
        glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
        glutInitWindowSize(900,900);
        glutCreateWindow("Linkcell algorithm");
        glutDisplayFunc(RenderScene);
        glutReshapeFunc(ChangeSize);

        SetupRC();
        init_system();
	init_algorithm();

        glutTimerFunc(1,TimerFunction,1);
        glutMainLoop();



	return 0;
}
