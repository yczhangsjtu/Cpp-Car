#include <iostream>
#include <vector>
#include <fstream>
#include "mathdraw.h"

using namespace std;

vector<float> points;

GLfloat fx, fy;

void compute()
{
	
}

void draw()
{
	mdrGrid(10.0);
	mdrAxes();
	mdrCoordinate(50.0);
	for(int i = 0; i+3 < points.size(); i+=2)
		mdrLine(points[i],points[i+1],points[i+2],points[i+3],0xff0000);
	if(points.size()>0)
		mdrLine(points[points.size()-2],points[points.size()-1],fx,fy,0xff0000);
}

void move()
{
}

void quit()
{
	ofstream fout("trace.txt");
	fout << points.size()/2 << endl;
	for(int i = 0; i+1 < points.size(); i+=2)
	{
		fout << points[i] << " " << points[i+1] << endl;
	}
	fout.close();
}

void click(int button, int state, int x, int y)
{
	if(button==GLUT_LEFT_BUTTON)
	{
		if(state==GLUT_DOWN && points.size()<256*2)
		{
			points.push_back(fx);
			points.push_back(fy);
		}
	}
	else if(button==GLUT_RIGHT_BUTTON)
	{
		if(!points.empty())
		{
			points.pop_back();
			points.pop_back();
		}
	}
	glutPostRedisplay();
}

void mousemove(int x, int y)
{
	fx = (GLfloat)x/ww*(realScaleRight-realScaleLeft) + realScaleLeft;
	fy = (GLfloat)(hh-y)/hh*(realScaleTop-realScaleBottom) + realScaleBottom;
	glutPostRedisplay();
}

int main(int argc, char *argv[])
{
	ifstream trace("trace.txt");
	if(trace)
	{
		GLint n;
		GLfloat x,y;
		trace >> n;
		while(!trace.eof())
		{
			trace >> x >> y;
			points.push_back(x);
			points.push_back(y);
		}
		trace.close();
	}
	mdrInit(&argc,argv,compute,draw,move,quit,-1000,1000,-500,500);
	glutMouseFunc(click);
	glutPassiveMotionFunc(mousemove);
	glutMainLoop();
}