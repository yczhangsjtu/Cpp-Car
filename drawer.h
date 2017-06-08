#ifndef _DRAWER_H
#define _DRAWER_H

#include <GL/freeglut.h>
#include <cmath>
#include <string>
#include <cstring>
#include <vector>
#include <sstream>
#include <fstream>
#include <ctime>

enum RotDir{RotX, RotY, RotZ};
enum ObjectType{Empty,Line,Square,Triangle,Ball,BallPart,Circle,
				Tetrahedron,Cylinder,BoldLine,Cone};

typedef struct materialStruct {
	GLfloat ambient[4];
	GLfloat diffuse[4];
	GLfloat specular[4];
	GLfloat shininess;
} materialStruct;

typedef struct Object {
	ObjectType type;
	GLfloat color[3];
	GLfloat points[4][3];
	GLfloat radius;
	GLfloat specular;
	GLfloat shininess;
} Object;

typedef struct ObjectCombine {
	std::vector<Object> objects;
} ObjectCombine;

typedef struct PlainObject {
	ObjectType type;
	GLfloat color[3];
	GLfloat points[4][2];
	GLfloat radius;
} PlainObject;

typedef struct Symbol {
	std::vector<PlainObject> objects;
} Symbol;

const int sbSize = 8192;

int ww = 1350, hh = 700;

int scaleLeft=-3;
int scaleRight=3;
int scaleBottom=-3;
int scaleTop=3;
int scaleNear=3;
int scaleFar=100;

GLfloat shininess = 42.0;
GLfloat specular = 0.4;

double eyetheta = 2.5;
double eyedist = 20;
double eyeheight = 5;
double maxeyeheight = 20;
double mineyeheight = -20;
const double rollSpeed = 0.05;

double lighttheta = 2.5;
double lightdist = 20;
double lightheight = 5;
double maxlightheight = 20;

int mousex;
int mousey;

int timer;
int timepassed;

GLfloat cursorx = 0;
GLfloat cursory = 0;
GLfloat cursorz = 0;
int fixCursor = 2;

////////////////////Routines////////////////////////////////////////
void updateLight()
{
	GLfloat light0_pos[4]={lightdist*cos(lighttheta),lightdist*sin(lighttheta),
						lightheight,1.0};
	GLfloat light0_spt[4]={-lightdist*cos(lighttheta),-lightdist*sin(lighttheta),
						-lightheight,1.0};
	
	glLightfv(GL_LIGHT0,GL_POSITION,light0_pos);
	glLightfv(GL_LIGHT0,GL_SPOT_DIRECTION,light0_spt);
	glLightf(GL_LIGHT0,GL_SPOT_CUTOFF,180);
	
}

void initRoutine(const char * title,int *argc,char *argv[])
{
	glutInit(argc,argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(ww,hh);
	glutInitWindowPosition(0,0);
	glutCreateWindow(title);
	
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glClearColor(0.0,0.0,0.0,1.0);
	
	updateLight();
	
	timer = std::clock();
	timepassed = 0;
}

void keepRatio()
{
	if(ww<hh)
		glFrustum(scaleLeft,scaleRight,scaleBottom*(GLfloat)hh/(GLfloat)ww,
				scaleTop*(GLfloat)hh/(GLfloat)ww,
				scaleNear,scaleFar);
	else
		glFrustum(scaleLeft*(GLfloat)ww/(GLfloat)hh,
				scaleRight*(GLfloat)ww/(GLfloat)hh,
				scaleBottom,scaleTop,scaleNear,scaleFar);
}

void keep2DRatio()
{
	GLfloat fw = ww;
	GLfloat fh = hh;
	if(ww>hh)
		gluOrtho2D(0,1000*fw/fh,0,1000);
	else
		gluOrtho2D(0,1000,0,1000*fh/fw);
}

void reshapeWithoutDeform(int w, int h)
{
	ww = w;
	hh = h;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
	keepRatio();
	
	glMatrixMode(GL_MODELVIEW);
	glViewport(0,0,w,h);
}

void displayRoutine()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyedist*cos(eyetheta),eyedist*sin(eyetheta),eyeheight,0,0,0,0,0,1);
}

void advanceTimer()
{
	timepassed += (clock() - timer)/100;
	timer = clock();
}

void trackMouse(int x, int y)
{
	mousex = x;
	mousey = y;
}

void drag(int x, int y)
{
	double dx=x-mousex;
	double dy=y-mousey;
	eyetheta -= dx/100.0;
	if(eyetheta > 2*M_PI) eyetheta -= 2*M_PI;
	if(eyetheta < 0) eyetheta += 2*M_PI;
	eyeheight += rollSpeed*dy;
	if(eyeheight > maxeyeheight) eyeheight = maxeyeheight;
	if(eyeheight < mineyeheight) eyeheight = mineyeheight;
}

Object newObject(ObjectType type, GLfloat R, GLfloat G, GLfloat B,
				 GLfloat specular, GLfloat shininess, GLfloat r = 0)
{
	Object object;
	object.type = type;
	object.color[0] = R;
	object.color[1] = G;
	object.color[2] = B;
	object.radius = r;
	object.specular = specular;
	object.shininess = shininess;
	return object;
}

PlainObject newPlainObject(ObjectType type, GLfloat R, GLfloat G, GLfloat B, GLfloat r = 0)
{
	PlainObject object;
	object.type = type;
	object.color[0] = R;
	object.color[1] = G;
	object.color[2] = B;
	object.radius = r;
	return object;
}

//////////////////////////////////////////////////////////////////////////
void useMaterial(GLfloat R,GLfloat G,GLfloat B,GLfloat specular,GLfloat shininess)
{
	materialStruct material = {
		{0.3*R, 0.3*G, 0.3*B, 1.0},
		{0.8*R, 0.8*G, 0.8*B, 1.0},
		{specular, specular, specular, 1.0},
		shininess
	};
	glMaterialfv(GL_FRONT,GL_AMBIENT, material.ambient);
	glMaterialfv(GL_FRONT,GL_DIFFUSE, material.diffuse);
	glMaterialfv(GL_FRONT,GL_SPECULAR,material.specular);
	glMaterialf(GL_FRONT,GL_SHININESS,material.shininess);
}

void crossProduct(GLfloat v1[],GLfloat v2[],GLfloat r[])
{
	r[0] = v1[1]*v2[2]-v1[2]*v2[1];
	r[1] = v1[2]*v2[0]-v1[0]*v2[2];
	r[2] = v1[0]*v2[1]-v1[1]*v2[0];
}
////////////////////////////////////////////////////////////////////////
///////////////////////////////3D Draw Functions////////////////////////
void drawLine(const GLfloat (*p)[3], const GLfloat *color,
			  GLfloat shininess = 42.0,GLfloat specular = 0.4)
{
	useMaterial(color[0],color[1],color[2],specular,shininess);
	glBegin(GL_LINES);
		glVertex3fv(p[0]);
		glVertex3fv(p[1]);
	glEnd();
}



void drawGrid()
{
	GLfloat cursor[] = {cursorx, cursory, cursorz};
	GLfloat points[4][3];
	for(int i = 0; i<4; i++)
		points[i][fixCursor] = cursor[fixCursor];
	for(GLfloat u = -50; u < 50; u++ )
	{
		if(fixCursor==0)
		{
			points[0][1] = -50;
			points[0][2] = u;
			points[1][1] = 50;
			points[1][2] = u;
			
			points[2][1] = u;
			points[2][2] = -50;
			points[3][1] = u;
			points[3][2] = 50;
		}
		else if(fixCursor==1)
		{
			points[0][0] = -50;
			points[0][2] = u;
			points[1][0] = 50;
			points[1][2] = u;
			
			points[2][0] = u;
			points[2][2] = -50;
			points[3][0] = u;
			points[3][2] = 50;
		}
		else if(fixCursor==2)
		{
			points[0][0] = -50;
			points[0][1] = u;
			points[1][0] = 50;
			points[1][1] = u;
			
			points[2][0] = u;
			points[2][1] = -50;
			points[3][0] = u;
			points[3][1] = 50;
		}
		GLfloat color[]={1,1,1};
		drawLine(points,color);
		drawLine(points+2,color);
	}
}

void drawSquare(const GLfloat (*p)[3], const GLfloat * color,
				GLfloat shininess = 42.0,GLfloat specular = 0.4,
				int id = 0, GLenum mode = GL_RENDER)
{
	if(mode==GL_SELECT) glLoadName(id);
	GLfloat v1[3]={p[1][0]-p[0][0],p[1][1]-p[0][1],p[1][2]-p[0][2]};
	GLfloat v2[3]={p[2][0]-p[1][0],p[2][1]-p[1][1],p[2][2]-p[1][2]};
	GLfloat normal[3];
	crossProduct(v1,v2,normal);
	GLfloat d = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
	for(int i=0; i<3; i++) normal[i]/=d;
	glNormal3fv(normal);
	useMaterial(color[0],color[1],color[2],specular,shininess);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	glBegin(GL_POLYGON);
		glVertex3fv(p[0]);
		glVertex3fv(p[1]);
		glVertex3fv(p[2]);
		glVertex3fv(p[3]);
	glEnd();
}




void drawTriangle(const GLfloat (*p)[3], const GLfloat * color,
				  GLfloat shininess = 42.0,GLfloat specular = 0.4,
				  int id = 0, GLenum mode = GL_RENDER)
{
	if(mode==GL_SELECT) glLoadName(id);
	GLfloat v1[3]={p[1][0]-p[0][0],p[1][1]-p[0][1],p[1][2]-p[0][2]};
	GLfloat v2[3]={p[2][0]-p[1][0],p[2][1]-p[1][1],p[2][2]-p[1][2]};
	GLfloat normal[3];
	crossProduct(v1,v2,normal);
	GLfloat d = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
	for(int i=0; i<3; i++) normal[i]/=d;
	glNormal3fv(normal);
	useMaterial(color[0],color[1],color[2],specular,shininess);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	glBegin(GL_POLYGON);
		glVertex3fv(p[0]);
		glVertex3fv(p[1]);
		glVertex3fv(p[2]);
	glEnd();
}

void drawBall(const GLfloat p[3], GLfloat r, const GLfloat * color,
			  GLfloat shininess = 42.0,GLfloat specular = 0.4,
			  int id = 0, GLenum mode = GL_RENDER)
{
	if(mode==GL_SELECT) glLoadName(id);
	GLUquadricObj * ball = gluNewQuadric();
	gluQuadricDrawStyle(ball,GLU_FILL);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslatef(p[0],p[1],p[2]);
	useMaterial(color[0],color[1],color[2],specular,shininess);
	//gluSphere(ball,r,(int)(r/(scaleRight-scaleLeft)*100),(int)(r/(scaleRight-scaleLeft)*100));
	gluSphere(ball,r,20,20);
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void drawLight()
{
	GLfloat p[3] = {lightdist*cos(lighttheta),lightdist*sin(lighttheta), lightheight};
	GLfloat color[]={1,1,1};
	GLfloat r = fabs(scaleRight-scaleLeft)/20;
	drawBall(p,r,color);
}

void drawCursor(GLfloat *color)
{
	GLfloat p[3] = {cursorx,cursory,cursorz};
	GLfloat pxoy[3] = {cursorx,cursory,0};
	GLfloat pyoz[3] = {0,cursory,cursorz};
	GLfloat pzox[3] = {cursorx,0,cursorz};
	GLfloat px[3] = {cursorx,0,0};
	GLfloat py[3] = {0,cursory,0};
	GLfloat pz[3] = {0,0,cursorz};
	GLfloat r = fabs(scaleRight-scaleLeft)/20;
	
	drawBall(p,r,color);
	drawBall(pxoy,r,color);
	drawBall(pyoz,r,color);
	drawBall(pzox,r,color);
	drawBall(px,r,color);
	drawBall(py,r,color);
	drawBall(pz,r,color);
}

void drawCircle(const GLfloat (*p)[3], GLfloat r, const GLfloat * color,
				  GLfloat shininess = 42.0,GLfloat specular = 0.4,
				  int id = 0, GLenum mode = GL_RENDER)
{
	GLfloat d[3] = {p[1][0]-p[0][0],
					p[1][1]-p[0][1],
					p[1][2]-p[0][2]};
	GLfloat dd = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
	GLfloat mu  = d[0]*d[0] + d[1]*d[1] + (d[2]-dd)*(d[2]-dd);
	GLfloat a11,a22,a33,a32,a31,a21;
	if(mu == 0)
	{
		a33 = 1; a22 = -1; a11 = 1;
		a32 = 0; a31 = 0; a21 = 0;
	}
	else
	{
		a22 = 1-2*d[1]*d[1]/mu;
		a11 = 1-2*d[0]*d[0]/mu;
		a32 = 2*d[1]*(dd-d[2])/mu;
		a31 = 2*d[0]*(dd-d[2])/mu;
		a21 = -2*d[0]*d[1]/mu;
	}
	int n = (int)(r/(scaleRight-scaleLeft)*1000);
	GLfloat dt = 2*M_PI/n;
	
	GLfloat normal[3];
	for(int i=0; i<3; i++) normal[i]=d[i]/dd;
	glNormal3fv(normal);
	useMaterial(color[0],color[1],color[2],specular,shininess);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	glBegin(GL_POLYGON);
	for(int i = 0; i < n; i++)
	{
		GLfloat t = i*dt;
		GLfloat st = sin(t);
		GLfloat ct = cos(t);
		
		GLfloat point[3] = 
		  {a11*r*ct+a21*r*st+p[0][0],
		   a21*r*ct+a22*r*st+p[0][1],
		   a31*r*ct+a32*r*st+p[0][2]};
		glVertex3fv(point);
	}
	glEnd();
}

void drawBallPart(const GLfloat (*p)[3], GLfloat r, const GLfloat * color,
				  GLfloat shininess = 42.0,GLfloat specular = 0.4,
				  int id = 0, GLenum mode = GL_RENDER)
{
	GLfloat dx = p[1][0]-p[0][0];
	GLfloat dy = p[1][1]-p[0][1];
	GLfloat dz = p[1][2]-p[0][2];
	GLfloat d2 = dx*dx+dy*dy+dz*dz;
	if(d2/r/r<0.0001)
	{
		drawCircle(p,r,color,shininess,specular,id,mode);
	}
	else
	{
		GLfloat l = (1-r*r/d2)/2;
		GLfloat br = sqrt(l*d2*l+r*r);
		GLfloat c[3] = {l*p[1][0]+(1-l)*p[0][0],
						l*p[1][1]+(1-l)*p[0][1],
						l*p[1][2]+(1-l)*p[0][2]};
		int n;
		if(l>0) n = (int)(br/(scaleRight-scaleLeft)*150);
		else    n = (int)(r/(scaleRight-scaleLeft)*150);
		GLfloat the = M_PI - acos(l/(1-l));
		GLfloat dlo = 2*M_PI/n;
		GLfloat dla = the/n;
		GLfloat d[3] = {p[1][0]-c[0],
						p[1][1]-c[1],
						p[1][2]-c[2]};
		GLfloat mu  = d[0]*d[0] + d[1]*d[1] + (d[2]-br)*(d[2]-br);
		GLfloat a11,a22,a33,a32,a31,a21;
		if(mu == 0)
		{
			a33 = 1; a22 = -1; a11 = 1;
			a32 = 0; a31 = 0; a21 = 0;
		}
		else
		{
			a33 = 1-2*(br-d[2])*(br-d[2])/mu;
			a22 = 1-2*d[1]*d[1]/mu;
			a11 = 1-2*d[0]*d[0]/mu;
			a32 = 2*d[1]*(br-d[2])/mu;
			a31 = 2*d[0]*(br-d[2])/mu;
			a21 = -2*d[0]*d[1]/mu;
		}
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
			{
				GLfloat t[] = {(i+1) * dla, i*dla};
				GLfloat st[] = {sin(t[0]),sin(t[1])};
				GLfloat ct[] = {cos(t[0]),cos(t[1])};
				GLfloat s[] = {j * dlo, (j+1)*dlo};
				GLfloat ss[] = {sin(s[0]),sin(s[1])};
				GLfloat cs[] = {cos(s[0]),cos(s[1])};
				
				GLfloat p[][3] = 
				 {{a11*br*st[0]*cs[0]+a21*br*st[0]*ss[0]+a31*br*ct[0]+c[0],
				   a21*br*st[0]*cs[0]+a22*br*st[0]*ss[0]+a32*br*ct[0]+c[1],
				   a31*br*st[0]*cs[0]+a32*br*st[0]*ss[0]+a33*br*ct[0]+c[2]},
				  {a11*br*st[1]*cs[0]+a21*br*st[1]*ss[0]+a31*br*ct[1]+c[0],
				   a21*br*st[1]*cs[0]+a22*br*st[1]*ss[0]+a32*br*ct[1]+c[1],
				   a31*br*st[1]*cs[0]+a32*br*st[1]*ss[0]+a33*br*ct[1]+c[2]},
				  {a11*br*st[1]*cs[1]+a21*br*st[1]*ss[1]+a31*br*ct[1]+c[0],
				   a21*br*st[1]*cs[1]+a22*br*st[1]*ss[1]+a32*br*ct[1]+c[1],
				   a31*br*st[1]*cs[1]+a32*br*st[1]*ss[1]+a33*br*ct[1]+c[2]},
				  {a11*br*st[0]*cs[1]+a21*br*st[0]*ss[1]+a31*br*ct[0]+c[0],
				   a21*br*st[0]*cs[1]+a22*br*st[0]*ss[1]+a32*br*ct[0]+c[1],
				   a31*br*st[0]*cs[1]+a32*br*st[0]*ss[1]+a33*br*ct[0]+c[2]}
				  };
				drawSquare(p,color,shininess,specular);
			}
	}
}

void drawCylinder(const GLfloat (*p)[3], GLfloat r, const GLfloat *color,
			  GLfloat shininess = 42.0,GLfloat specular = 0.4,
			  int id = 0, GLenum mode = GL_RENDER)
{
	GLfloat p1[2][3] = {{p[0][0],p[0][1],p[0][2]},{2*p[0][0]-p[1][0],2*p[0][1]-p[1][1],2*p[0][2]-p[1][2]}};
	GLfloat p2[2][3] = {{p[1][0],p[1][1],p[1][2]},{2*p[1][0]-p[0][0],2*p[1][1]-p[0][1],2*p[1][2]-p[0][2]}};
	drawCircle(p1,r,color,shininess,specular);
	drawCircle(p2,r,color,shininess,specular);
	
	GLfloat d[3] = {p[1][0]-p[0][0],
					p[1][1]-p[0][1],
					p[1][2]-p[0][2]};
	GLfloat l = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
	GLfloat mu  = d[0]*d[0] + d[1]*d[1] + (d[2]-l)*(d[2]-l);
	GLfloat a11,a22,a33,a32,a31,a21;
	if(mu == 0)
	{
		a33 = 1; a22 = -1; a11 = 1;
		a32 = 0; a31 = 0; a21 = 0;
	}
	else
	{
		a33 = 1-2*(l-d[2])*(l-d[2])/mu;
		a22 = 1-2*d[1]*d[1]/mu;
		a11 = 1-2*d[0]*d[0]/mu;
		a32 = 2*d[1]*(l-d[2])/mu;
		a31 = 2*d[0]*(l-d[2])/mu;
		a21 = -2*d[0]*d[1]/mu;
	}
			
	int n = (int)(r/(scaleRight-scaleLeft)*1000);
	GLfloat dt = 2*M_PI/n;
	
	for(int i = 0; i < n; i++)
	{
		GLfloat t[] = {i*dt,(i+1)*dt};
		GLfloat st[] = {sin(t[0]),sin(t[1])};
		GLfloat ct[] = {cos(t[0]),cos(t[1])};
		
		GLfloat normal[3];
		for(int i=0; i<3; i++) normal[i]=d[i]/l;
		
		GLfloat points[4][3]=
		{{a11*r*ct[0]+a21*r*st[0]+a31*l+p[0][0],
		  a21*r*ct[0]+a22*r*st[0]+a32*l+p[0][1],
		  a31*r*ct[0]+a32*r*st[0]+a33*l+p[0][2]},
		 {a11*r*ct[1]+a21*r*st[1]+a31*l+p[0][0],
		  a21*r*ct[1]+a22*r*st[1]+a32*l+p[0][1],
		  a31*r*ct[1]+a32*r*st[1]+a33*l+p[0][2]},
		 {a11*r*ct[1]+a21*r*st[1]+p[0][0],
		  a21*r*ct[1]+a22*r*st[1]+p[0][1],
		  a31*r*ct[1]+a32*r*st[1]+p[0][2]},
		 {a11*r*ct[0]+a21*r*st[0]+p[0][0],
		  a21*r*ct[0]+a22*r*st[0]+p[0][1],
		  a31*r*ct[0]+a32*r*st[0]+p[0][2]}};
		 drawSquare(points,color,shininess,specular,id,mode);
	}
}

void drawCone(const GLfloat (*p)[3], GLfloat r, const GLfloat * color,
			  GLfloat shininess = 42.0,GLfloat specular = 0.4,
			  int id = 0, GLenum mode = GL_RENDER)
{
	GLfloat p1[2][3] = {{p[0][0],p[0][1],p[0][2]},{2*p[0][0]-p[1][0],2*p[0][1]-p[1][1],2*p[0][2]-p[1][2]}};
	drawCircle(p1,r,color,shininess,specular);
	
	GLfloat d[3] = {p[1][0]-p[0][0],
					p[1][1]-p[0][1],
					p[1][2]-p[0][2]};
	GLfloat l = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
	GLfloat mu  = d[0]*d[0] + d[1]*d[1] + (d[2]-l)*(d[2]-l);
	GLfloat a11,a22,a33,a32,a31,a21;
	if(mu == 0)
	{
		a33 = 1; a22 = -1; a11 = 1;
		a32 = 0; a31 = 0; a21 = 0;
	}
	else
	{
		a33 = 1-2*(l-d[2])*(l-d[2])/mu;
		a22 = 1-2*d[1]*d[1]/mu;
		a11 = 1-2*d[0]*d[0]/mu;
		a32 = 2*d[1]*(l-d[2])/mu;
		a31 = 2*d[0]*(l-d[2])/mu;
		a21 = -2*d[0]*d[1]/mu;
	}
			
	int n = (int)(r/(scaleRight-scaleLeft)*1000);
	GLfloat dt = 2*M_PI/n;
	
	for(int i = 0; i < n; i++)
	{
		GLfloat t[] = {i*dt,(i+1)*dt};
		GLfloat st[] = {sin(t[0]),sin(t[1])};
		GLfloat ct[] = {cos(t[0]),cos(t[1])};
		
		GLfloat normal[3];
		for(int i=0; i<3; i++) normal[i]=d[i]/l;
		
		GLfloat points[3][3]=
		{{a31*l+p[0][0],
		  a32*l+p[0][1],
		  a33*l+p[0][2]},
		 {a11*r*ct[1]+a21*r*st[1]+p[0][0],
		  a21*r*ct[1]+a22*r*st[1]+p[0][1],
		  a31*r*ct[1]+a32*r*st[1]+p[0][2]},
		 {a11*r*ct[0]+a21*r*st[0]+p[0][0],
		  a21*r*ct[0]+a22*r*st[0]+p[0][1],
		  a31*r*ct[0]+a32*r*st[0]+p[0][2]}};
		 drawTriangle(points,color,shininess,specular,id,mode);
	}
}

void drawBoldLine(const GLfloat (*p)[3], const GLfloat *color,
			  GLfloat shininess = 42.0,GLfloat specular = 0.4,
			  int id = 0, GLenum mode = GL_RENDER)
{
	drawCylinder(p,0.1,color,shininess,specular,id,mode);
}

void drawAxes()
{
	GLfloat points[3][2][3];
	for(int i = 0; i<3; i++)
		for(int j = 0; j<2; j++)
			for(int k = 0; k<3; k++)
				points[i][j][k] = 0;
	points[0][0][0] = -50;
	points[0][1][0] = 50;
	points[1][0][1] = -50;
	points[1][1][1] = 50;
	points[2][0][2] = -50;
	points[2][1][2] = 50;
	{
	GLfloat color[]={1,0,0};
	drawCylinder(points[0],0.1,color);
	}
	{
	GLfloat color[]={0,1,0};
	drawCylinder(points[1],0.1,color);
	}
	{
	GLfloat color[]={0,0,1};
	drawCylinder(points[2],0.1,color);
	}
}

void drawTetrahedron(const GLfloat (*p)[3], const GLfloat * color,
				GLfloat shininess = 42.0,GLfloat specular = 0.4,
				int id = 0, GLenum mode = GL_RENDER)
{
	if(mode==GL_SELECT) glLoadName(id);
	//Draw first triangle
	{
		GLfloat v1[3]={p[1][0]-p[0][0],p[1][1]-p[0][1],p[1][2]-p[0][2]};
		GLfloat v2[3]={p[2][0]-p[1][0],p[2][1]-p[1][1],p[2][2]-p[1][2]};
		GLfloat normal[3];
		
		crossProduct(v1,v2,normal);
		GLfloat d = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
		for(int i=0; i<3; i++) normal[i]/=d;
		
		GLfloat orient[] = {p[3][0]-p[0][0],p[3][1]-p[0][1],p[3][2]-p[0][2]};
		GLfloat inner = orient[0]*normal[0] + orient[1]*normal[1] + orient[2]*normal[2];
		if(inner>0) for(int i=0; i<3; i++) normal[i]*=-1;
		
		glNormal3fv(normal);
		useMaterial(color[0],color[1],color[2],specular,shininess);
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		glBegin(GL_POLYGON);
			glVertex3fv(p[0]);
			glVertex3fv(p[1]);
			glVertex3fv(p[2]);
		glEnd();
	}
	//Draw second triangle
	{
		GLfloat v1[3]={p[2][0]-p[1][0],p[2][1]-p[1][1],p[2][2]-p[1][2]};
		GLfloat v2[3]={p[3][0]-p[2][0],p[3][1]-p[2][1],p[3][2]-p[2][2]};
		GLfloat normal[3];
		
		crossProduct(v1,v2,normal);
		GLfloat d = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
		for(int i=0; i<3; i++) normal[i]/=d;
		
		GLfloat orient[] = {p[0][0]-p[1][0],p[0][1]-p[1][1],p[0][2]-p[1][2]};
		GLfloat inner = orient[0]*normal[0] + orient[1]*normal[1] + orient[2]*normal[2];
		if(inner>0) for(int i=0; i<3; i++) normal[i]*=-1;
		
		glNormal3fv(normal);
		useMaterial(color[0],color[1],color[2],specular,shininess);
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		glBegin(GL_POLYGON);
			glVertex3fv(p[1]);
			glVertex3fv(p[2]);
			glVertex3fv(p[3]);
		glEnd();
	}
	//Draw third triangle
	{
		GLfloat v1[3]={p[3][0]-p[2][0],p[3][1]-p[2][1],p[3][2]-p[2][2]};
		GLfloat v2[3]={p[0][0]-p[3][0],p[0][1]-p[3][1],p[0][2]-p[3][2]};
		GLfloat normal[3];
		
		crossProduct(v1,v2,normal);
		GLfloat d = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
		for(int i=0; i<3; i++) normal[i]/=d;
		
		GLfloat orient[] = {p[1][0]-p[2][0],p[1][1]-p[2][1],p[1][2]-p[2][2]};
		GLfloat inner = orient[0]*normal[0] + orient[1]*normal[1] + orient[2]*normal[2];
		if(inner>0) for(int i=0; i<3; i++) normal[i]*=-1;
		
		glNormal3fv(normal);
		useMaterial(color[0],color[1],color[2],specular,shininess);
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		glBegin(GL_POLYGON);
			glVertex3fv(p[2]);
			glVertex3fv(p[3]);
			glVertex3fv(p[0]);
		glEnd();
	}
	//Draw fourth triangle
	{
		GLfloat v1[3]={p[0][0]-p[3][0],p[0][1]-p[3][1],p[0][2]-p[3][2]};
		GLfloat v2[3]={p[1][0]-p[0][0],p[1][1]-p[0][1],p[1][2]-p[0][2]};
		GLfloat normal[3];
		
		crossProduct(v1,v2,normal);
		GLfloat d = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
		for(int i=0; i<3; i++) normal[i]/=d;
		
		GLfloat orient[] = {p[2][0]-p[3][0],p[2][1]-p[3][1],p[2][2]-p[3][2]};
		GLfloat inner = orient[0]*normal[0] + orient[1]*normal[1] + orient[2]*normal[2];
		if(inner>0) for(int i=0; i<3; i++) normal[i]*=-1;
		
		glNormal3fv(normal);
		useMaterial(color[0],color[1],color[2],specular,shininess);
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		glBegin(GL_POLYGON);
			glVertex3fv(p[3]);
			glVertex3fv(p[0]);
			glVertex3fv(p[1]);
		glEnd();
	}
}

void drawFunction(double(*z)(double,double),
				  GLfloat minx, GLfloat miny, GLfloat maxx, GLfloat maxy,
				  int nxpieces,int nypieces,GLfloat *color,
				  int id = 0, GLenum mode = GL_RENDER)
{
	if(mode==GL_SELECT) glLoadName(id);
	GLfloat p[4][3];
	for(int i = 0; i<nxpieces; i++)
	{
		for(int j = 0; j<nypieces; j++)
		{
			GLfloat x1 = i*(maxx-minx)/nxpieces+minx;
			GLfloat x2 = (i+1)*(maxx-minx)/nxpieces+minx;
			GLfloat y1 = j*(maxy-miny)/nypieces+miny;
			GLfloat y2 = (j+1)*(maxy-miny)/nypieces+miny;
			p[0][0] = x1; p[0][1] = y1; p[0][2] = z(p[0][0],p[0][1]);
			p[1][0] = x2; p[1][1] = y1; p[1][2] = z(p[1][0],p[1][1]);
			p[2][0] = x2; p[2][1] = y2; p[2][2] = z(p[2][0],p[2][1]);
			p[3][0] = x1; p[3][1] = y2; p[3][2] = z(p[3][0],p[3][1]);
			drawSquare(p,color);
		}
	}
}
//////////////////////////////////////////////////////////////////////////////
void drawObject(const Object &object)
{
	if(object.type == Square)
	{
		drawSquare(object.points,object.color,object.shininess,object.specular);
	}
	else if(object.type == Line)
	{
		drawLine(object.points,object.color,object.shininess,object.specular);
	}
	else if(object.type == Triangle)
	{
		drawTriangle(object.points,object.color,object.shininess,object.specular);
	}
	else if(object.type == Ball)
	{
		drawBall(object.points[0],object.radius,object.color,object.shininess,object.specular);
	}
	else if(object.type == BallPart)
	{
		drawBallPart(object.points,object.radius,object.color,object.shininess,object.specular);
	}
	else if(object.type == Circle)
	{
		drawCircle(object.points,object.radius,object.color,object.shininess,object.specular);
	}
	else if(object.type == Tetrahedron)
	{
		drawTetrahedron(object.points,object.color,object.shininess,object.specular);
	}
	else if(object.type == Cylinder)
	{
		drawCylinder(object.points,object.radius,object.color,object.shininess,object.specular);
	}
	else if(object.type == BoldLine)
	{
		drawBoldLine(object.points,object.color,object.shininess,object.specular);
	}
	else if(object.type == Cone)
	{
		drawCone(object.points,object.radius,object.color,object.shininess,object.specular);
	}
}

int numOfPoints(ObjectType type)
{
	switch(type)
	{
	case Ball:
		return 1;
	case Line:
	case BoldLine:
	case BallPart:
	case Circle:
	case Cylinder:
	case Cone:
		return 2;
	case Triangle:
		return 3;
	case Square:
	case Tetrahedron:
		return 4;
	default:
		return 0;
	}
}

void setMatrixIdentity(GLfloat (*matrix)[3])
{
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			matrix[i][j] = (i==j?1:0);
}

void getMatrix(RotDir rotDir, GLfloat rot, GLfloat ratio, GLfloat (*matrix)[3])
{
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			matrix[i][j] = (i==j?ratio:0);
	switch(rotDir)
	{
	case RotX:
		matrix[1][1] = ratio*cos(rot); matrix[1][2] = ratio*sin(rot);
		matrix[2][1] =-ratio*sin(rot); matrix[2][2] = ratio*cos(rot);
		break;
	case RotY:
		matrix[2][2] = ratio*cos(rot); matrix[2][0] = ratio*sin(rot);
		matrix[0][2] =-ratio*sin(rot); matrix[0][0] = ratio*cos(rot);
		break;
	case RotZ:
		matrix[0][0] = ratio*cos(rot); matrix[0][1] = ratio*sin(rot);
		matrix[1][0] =-ratio*sin(rot); matrix[1][1] = ratio*cos(rot);
		break;
	}
}

void linearTransform(const GLfloat oldP[3], GLfloat newP[3], const GLfloat (*matrix)[3])
{
	for(int i = 0; i < 3; i++)
	{
		newP[i] = 0;
		for(int j = 0; j < 3; j++)
			newP[i] += matrix[i][j]*oldP[j];
	}
}

void movePoint(GLfloat position[3], const GLfloat move[3])
{
	for(int i = 0; i < 3; i++)
		position[i] += move[i];
}

void matrixMult(const GLfloat (*A)[3], const GLfloat (*B)[3], GLfloat (*R)[3])
{
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			R[i][j] = 0;
			for(int k = 0; k < 3; k++)
				R[i][j] += A[i][k] * B[k][j];
		}
	}
}

void drawObjectWithMatrix(const Object &object, GLfloat (*matrix)[3], GLfloat ratio, const GLfloat position[3])
{
	Object nobj = object;
	nobj.radius *= ratio;
	int npoint = numOfPoints(object.type);
	for(int i = 0; i < npoint; i++)
	{
		linearTransform(object.points[i],nobj.points[i],matrix);
		movePoint(nobj.points[i],position);
	}
	drawObject(nobj);
}

void drawPartialObject(const Object & object,int step)
{
	switch(object.type)
	{
	case Line:
	case BoldLine:
	case Ball:
		if(step==1)
			drawObject(object);
		break;
	case Triangle:
		if(step==1)
		{
			drawLine(object.points,object.color,object.shininess,object.specular);
		}
		else if(step==2)
		{
			drawObject(object);
		}
		break;
	case Square:
		if(step==1)
		{
			drawLine(object.points,object.color,object.shininess,object.specular);
		}
		else if(step==2)
		{
			drawTriangle(object.points,object.color,object.shininess,object.specular);
		}
		else if(step==3)
		{
			drawObject(object);
		}
		break;
	case BallPart:
		if(step==1)
			drawBall(object.points[0],object.radius,object.color,object.shininess,object.specular);
		else if(step==2)
			drawObject(object);
		 break;
	case Circle:
		if(step==2)
			drawObject(object);
		 break;
	case Tetrahedron:
		if(step==1)
		{
			drawLine(object.points,object.color,object.shininess,object.specular);
		}
		else if(step==2)
		{
			drawTriangle(object.points,object.color,object.shininess,object.specular);
		}
		else if(step==3)
		{
			drawObject(object);
		}
		break;
	case Cylinder:
	case Cone:
		if(step==1)
		{
			drawLine(object.points,object.color,object.shininess,object.specular);
		}
		else if(step==2)
		{
			drawObject(object);
		}
		break;
	}
}

bool intepret(const std::string &str, Object *object)
{
	std::istringstream is(str);
	std::string type;
	is >> type;
	if(type=="Square")
	{
		object->type = Square;
		is >> (object->color)[0] >> (object->color)[1] >> (object->color)[2]
		   >> (object->points)[0][0] >> (object->points)[0][1] >> (object->points)[0][2]
		   >> (object->points)[1][0] >> (object->points)[1][1] >> (object->points)[1][2]
		   >> (object->points)[2][0] >> (object->points)[2][1] >> (object->points)[2][2]
		   >> (object->points)[3][0] >> (object->points)[3][1] >> (object->points)[3][2];
		if(!is.eof()) is >> (object->specular) >> (object->shininess);
		else
		{
			object->specular = specular;
			object->shininess = shininess;
		}
	}
	else if(type=="Line")
	{
		object->type = Line;
		is >> (object->color)[0] >> (object->color)[1] >> (object->color)[2]
		   >> (object->points)[0][0] >> (object->points)[0][1] >> (object->points)[0][2]
		   >> (object->points)[1][0] >> (object->points)[1][1] >> (object->points)[1][2];
		if(!is.eof()) is >> (object->specular) >> (object->shininess);
		else
		{
			object->specular = specular;
			object->shininess = shininess;
		}
	}
	else if(type=="BoldLine")
	{
		object->type = BoldLine;
		is >> (object->color)[0] >> (object->color)[1] >> (object->color)[2]
		   >> (object->points)[0][0] >> (object->points)[0][1] >> (object->points)[0][2]
		   >> (object->points)[1][0] >> (object->points)[1][1] >> (object->points)[1][2];
		if(!is.eof()) is >> (object->specular) >> (object->shininess);
		else
		{
			object->specular = specular;
			object->shininess = shininess;
		}
	}
	else if(type=="Triangle")
	{
		object->type = Triangle;
		is >> (object->color)[0] >> (object->color)[1] >> (object->color)[2]
		   >> (object->points)[0][0] >> (object->points)[0][1] >> (object->points)[0][2]
		   >> (object->points)[1][0] >> (object->points)[1][1] >> (object->points)[1][2]
		   >> (object->points)[2][0] >> (object->points)[2][1] >> (object->points)[2][2];
		if(!is.eof()) is >> (object->specular) >> (object->shininess);
		else
		{
			object->specular = specular;
			object->shininess = shininess;
		}
	}
	else if(type=="Ball")
	{
		object->type = Ball;
		is >> (object->color)[0] >> (object->color)[1] >> (object->color)[2] >> (object->radius)
		   >> (object->points)[0][0] >> (object->points)[0][1] >> (object->points)[0][2];
		if(!is.eof()) is >> (object->specular) >> (object->shininess);
		else
		{
			object->specular = specular;
			object->shininess = shininess;
		}
	}
	else if(type=="BallPart")
	{
		object->type = BallPart;
		is >> (object->color)[0] >> (object->color)[1] >> (object->color)[2] >> (object->radius)
		   >> (object->points)[0][0] >> (object->points)[0][1] >> (object->points)[0][2]
		   >> (object->points)[1][0] >> (object->points)[1][1] >> (object->points)[1][2];
		if(!is.eof()) is >> (object->specular) >> (object->shininess);
		else
		{
			object->specular = specular;
			object->shininess = shininess;
		}
	}
	else if(type=="Circle")
	{
		object->type = Circle;
		is >> (object->color)[0] >> (object->color)[1] >> (object->color)[2] >> (object->radius)
		   >> (object->points)[0][0] >> (object->points)[0][1] >> (object->points)[0][2]
		   >> (object->points)[1][0] >> (object->points)[1][1] >> (object->points)[1][2];
		if(!is.eof()) is >> (object->specular) >> (object->shininess);
		else
		{
			object->specular = specular;
			object->shininess = shininess;
		}
	}
	else if(type=="Tetrahedron")
	{
		object->type = Tetrahedron;
		is >> (object->color)[0] >> (object->color)[1] >> (object->color)[2]
		   >> (object->points)[0][0] >> (object->points)[0][1] >> (object->points)[0][2]
		   >> (object->points)[1][0] >> (object->points)[1][1] >> (object->points)[1][2]
		   >> (object->points)[2][0] >> (object->points)[2][1] >> (object->points)[2][2]
		   >> (object->points)[3][0] >> (object->points)[3][1] >> (object->points)[3][2];
		if(!is.eof()) is >> (object->specular) >> (object->shininess);
		else
		{
			object->specular = specular;
			object->shininess = shininess;
		}
	}
	else if(type=="Cylinder")
	{
		object->type = Cylinder;
		is >> (object->color)[0] >> (object->color)[1] >> (object->color)[2] >> (object->radius)
		   >> (object->points)[0][0] >> (object->points)[0][1] >> (object->points)[0][2]
		   >> (object->points)[1][0] >> (object->points)[1][1] >> (object->points)[1][2];
		if(!is.eof()) is >> (object->specular) >> (object->shininess);
		else
		{
			object->specular = specular;
			object->shininess = shininess;
		}
	}
	else if(type=="Cone")
	{
		object->type = Cone;
		is >> (object->color)[0] >> (object->color)[1] >> (object->color)[2] >> (object->radius)
		   >> (object->points)[0][0] >> (object->points)[0][1] >> (object->points)[0][2]
		   >> (object->points)[1][0] >> (object->points)[1][1] >> (object->points)[1][2];
		if(!is.eof()) is >> (object->specular) >> (object->shininess);
		else
		{
			object->specular = specular;
			object->shininess = shininess;
		}
	}
	else
		return false;
	return true;
}

std::string express(const Object & object)
{
	std::ostringstream os;
	if(object.type==Line) os << "Line ";
	else if(object.type==BoldLine) os << "BoldLine ";
	else if(object.type==Square) os << "Square ";
	else if(object.type==Triangle) os << "Triangle ";
	else if(object.type==Ball) os << "Ball ";
	else if(object.type==BallPart) os << "BallPart ";
	else if(object.type==Circle) os << "Circle ";
	else if(object.type==Tetrahedron) os << "Tetrahedron ";
	else if(object.type==Cylinder) os << "Cylinder ";
	else if(object.type==Cone) os << "Cone ";
	else return os.str();
	
	for(int i = 0; i<3; i++)
		os << object.color[i] << " ";
	
	if(object.type==Ball || object.type==BallPart || object.type==Circle ||
	   object.type==Cylinder || object.type==Cone)
		os << object.radius << " ";
	
	for(int i = 0; i<3; i++)
		os << object.points[0][i] << " ";
		
	if(object.type==Line || object.type==BoldLine || object.type==Square ||
	   object.type==Tetrahedron || object.type==Triangle || object.type==BallPart ||
	   object.type==Circle || object.type==Cylinder || object.type==Cone)
		for(int i = 0; i<3; i++)
			os << object.points[1][i] << " ";
	
	if(object.type==Circle)
		for(int i = 0; i<3; i++)
			os << object.points[1][i] << " ";
	
	if(object.type==Square || object.type==Tetrahedron || object.type==Triangle)
		for(int i = 0; i<3; i++)
			os << object.points[2][i] << " ";
			
	if(object.type==Square || object.type==Tetrahedron)
		for(int i = 0; i<3; i++)
			os << object.points[3][i] << " ";
	
	os << object.specular << " " << object.shininess;
	return os.str();
}

void saveFile(const char* filename, const std::vector<Object> & objects)
{
	std::ofstream fout(filename);
	for(int i = 0; i < objects.size(); i++)
	{
		fout << express(objects[i]) << std::endl;
	}
	fout.close();
}

void openFile(const char* filename, std::vector<Object> & objects)
{
	std::ifstream fin(filename);
	if(!fin) return;
	
	objects.clear();
	while(!fin.eof())
	{
		objects.push_back(newObject(Empty,0,0,0,0,0));
		char str[256];
		fin.getline(str,256);
		if(!intepret(str,&objects.back()))
			objects.pop_back();
	}
	
	fin.close();
}

////////////////////////////////////////////////////////////////////////////////


//////////////////////////2D Draw Functions/////////////////////////////////
void prepareFor2D()
{
	glPushAttrib(GL_LIGHTING_BIT);
	glDisable(GL_LIGHTING);
	glDisable(GL_DEPTH_TEST);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	GLfloat fw = ww;
	GLfloat fh = hh;
	if(ww>hh)
		gluOrtho2D(0,1000*fw/fh,0,1000);
	else
		gluOrtho2D(0,1000,0,1000*fh/fw);
	
}

void goBackTo3D()
{
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glEnable(GL_DEPTH_TEST);
	glPopAttrib();
}

void printLetter(
	char letter,
	int x,
	int y,
	GLfloat letterColor[],
	GLfloat letterSize
)
{
	prepareFor2D();
	glTranslatef(x,y,0);
	GLfloat ratio = 
		letterSize/glutStrokeWidth(GLUT_STROKE_MONO_ROMAN,letter);
	glScalef(ratio,ratio,1);
	glColor3fv(letterColor);
	glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN,letter);
	goBackTo3D();
}

void printString(const char * str, GLfloat letterColor[], GLfloat letterSize, 
				 int x, int y)
{
	GLfloat curx, cury;
	curx = x; cury = y;
	for(int i = 0; str[i]!='\0'; i++, curx+=letterSize)
	{
		printLetter(str[i],curx,cury,letterColor,letterSize);
		if(str[i]=='\n')
		{
			cury-= 2*letterSize;
			curx = x-letterSize;
		}
	}
}

void printCoordinate()
{
	GLfloat color[3] = {1,1,1};
	std::ostringstream os;
	os << "x: " << cursorx << "\ny: " << cursory << "\nz: " << cursorz;
	printString(os.str().c_str(),color,20,50,950);
}

void printFloat(GLfloat number, GLfloat letterColor[], GLfloat letterSize, 
				 int x, int y)
{
	std::ostringstream os;
	os << number;
	printString(os.str().c_str(),letterColor,letterSize,x,y);
}

void printInt(int number, GLfloat letterColor[], GLfloat letterSize, 
				 int x, int y)
{
	std::ostringstream os;
	os << number;
	printString(os.str().c_str(),letterColor,letterSize,x,y);
}

void printLine(const GLfloat (*p)[2], const GLfloat * color,
			   int id = 0, GLenum mode = GL_RENDER)
{
	glColor3fv(color);
	glBegin(GL_LINES);
		glVertex2fv(p[0]);
		glVertex2fv(p[1]);
	glEnd();
}

void printTriangle(const GLfloat (*p)[2], const GLfloat * color,
				   int id = 0, GLenum mode = GL_RENDER)
{
	if(mode==GL_SELECT) glLoadName(id);
	glColor3fv(color);
	glBegin(GL_POLYGON);
		glVertex2fv(p[0]);
		glVertex2fv(p[1]);
		glVertex2fv(p[2]);
	glEnd();
}

void printCircle(const GLfloat p[2], GLfloat r, const GLfloat * color,
				 int id = 0, GLenum mode = GL_RENDER)
{
	if(mode==GL_SELECT) glLoadName(id);
	glColor3fv(color);
	int n = (int)(r/(scaleRight-scaleLeft)*100);
	GLfloat point[2];
	
	glBegin(GL_POLYGON);
	for(int i = 0; i<n; i++)
	{
		point[0] = p[0] + r*cos(2*M_PI/n*i);
		point[1] = p[1] + r*sin(2*M_PI/n*i);
		glVertex2fv(point);
	}
	glEnd();
}

void printSquare(const GLfloat (*p)[2], const GLfloat * color,
				 int id = 0, GLenum mode = GL_RENDER)
{
	if(mode==GL_SELECT) glLoadName(id);
	glColor3fv(color);
	glBegin(GL_POLYGON);
		glVertex2fv(p[0]);
		glVertex2fv(p[1]);
		glVertex2fv(p[2]);
		glVertex2fv(p[3]);
	glEnd();
}

void printObject(const PlainObject &object, int id=0, GLenum mode=GL_RENDER)
{
	if(object.type==Line)
	{
		printLine(object.points,object.color,id,mode);
	}
	else if(object.type==Triangle)
	{
		printTriangle(object.points,object.color,id,mode);
	}
	else if(object.type==Square)
	{
		printSquare(object.points,object.color,id,mode);
	}
	else if(object.type==Circle)
	{
		printCircle(object.points[0],object.radius,object.color,id,mode);
	}
}

void printSymbol(const Symbol &symbol, int id=0, GLenum mode=GL_RENDER)
{
	for(int i = 0; i<symbol.objects.size(); i++)
		printObject(symbol.objects[i],id,mode);
}

PlainObject newSquare(GLfloat x0,GLfloat y0,GLfloat x1,GLfloat y1,
					  GLfloat x2,GLfloat y2,GLfloat x3,GLfloat y3,
					  GLfloat R, GLfloat G, GLfloat B)
{
	PlainObject obj = newPlainObject(Square,R,G,B);
	obj.points[0][0] = x0; obj.points[0][1] = y0;
	obj.points[1][0] = x1; obj.points[1][1] = y1;
	obj.points[2][0] = x2; obj.points[2][1] = y2;
	obj.points[3][0] = x3; obj.points[3][1] = y3;
	return obj;
}

PlainObject newTriangle(GLfloat x0,GLfloat y0,GLfloat x1,GLfloat y1,
					    GLfloat x2,GLfloat y2,GLfloat R, GLfloat G, GLfloat B)
{
	PlainObject obj = newPlainObject(Triangle,R,G,B);
	obj.points[0][0] = x0; obj.points[0][1] = y0;
	obj.points[1][0] = x1; obj.points[1][1] = y1;
	obj.points[2][0] = x2; obj.points[2][1] = y2;
	return obj;
}

PlainObject newLine(GLfloat x0,GLfloat y0,GLfloat x1,GLfloat y1,GLfloat R, GLfloat G, GLfloat B)
{
	PlainObject obj = newPlainObject(Line,R,G,B);
	obj.points[0][0] = x0; obj.points[0][1] = y0;
	obj.points[1][0] = x1; obj.points[1][1] = y1;
	return obj;
}

PlainObject newCircle(GLfloat x0,GLfloat y0,GLfloat R, GLfloat G, GLfloat B,GLfloat r)
{
	PlainObject obj = newPlainObject(Circle,R,G,B,r);
	obj.points[0][0] = x0; obj.points[0][1] = y0;
	return obj;
}

Symbol newSymbol()
{
	Symbol symbol;
	return symbol;
}

bool symbolClicked(const Symbol &symbol)
{
	prepareFor2D();
	
	GLint hits;
	GLuint selectBuf[sbSize];
	GLint viewport[4];
	
	glGetIntegerv(GL_VIEWPORT,viewport);
	glSelectBuffer(sbSize,selectBuf);
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);
	
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	gluPickMatrix ((GLdouble) mousex, (GLdouble)(viewport[3]- mousey), 
		0.1, 0.1, viewport);
	keep2DRatio();
	
	printSymbol(symbol,0,GL_SELECT);
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	
	hits=glRenderMode(GL_RENDER);
	
	goBackTo3D();
	return hits > 0;
}
//////////////////////////////////////////////////////////////////////////


//////////////////////Other functions///////////////////////////////////////////
void locateMouse(GLfloat *x, GLfloat *y, GLfloat *z,
				 GLfloat u0, GLfloat v0, GLfloat u1, GLfloat v1,
				 int index)
{
	if(fabs(u1-u0) < 0.01 && fabs(v1-v0) < 0.01)
	{
		if(index==0)
		{
			*y = (u0+u1)/2;
			*z = (v0+v1)/2;
		}
		else if(index==1)
		{
			*x = (u0+u1)/2;
			*z = (v0+v1)/2;
		}
		else
		{
			*x = (u0+u1)/2;
			*y = (v0+v1)/2;
		}
		return;
	}
	GLint hits;
	GLuint selectBuf[sbSize];
	GLint viewport[4];
	
	glGetIntegerv(GL_VIEWPORT,viewport);
	glSelectBuffer(sbSize,selectBuf);
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);
	
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	gluPickMatrix ((GLdouble) mousex, (GLdouble)(viewport[3]- mousey), 
		0.1, 0.1, viewport);
		
	if(ww<hh)
		glFrustum(scaleLeft,scaleRight,scaleBottom*(GLfloat)hh/(GLfloat)ww,
				scaleTop*(GLfloat)hh/(GLfloat)ww,
				scaleNear,scaleFar);
	else
		glFrustum(scaleLeft*(GLfloat)ww/(GLfloat)hh,scaleRight*(GLfloat)ww/(GLfloat)hh,
				scaleBottom,scaleTop,scaleNear,scaleFar);glMatrixMode(GL_PROJECTION);
				
	GLfloat p[4][3];
	if(index==0)
	{
		p[0][0] = *x; p[1][0] = *x; p[2][0] = *x; p[3][0] = *x;
		p[0][1] = u0; p[1][1] = u0; p[2][1] = u1; p[3][1] = u1;
		p[0][2] = v0; p[1][2] = v1; p[2][2] = v1; p[3][2] = v0;
	}
	else if(index==1)
	{
		p[0][0] = u0; p[1][0] = u0; p[2][0] = u1; p[3][0] = u1;
		p[0][1] = *y; p[1][1] = *y; p[2][1] = *y; p[3][1] = *y;
		p[0][2] = v0; p[1][2] = v1; p[2][2] = v1; p[3][2] = v0;
	}
	else
	{
		p[0][0] = u0; p[1][0] = u0; p[2][0] = u1; p[3][0] = u1;
		p[0][1] = v0; p[1][1] = v1; p[2][1] = v1; p[3][1] = v0;
		p[0][2] = *z; p[1][2] = *z; p[2][2] = *z; p[3][2] = *z;
	}
	GLfloat color[3] = {1,1,1};
	
	drawSquare(p,color,0,0,0,GL_SELECT);
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	
	hits=glRenderMode(GL_RENDER);
	
	if(hits>0)
	{
		GLfloat um = (u0+u1)/2;
		GLfloat vm = (v0+v1)/2;
		locateMouse(x,y,z,u0,v0,um,vm,index);
		locateMouse(x,y,z,um,v0,u1,vm,index);
		locateMouse(x,y,z,u0,vm,um,v1,index);
		locateMouse(x,y,z,um,vm,u1,v1,index);
	}
}

int requiredNumberOfSteps(ObjectType type)
{
	switch(type)
	{
	case Line:
	case BoldLine:
	case Ball:
		return 2;
	case Triangle:
	case BallPart:
	case Circle:
	case Cylinder:
	case Cone:
		return 3;
	case Square:
	case Tetrahedron:
		return 4;
	default:
		return 0;
	}
}

void saveObjects(const char* filename, const std::vector<Object> & objects)
{
	std::ofstream fout(filename,std::ios::binary);
	int nbyte = sizeof (Object);
	unsigned char *buffer = new unsigned char[nbyte];
	
	for(int i = 0; i < objects.size(); i++)
	{
		std::memcpy(buffer,&objects[i],nbyte);
		for(int j = 0; j < nbyte; j++)
			fout << buffer[j];
	}
	
	delete [] buffer;
	fout.close();
}

void loadObjects(const char* filename, std::vector<Object> & objects)
{
	std::ifstream fin(filename,std::ios::binary);
	if(!fin) return;
	
	int nbyte = sizeof (Object);
	char *buffer = new char[nbyte];
	
	objects.clear();
	while(!fin.eof())
	{
		int i = 0;
		for(; i < nbyte; i++)
		{
			if(fin.eof()) break;
			fin.get(buffer[i]);
		}
		if(i!=nbyte) break;
		objects.push_back(newObject(Empty,0,0,0,0,0));
		std::memcpy(&objects.back(),buffer,nbyte);
	}
	
	delete [] buffer;
	fin.close();
}

#endif
