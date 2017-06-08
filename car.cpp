#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>
#include <GL/freeglut.h>
#include "drawer.h"
#include "lib3ds.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb-master/stb_image.h"

using namespace std;

const char filepathcar[] = "car.3ds";
const char filepathtree[] = "tree.3ds";
Lib3dsFile * carfile, * treefile;

GLfloat speed = 0;

GLfloat pos[3] = {0,0,0};
GLfloat rot = 0.0;

GLfloat eyerot = rot;

GLfloat frictionOnRoad = 0.01;
GLfloat frictionOffRoad = 0.03;
GLfloat acceleration = 0.05;
GLfloat brake = 0.15;
GLfloat turnRate = 0.030;
GLfloat roadWidth = 35;
bool drawTrees = true;

bool accelerating = false;
bool turningLeft = false;
bool turningRight = false;
bool braking = false;

bool crash = false;
bool onroad = false;
int framePerSecond = 24;

const GLfloat southmost = -600;
const GLfloat northmost = 600;
const GLfloat eastmost = 1100;
const GLfloat westmost = -1100;

const int maxnroads = 256;

GLfloat roads[maxnroads][2] =
{{-524.444,-360},{582.222,-360},{746.667,-335.714},{834.074,-244.286},
    {878.518,-105.714},{877.037,147.143},{829.63,258.571},{736.296,341.429},
    {591.111,380},{428.148,382.857},{266.667,332.857},{200,237.143},
    {188.148,140},{143.704,75.7143},{63.7036,35.7143},{-84.4445,34.2857},
    {-271.111,35.7143},{-434.074,10},{-571.852,-64.2857},{-629.63,-160},
    {-641.481,-272.857},{-595.556,-344.286},{-524.444,-360}
};

const int roadTexSize = 128;
GLubyte roadtex[roadTexSize][roadTexSize][3];

const int groundTextureSize = 512;
GLubyte groundTexture[groundTextureSize][groundTextureSize][3];

int nroads = 22;

GLfloat (*road)[2] = roads;


static void
load_model(void)
{
    carfile=lib3ds_file_load(filepathcar);
	treefile=lib3ds_file_load(filepathtree);
	treefile->meshes->next = NULL;
}

void
render_mesh(Lib3dsMesh *mesh, int skip = 0)
{
    ASSERT(carfile);
	if(!mesh) return;
	if (!mesh->user.d)
	{
		// Arrange a free integer as the handle of this list
		mesh->user.d=glGenLists(1);
		// Assign the list to this handle
		glNewList(mesh->user.d, GL_COMPILE);

		{
			// Used to store the normal vectors of all faces in this mesh
			Lib3dsVector *normalL =
				(float (*)[3])malloc(3*sizeof(Lib3dsVector)*mesh->faces);
			Lib3dsMaterial *oldmat = (Lib3dsMaterial *) (-1);
			// Multiply the current matrix by the inverse of the matrix of this
			// mesh
			{
				Lib3dsMatrix M;
				lib3ds_matrix_copy(M, mesh->matrix);
				lib3ds_matrix_inv(M);
				glMultMatrixf(&M[0][0]);
			}
			// Get the normal of the current mesh
			lib3ds_mesh_calculate_normals(mesh, normalL);

			for (int p=0; p<mesh->faces; p+=(skip+1))
			{
				// Go through each face in the face list of the current mesh
				Lib3dsFace *f=&mesh->faceL[p];
				// Get the meterial of this face
				Lib3dsMaterial *mat=0;
				// material is a string, and if it's not empty, read the details
				// of this material from the file
				if (f->material[0])
					mat=lib3ds_file_material_by_name(carfile, f->material);
				// If mat != previous material(which means the OpenGL material needs
				// reset)
				if ( mat != oldmat && mat )
				{
					//glDisable(GL_CULL_FACE);
					for(int i = 0; i < 3; i++)
						mat->specular[i] *= 0.6;
					// Set the material of the OpenGL
					glMaterialfv(GL_FRONT, GL_AMBIENT, mat->ambient);
					glMaterialfv(GL_FRONT, GL_DIFFUSE, mat->diffuse);
					glMaterialfv(GL_FRONT, GL_SPECULAR, mat->specular);
					glMaterialf(GL_FRONT, GL_SHININESS, pow(2, 10.0*mat->shininess));
				}
				oldmat = mat;
				// Drawing
				{
					int i;
					glBegin(GL_TRIANGLES);
					glNormal3fv(f->normal);
					for (i=0; i<3; ++i)
					{
						glNormal3fv(normalL[3*p+i]);
						glVertex3fv(mesh->pointL[f->points[i]].pos);
					}
					glEnd();
				}
			}

			free(normalL);
		}

		glEndList();
	}
	if (mesh->user.d)
	{
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glCallList(mesh->user.d);
		glMatrixMode(GL_MODELVIEW);
		glPopMatrix();
	}
    if(mesh->next!=NULL) render_mesh(mesh->next,skip);
}


void initCar()
{
    pos[0] = road[0][0];
    pos[1] = road[0][1];
    pos[2] = 1;
    rot = 0;
}

void initGroundTexture()
{
    for (int i = 0; i < groundTextureSize; i++)
        for (int j = 0; j < groundTextureSize; j++)
        {
            groundTexture[i][j][0] = 0;
            groundTexture[i][j][1] = (GLfloat)rand()/RAND_MAX*127+128;
            groundTexture[i][j][2] = 0;
        }
    glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
}

void initRoadTexture()
{
	
	int x,y,n;
	unsigned char *data = stbi_load("road.bmp", &x, &y, &n, 0);
	if(data==NULL) return;
    for (int i = 0; i < roadTexSize; i++)
        for (int j = 0; j < roadTexSize; j++)
        {
            roadtex[i][j][0] = data[(j*y/roadTexSize*x+i*x/roadTexSize)*n];
            roadtex[i][j][1] = data[(j*y/roadTexSize*x+i*x/roadTexSize)*n+1];
            roadtex[i][j][2] = data[(j*y/roadTexSize*x+i*x/roadTexSize)*n+2];
        }
	stbi_image_free(data);
}

void initRoad()
{
	ifstream roadfile("road.txt");
	if(roadfile)
	{
		int n = 0;
		roadfile >> n;
		if(n > 0 && n < maxnroads)
		{
			nroads = n-1;
			for(int i = 0; i < n; i++)
				roadfile >> roads[i][0] >> roads[i][1];
		}
		roadfile.close();
	}
}

void init(const char * title,int *argc,char *argv[])
{
    initRoutine(title,argc,argv);
    mineyeheight = 7;
    maxeyeheight = 7;
    eyeheight = 7;
    eyedist = 10;
    glClearColor(0,0.5,0.5,0.5);
	initRoad();
    initCar();
	load_model();
    initGroundTexture();
    initRoadTexture();
	ifstream setting("setting.txt");
	if(setting)
	{
		char rule[128];
		char equal;
		while(!setting.eof())
		{
			setting >> rule;
			if(strncmp(rule,"Acceleration",128)==0)
			{
				setting >> equal;
				if(equal=='=')
					setting >> acceleration;
			}
			if(strncmp(rule,"Brake",128)==0)
			{
				setting >> equal;
				if(equal=='=')
					setting >> brake;
			}
			if(strncmp(rule,"TurnRate",128)==0)
			{
				setting >> equal;
				if(equal=='=')
					setting >> turnRate;
			}
			if(strncmp(rule,"RoadWidth",128)==0)
			{
				setting >> equal;
				if(equal=='=')
					setting >> roadWidth;
			}
			if(strncmp(rule,"FrictionOnRoad",128)==0)
			{
				setting >> equal;
				if(equal=='=')
					setting >> frictionOnRoad;
			}
			if(strncmp(rule,"FrictionOffRoad",128)==0)
			{
				setting >> equal;
				if(equal=='=')
					setting >> frictionOffRoad;
			}
			if(strncmp(rule,"DrawTrees",128)==0)
			{
				setting >> equal;
				if(equal=='=')
				{
					setting >> rule;
					if(strncmp(rule,"True",128)==0)
						drawTrees = true;
					else if(strncmp(rule,"False",128)==0)
						drawTrees = false;
				}
			}
			
		}
		setting.close();
	}
}

void drawBackground()
{
    glEnable(GL_TEXTURE_2D);
    GLfloat ground[4][3] = {{westmost,southmost,0},{eastmost,southmost,0},
        {eastmost,northmost,0},{westmost,northmost,0}
    };
    glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,groundTextureSize,groundTextureSize,
                 0,GL_RGB,GL_UNSIGNED_BYTE,groundTexture);
	useMaterial(0,0,0,0.6,0.5);
	glNormal3f(0,0,1);
    glColor3f(1,1,1);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0,0.0);
    glVertex3fv(ground[0]);
    glTexCoord2f(0.0,1.0);
    glVertex3fv(ground[1]);
    glTexCoord2f(1.0,1.0);
    glVertex3fv(ground[2]);
    glTexCoord2f(1.0,0.0);
    glVertex3fv(ground[3]);
    glEnd();
    glDisable(GL_TEXTURE_2D);
}

void drawTree(GLfloat x0, GLfloat y0, GLfloat h, GLfloat w, GLfloat size)
{
	if(fabs(x0-pos[0])>100 || fabs(y0-pos[1])>100) return;
	
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	
	useMaterial(0.0,1.0,0.0,0.3,0.3);
	glTranslatef(x0,y0,0);
	glScalef(0.1,0.1,0.025);
	if(treefile)
	{
		render_mesh(treefile->meshes,2);
	}
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	GLfloat p[][3] = {{x0,y0,0},{x0,y0,h}};
	GLfloat c[3] = {0.8,0.5,0.3};
	drawCylinder(p,w,c);
    //Judge if crash the tree
    GLfloat dx = x0 - pos[0];
    GLfloat dy = y0 - pos[1];
    if (dx*dx+dy*dy < (w+2)*(w+2))
        crash = true;
}

void drawTreeLine(GLfloat x0, GLfloat y0, GLfloat x1, GLfloat y1, GLfloat distance)
{
    GLfloat h = 8, size = 2, w = 0.5;
    int n = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))/distance;
    GLfloat dx = (x1-x0)/n;
    GLfloat dy = (y1-y0)/n;
    for (int i = 0; i < n; i++)
        drawTree(x0+i*dx,y0+i*dy,h,w,size);
}

void drawRoadPiece(GLfloat x0, GLfloat y0, GLfloat x1, GLfloat y1,
                   GLfloat width, GLfloat treeInterval)
{
    glEnable(GL_TEXTURE_2D);
    glTexImage2D
		(GL_TEXTURE_2D,0,GL_RGB,128,128,0,GL_RGB,GL_UNSIGNED_BYTE,roadtex);
    GLfloat dx = y1 - y0;
    GLfloat dy = x0 - x1;
    GLfloat l  = sqrt(dx*dx+dy*dy);
    dx = dx/l*width/2;
    dy = dy/l*width/2;
    GLfloat p[][3] =
		{{x0-dx+dy,y0-dy-dx,0.05},{x0+dx+dy,y0+dy-dx,0.05},
		 {x1+dx-dy,y1+dy+dx,0.1},{x1-dx-dy,y1-dy+dx,0.1}};
	useMaterial(1,1,1,1,1);
	glNormal3f(0,0,1);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0,0.0);
    glVertex3fv(p[0]);
    glTexCoord2f(0.0,1.0);
    glVertex3fv(p[1]);
    glTexCoord2f(1.0,1.0);
    glVertex3fv(p[2]);
    glTexCoord2f(1.0,0.0);
    glVertex3fv(p[3]);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    if (treeInterval > 5 && drawTrees)
    {
        drawTreeLine(x0-dx,y0-dy,x1-dx,y1-dy,treeInterval);
        drawTreeLine(x0+dx,y0+dy,x1+dx,y1+dy,treeInterval);
    }
	GLfloat len = sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
	GLfloat dl = (y1-y0)*pos[0]-(x1-x0)*pos[1]-(y1-y0)*x0+(x1-x0)*y0;
	dl = fabs(dl);
	GLfloat l0 = sqrt((pos[0]-x0)*(pos[0]-x0)+(pos[1]-y0)*(pos[1]-y0));
	GLfloat l1 = sqrt((pos[0]-x1)*(pos[0]-x1)+(pos[1]-y1)*(pos[1]-y1));
	if(dl < width*len && l0 < len+width && l1<len+width)
		onroad = true;
}

void drawPath(GLfloat (*p)[2], int n, GLfloat width, GLfloat treeInterval)
{
	onroad = false;
    for (int i = 0; i < n; i++)
        drawRoadPiece(p[i][0],p[i][1],p[i+1][0],p[i+1][1],width,treeInterval);
}

void drawRoad()
{
    GLfloat w = roadWidth;
    drawPath(road,nroads,w,40);
}

GLfloat mapx(GLfloat x)
{
    return 250/(eastmost-westmost)*(x-westmost);
}

GLfloat mapy(GLfloat y)
{
    return 120/(northmost-southmost)*(y-southmost) + 800;
}

void drawMap(GLfloat (*p)[2], int n)
{
    prepareFor2D();
    glColor3f(0,0,0);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < n; i++)
        glVertex2f(mapx(p[i][0]),mapy(p[i][1]));
    glEnd();
    GLfloat point[2] = {mapx(pos[0]),mapy(pos[1])};
    GLfloat c[3] = {1,0,0};
    printCircle(point,3,c);
    goBackTo3D();
}

void drawCar()
{
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslatef(pos[0],pos[1],pos[2]);
	glRotatef(-rot*180/M_PI+90,0,0,1);
	glScalef(0.06,0.06,0.06);
	
	if(carfile)
		render_mesh(carfile->meshes);
			
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
	
    gluLookAt(-eyedist*cos(-eyerot)+pos[0],
              -eyedist*sin(-eyerot)+pos[1],
              eyeheight,pos[0],pos[1],0,0,0,1);
	updateLight();
    drawBackground();
    drawRoad();
	drawCar();
    drawMap(road,nroads+1);
    GLfloat c[] = {0,0,0};
    printString("V: ",c,20,5,950);
    printInt(speed*30,c,20,50,950);
    printString("km/h",c,20,110,950);
    glFlush();
}

void reshape(int w, int h)
{
    reshapeWithoutDeform(w,h);
}

void handleMouseMoveEvent(int x, int y)
{
    trackMouse(x,y);
    glutPostRedisplay();
}

void handlePassiveMouseMoveEvent(int x, int y)
{
    trackMouse(x,y);
    glutPostRedisplay();
}

void handleKeyEvent(unsigned char keyCode, int x, int y)
{
    switch (keyCode)
    {
    case 'w':
        accelerating = true;
        break;
    case 'a':
        turningLeft = true;
        break;
    case 's':
        braking = true;
        break;
    case 'd':
        turningRight = true;
        break;
    }
    glutPostRedisplay();
}

void handleKeyUpEvent(unsigned char keyCode, int x, int y)
{
    switch (keyCode)
    {
    case 'w':
        accelerating = false;
        break;
    case 'a':
        turningLeft = false;
        break;
    case 's':
        braking = false;
        break;
    case 'd':
        turningRight = false;
        break;
    }
    glutPostRedisplay();
}

void special(int keyCode, int x, int y)
{
    switch (keyCode)
    {
    case GLUT_KEY_UP:
        accelerating = true;
        break;
    case GLUT_KEY_LEFT:
        turningLeft = true;
        break;
    case GLUT_KEY_DOWN:
        braking = true;
        break;
    case GLUT_KEY_RIGHT:
        turningRight = true;
        break;
    }
    glutPostRedisplay();
	
}

void specialUp(int keyCode, int x, int y)
{
    switch (keyCode)
    {
    case GLUT_KEY_UP:
        accelerating = false;
        break;
    case GLUT_KEY_LEFT:
        turningLeft = false;
        break;
    case GLUT_KEY_DOWN:
        braking = false;
        break;
    case GLUT_KEY_RIGHT:
        turningRight = false;
        break;
    }
    glutPostRedisplay();
	
}

void idle()
{
    advanceTimer();
    while (timepassed>1000/framePerSecond)
    {
        pos[0] += speed*cos(-rot);
        pos[1] += speed*sin(-rot);
		GLfloat friction = fabs((onroad? frictionOnRoad : frictionOffRoad)*speed);
        if (speed > 0) speed -= friction;
        if (speed < 0) speed += friction;
        if (speed < friction && speed > -friction) speed = 0;
        if (pos[0] > eastmost)
        {
            speed = -speed/2;
            pos[0] = eastmost;
        }
        if (pos[0] < westmost)
        {
            speed = -speed/2;
            pos[0] = westmost;
        }
        if (pos[1] > northmost)
        {
            speed = -speed/2;
            pos[1] = northmost;
        }
        if (pos[1] < southmost)
        {
            speed = -speed/2;
            pos[1] = southmost;
        }

        if (accelerating) speed += acceleration;
        if (turningLeft) rot -= turnRate;
        if (turningRight) rot += turnRate;
		if(rot > 2*M_PI) rot -= 2*M_PI;
		if(rot < 0) rot += 2*M_PI;
		
		if(fabs(eyerot-rot)<turnRate*3/4)
			eyerot = rot;
		else if(eyerot < rot ^ fabs(eyerot-rot) < M_PI)
			eyerot -= turnRate*3/4;
		else
			eyerot += turnRate*3/4;
		
		if(eyerot > 2*M_PI) eyerot -= 2*M_PI;
		if(eyerot < 0) eyerot += 2*M_PI;
		
        if (braking)
        {
            if (speed > 0) speed -= brake;
            if (speed < 0) speed += brake;
            if (speed < brake && speed > -brake) speed = 0;
        }

        if (crash)
        {
            speed = -speed/3;
            pos[0] += 10*speed*cos(-rot);
            pos[1] += 10*speed*sin(-rot);
            crash = false;
        }
        timepassed -= 1000/framePerSecond;
    }
    glutPostRedisplay();
}

int main(int argc, char *argv[])
{
    init("Car",&argc,argv);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMotionFunc(handleMouseMoveEvent);
    glutPassiveMotionFunc(handlePassiveMouseMoveEvent);
    glutKeyboardFunc(handleKeyEvent);
    glutKeyboardUpFunc(handleKeyUpEvent);
		glutSpecialFunc(special);
		glutSpecialUpFunc(specialUp);
    glutIdleFunc(idle);
    glutMainLoop();
    return 0;
}
