//------------------------------------------------------------------------------
#include "chai3d.h"
#include "GEL3D.h"
#include "tetgen.h"
#include <set>
//------------------------------------------------------------------------------
using namespace chai3d;
using namespace std;
//------------------------------------------------------------------------------
#ifndef MACOSX
#include "GL/glut.h"
#else
#include "GLUT/glut.h"
#endif
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// 常规设置
//------------------------------------------------------------------------------

cStereoMode stereoMode = C_STEREO_DISABLED;

// 全屏
bool fullscreen = false;

// 镜像模式
bool mirroredDisplay = false;

// Tetgen参数
char TETGEN_SWITCHES[] = "pq445.000a0.002";


//---------------------------------------------------------------------------
// 变量声明
//---------------------------------------------------------------------------

cWorld* world;
cCamera* camera;
cDirectionalLight *light;
cHapticDeviceHandler* handler;
shared_ptr<cGenericHapticDevice> hapticDevice;
double deviceForceScale;
double workspaceScaleFactor;
double cursorWorkspaceRadius;
cLabel* labelHapticRate;
bool simulationRunning = false;
bool simulationFinished = true;
cFrequencyCounter frequencyCounter;

int screenW;
int screenH;
int windowW;
int windowH;
int windowPosX;
int windowPosY;

string resourceRoot;

double groundLevel = -0.4;


//---------------------------------------------------------------------------
// 宏定义
//---------------------------------------------------------------------------

// 资源路径转化
#define RESOURCE_PATH(p)    (char*)((resourceRoot+string(p)).c_str())


//---------------------------------------------------------------------------
// 3d
//---------------------------------------------------------------------------

cGELWorld* defWorld;
cGELMesh* ground;
cGELMesh* defObject;
cShapeSphere* device;

double deviceRadius;
double radius;
double stiffness;


//---------------------------------------------------------------------------
// 函数声明
//---------------------------------------------------------------------------

void resizeWindow(int w, int h);
void keySelect(unsigned char key, int x, int y);
void updateGraphics(void);
void graphicsTimer(int data);
void close(void);
void updateHaptics(void);

cVector3d computeForce(const cVector3d& a_cursor,
	double a_cursorRadius,
	const cVector3d& a_spherePos,
	double a_radius,
	double a_stiffness);

bool createTetGenMesh(cGELMesh *a_object, char *a_filename, char *a_filenameHighRes);

bool createSkeletonMesh(cGELMesh *a_object, char *a_filename, char *a_filenameHighRes);


//===========================================================================
/*
tetgen
*/
//===========================================================================

int main(int argc, char* argv[])
{
	//-----------------------------------------------------------------------
	// 初始化工作
	//-----------------------------------------------------------------------

	cout << endl;
	cout << "-----------------------------------" << endl;
	cout << "tetgen三维显示测试" << endl;
	cout << "for 四面体窌分" << endl;
	cout << "-----------------------------------" << endl << endl << endl;
	cout << "快捷键:" << endl << endl;
	cout << "[s] - 显示或者隐藏蒙皮" << endl;
	cout << "[m] - 上下颠倒" << endl;
	cout << "[x] - 退出程序" << endl;
	cout << endl << endl;

	resourceRoot = string(argv[0]).substr(0, string(argv[0]).find_last_of("/\\") + 1);

	glutInit(&argc, argv);

	screenW = glutGet(GLUT_SCREEN_WIDTH);
	screenH = glutGet(GLUT_SCREEN_HEIGHT);
	windowW = 0.8 * screenH;
	windowH = 0.5 * screenH;
	windowPosY = (screenH - windowH) / 2;
	windowPosX = windowPosY;

	glutInitWindowPosition(windowPosX, windowPosY);
	glutInitWindowSize(windowW, windowH);

	if (stereoMode == C_STEREO_ACTIVE)
		glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE | GLUT_STEREO);
	else
		glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);

	glutCreateWindow(argv[0]);

#ifdef GLEW_VERSION
	glewInit();
#endif

	glutDisplayFunc(updateGraphics);
	glutKeyboardFunc(keySelect);
	glutReshapeFunc(resizeWindow);
	glutSetWindowTitle("tetgen测试");

	if (fullscreen)
	{
		glutFullScreen();
	}


	//-----------------------------------------------------------------------
	//-----------------------------------------------------------------------

	world = new cWorld();

	world->setBackgroundColor(1.0, 1.0, 1.0);

	camera = new cCamera(world);
	world->addChild(camera);

	camera->set(cVector3d(1.5, 0.0, 1.6),    
		cVector3d(0.5, 0.0, 0.0),    
		cVector3d(0.0, 0.0, 1.0));  

	camera->setClippingPlanes(0.01, 10.0);
	camera->setStereoMode(stereoMode);
	camera->setStereoEyeSeparation(0.02);
	camera->setStereoFocalLength(3.0);
	camera->setMirrorVertical(mirroredDisplay);
	camera->setUseMultipassTransparency(true);

	light = new cDirectionalLight(world);

	world->addChild(light);

	light->setEnabled(true);

	light->setDir(-1.0, -1.0, -1.0);


	handler = new cHapticDeviceHandler();

	handler->getDevice(hapticDevice, 0);

	cHapticDeviceInfo hapticDeviceInfo = hapticDevice->getSpecifications();

	hapticDevice->open();

	cursorWorkspaceRadius = 0.8;

	workspaceScaleFactor = cursorWorkspaceRadius / hapticDeviceInfo.m_workspaceRadius;

	deviceForceScale = 0.05 * hapticDeviceInfo.m_maxLinearForce;

	deviceRadius = 0.12;
	device = new cShapeSphere(deviceRadius);
	world->addChild(device);
	device->m_material->setWhite();
	device->m_material->setShininess(100);

	stiffness = 10;

	defWorld = new cGELWorld();
	world->addChild(defWorld);

	ground = new cGELMesh();

	defWorld->m_gelMeshes.push_back(ground);

	ground->m_useMassParticleModel = true;

	cGELMassParticle::s_default_mass = 0.010;
	cGELMassParticle::s_default_kDampingPos = 0.4;
	cGELMassParticle::s_default_gravity.set(0, 0, 0);

	cMesh* mesh = ground->newMesh();

	int RESOLUTION = 15;
	double SIZE = 5.0;
	for (int v = 0; v<RESOLUTION; v++)
	{
		for (int u = 0; u<RESOLUTION; u++)
		{
			double px, py, tu, tv;

			px = SIZE / (double)RESOLUTION * (double)u - (SIZE / 2.0);
			py = SIZE / (double)RESOLUTION * (double)v - (SIZE / 2.0);

			unsigned int index = mesh->newVertex(px, py, groundLevel);

			tu = (double)u / (double)RESOLUTION;
			tv = (double)v / (double)RESOLUTION;
			mesh->m_vertices->setTexCoord(index, tu, tv);
			mesh->m_vertices->setColor(index, cColorf(1.0, 0.0, 0.1));
		}
	}

	ground->buildVertices();

	for (int v = 0; v<RESOLUTION; v++)
	{
		for (int u = 0; u<RESOLUTION; u++)
		{
			if ((u == 0) || (v == 0) || (u == (RESOLUTION - 1)) || (v == (RESOLUTION - 1)))
			{
				unsigned int index = ((v + 0) * RESOLUTION) + (u + 0);
				ground->m_gelVertices[index].m_massParticle->m_fixed = true;
			}
		}
	}

	cGELLinearSpring::s_default_kSpringElongation = 10.0; // [N/m]

	for (int v = 0; v<(RESOLUTION - 1); v++)
	{
		for (int u = 0; u<(RESOLUTION - 1); u++)
		{
			unsigned int index00 = ((v + 0) * RESOLUTION) + (u + 0);
			unsigned int index01 = ((v + 0) * RESOLUTION) + (u + 1);
			unsigned int index10 = ((v + 1) * RESOLUTION) + (u + 0);
			unsigned int index11 = ((v + 1) * RESOLUTION) + (u + 1);

			mesh->newTriangle(index00, index01, index10);
			mesh->newTriangle(index10, index01, index11);

			cGELMassParticle* m0 = ground->m_gelVertices[index00].m_massParticle;
			cGELMassParticle* m1 = ground->m_gelVertices[index01].m_massParticle;
			cGELMassParticle* m2 = ground->m_gelVertices[index10].m_massParticle;
			cGELMassParticle* m3 = ground->m_gelVertices[index11].m_massParticle;

			cGELLinearSpring* spring0 = new cGELLinearSpring(m0, m1);
			cGELLinearSpring* spring1 = new cGELLinearSpring(m0, m2);
			ground->m_linearSprings.push_back(spring0);
			ground->m_linearSprings.push_back(spring1);

			if ((u == (RESOLUTION - 2)) || (v == (RESOLUTION - 2)))
			{
				cGELLinearSpring* spring2 = new cGELLinearSpring(m3, m1);
				cGELLinearSpring* spring3 = new cGELLinearSpring(m3, m2);
				ground->m_linearSprings.push_back(spring2);
				ground->m_linearSprings.push_back(spring3);
			}
		}
	}

	ground->setUseMaterial(true);
	ground->m_material->setGrayLevel(0.1);
	ground->setTransparencyLevel(0.7);
	ground->setUseTransparency(true);

	shared_ptr<cTexture2d> textureGround(new cTexture2d());
	ground->setTexture(textureGround);
	ground->setUseTexture(true, true);

	bool fileload;
	fileload = textureGround->loadFromFile(RESOURCE_PATH("../resources/images/water.jpg"));
	if (!fileload)
	{
#if defined(_MSVC)
		fileload = textureGround->loadFromFile("../../../bin/resources/images/water.jpg");
#endif
		if (!fileload)
		{
			cout << "加载3d模型错误." << endl;
			close();
			return (-1);
		}
	}

	textureGround->setEnvironmentMode(GL_DECAL);
	textureGround->setSphericalMappingEnabled(true);

	defObject = new cGELMesh();

	defWorld->m_gelMeshes.push_back(defObject);

	fileload = createSkeletonMesh(defObject, RESOURCE_PATH("../resources/models/ducky/duck.off"), RESOURCE_PATH("../resources/models/ducky/duck.obj"));
	if (!fileload)
	{
#if defined(_MSVC)
		fileload = createSkeletonMesh(defObject, "../../../bin/resources/models/ducky/duck.off", "../../../bin/resources/models/ducky/duck.obj");
#endif
		if (!fileload)
		{
			cout << "加载模型错误." << endl;
			close();
			return (-1);
		}
	}

	cFont *font = NEW_CFONTCALIBRI20();

	labelHapticRate = new cLabel(font);
	camera->m_frontLayer->addChild(labelHapticRate);
	labelHapticRate->m_fontColor.setWhite();

	cBackground* background = new cBackground();
	camera->m_backLayer->addChild(background);

	fileload = background->loadFromFile(RESOURCE_PATH("../resources/images/stone.jpg"));
	if (!fileload)
	{
#if defined(_MSVC)
		fileload = background->loadFromFile("../../../bin/resources/images/stone.jpg");
#endif
	}
	if (!fileload)
	{
		cout << "载入图像错误" << endl;
		close();
		return (-1);
	}

	background->setFixedAspectRatio(true);

	cThread* hapticsThread = new cThread();
	hapticsThread->start(updateHaptics, CTHREAD_PRIORITY_HAPTICS);

	atexit(close);

	glutTimerFunc(50, graphicsTimer, 0);
	glutMainLoop();

	return (0);
}

//---------------------------------------------------------------------------

bool createSkeletonMesh(cGELMesh *a_object, char *a_filename, char *a_filenameHighRes)
{
	a_object->m_useSkeletonModel = true;
	a_object->m_useMassParticleModel = false;
	a_object->loadFromFile(a_filenameHighRes);

	cGELMesh* model = new cGELMesh();
	cMesh* mesh = model->newMesh();

	tetgenio input;
	if (input.load_off(a_filename))
	{
		// 使用TetGen来刨分我们的网格
		tetgenio output;
		tetrahedralize(TETGEN_SWITCHES, &input, &output);

		for (int p = 0, pi = 0; p < output.numberofpoints; ++p, pi += 3)
		{
			cVector3d point;
			point.x(output.pointlist[pi + 0]);
			point.y(output.pointlist[pi + 1]);
			point.z(output.pointlist[pi + 2]);
			mesh->newVertex(point);
		}

		for (int t = 0, ti = 0; t < output.numberoftrifaces; ++t, ti += 3)
		{
			cVector3d p[3];
			unsigned int vi[3];
			for (int i = 0; i < 3; ++i)
			{
				int tc = output.trifacelist[ti + i];
				vi[i] = tc;
				int pi = tc * 3;
				p[i].x(output.pointlist[pi + 0]);
				p[i].y(output.pointlist[pi + 1]);
				p[i].z(output.pointlist[pi + 2]);
			}
			mesh->newTriangle(vi[1], vi[0], vi[2]);
		}

		set<int> inside, outside;
		for (int t = 0; t < output.numberoftrifaces * 3; ++t)
		{
			outside.insert(output.trifacelist[t]);
		}
		for (int p = 0; p < output.numberofpoints; ++p)
		{
			if (outside.find(p) == outside.end())
				inside.insert(p);
		}

		model->computeAllNormals();

		model->computeBoundaryBox(true);

		double size = cSub(model->getBoundaryMax(), model->getBoundaryMin()).length();

		if (size > 0)
		{
			model->scale(1.5 / size);
			a_object->scale(1.5 / size);
		}

		cGELSkeletonNode::s_default_radius = 0.05;
		cGELSkeletonNode::s_default_kDampingPos = 0.3;
		cGELSkeletonNode::s_default_kDampingRot = 0.1;
		cGELSkeletonNode::s_default_mass = 0.002;  // [kg]
		cGELSkeletonNode::s_default_showFrame = false;
		cGELSkeletonNode::s_default_color.set(1.0, 0.6, 0.6);
		cGELSkeletonNode::s_default_useGravity = true;
		cGELSkeletonNode::s_default_gravity.set(0.00, 0.00, -3.45);
		radius = cGELSkeletonNode::s_default_radius;

		a_object->buildVertices();
		model->buildVertices();

		vector<cGELSkeletonNode*> nodes;
		int i = 0;
		for (set<int>::iterator it = inside.begin(); it != inside.end(); ++it)
		{
			cGELSkeletonNode* newNode = new cGELSkeletonNode();
			a_object->m_nodes.push_front(newNode);

			unsigned int vertexIndex = 0;
			cMesh* mesh = NULL;

			if (model->getVertex(*it, mesh, vertexIndex))
			{
				newNode->m_pos = mesh->m_vertices->getLocalPos(vertexIndex);
				newNode->m_rot.identity();
				newNode->m_radius = 0.1;
				newNode->m_fixed = false;
				mesh->m_vertices->setUserData(vertexIndex, i);
				i++;
				nodes.push_back(newNode);
			}
		}

		set< pair<int, int> > springs;
		for (int t = 0, ti = 0; t < output.numberoftetrahedra; ++t, ti += 4)
		{
			for (int i = 0; i < 4; ++i) {
				int v0 = output.tetrahedronlist[ti + i];
				for (int j = i + 1; j < 4; ++j) {
					int v1 = output.tetrahedronlist[ti + j];

					if (inside.find(v0) != inside.end() && inside.find(v1) != inside.end())
						springs.insert(pair<int, int>(min(v0, v1), max(v0, v1)));
				}
			}
		}

		cGELSkeletonLink::s_default_kSpringElongation = 100.0; // [N/m]
		cGELSkeletonLink::s_default_kSpringFlexion = 0.1;   // [Nm/RAD]
		cGELSkeletonLink::s_default_kSpringTorsion = 0.1;   // [Nm/RAD]
		cGELSkeletonLink::s_default_color.set(0.2, 0.2, 1.0);

		for (set< pair<int, int> >::iterator it = springs.begin(); it != springs.end(); ++it)
		{
			unsigned int vertexIndex0 = 0;
			unsigned int vertexIndex1 = 0;
			cMesh* mesh0 = NULL;
			cMesh* mesh1 = NULL;

			model->getVertex(it->first, mesh0, vertexIndex0);
			model->getVertex(it->second, mesh1, vertexIndex1);
			cGELSkeletonNode* n0 = nodes[mesh0->m_vertices->getUserData(vertexIndex0)];
			cGELSkeletonNode* n1 = nodes[mesh1->m_vertices->getUserData(vertexIndex1)];
			cGELSkeletonLink* newLink = new cGELSkeletonLink(n0, n1);
			a_object->m_links.push_front(newLink);
		}

		a_object->connectVerticesToSkeleton(false);

		cMaterial mat;
		mat.m_ambient.set(0.7, 0.7, 0.7);
		mat.m_diffuse.set(0.8, 0.8, 0.8);
		mat.m_specular.set(0.0, 0.0, 0.0);
		a_object->setMaterial(mat, true);

		return (true);
	}
	return (false);
}

//---------------------------------------------------------------------------

void resizeWindow(int w, int h)
{
	windowW = w;
	windowH = h;
}

//---------------------------------------------------------------------------

void keySelect(unsigned char key, int x, int y)
{
	if ((key == 27) || (key == 'x'))
	{
		exit(0);
	}

	if (key == 's')
	{
		defObject->m_showSkeletonModel = !defObject->m_showSkeletonModel;

		if (defObject->m_showSkeletonModel)
		{
			defObject->setWireMode(true, true);
		}
		else
		{
			defObject->setWireMode(false, true);
		}
	}

	if (key == 'f')
	{
		if (fullscreen)
		{
			windowPosX = glutGet(GLUT_INIT_WINDOW_X);
			windowPosY = glutGet(GLUT_INIT_WINDOW_Y);
			windowW = glutGet(GLUT_INIT_WINDOW_WIDTH);
			windowH = glutGet(GLUT_INIT_WINDOW_HEIGHT);
			glutPositionWindow(windowPosX, windowPosY);
			glutReshapeWindow(windowW, windowH);
			fullscreen = false;
		}
		else
		{
			glutFullScreen();
			fullscreen = true;
		}
	}

	if (key == 'm')
	{
		mirroredDisplay = !mirroredDisplay;
		camera->setMirrorVertical(mirroredDisplay);
	}
}


//---------------------------------------------------------------------------

void close(void)
{
	simulationRunning = false;

	while (!simulationFinished) { cSleepMs(100); }

	hapticDevice->close();
}

//---------------------------------------------------------------------------

void graphicsTimer(int data)
{
	if (simulationRunning)
	{
		glutPostRedisplay();
	}

	glutTimerFunc(50, graphicsTimer, 0);
}

//---------------------------------------------------------------------------

void updateGraphics(void)
{
	labelHapticRate->setText(cStr(frequencyCounter.getFrequency(), 0) + " Hz");

	labelHapticRate->setLocalPos((int)(0.5 * (windowW - labelHapticRate->getWidth())), 15);

	defWorld->updateSkins(true);
	world->updateShadowMaps(false, mirroredDisplay);

	camera->renderView(windowW, windowH);

	glutSwapBuffers();

	glFinish();

	GLenum err;
	err = glGetError();
	if (err != GL_NO_ERROR) cout << "错误: " << gluErrorString(err) << endl;
}

//---------------------------------------------------------------------------

void updateHaptics(void)
{
	frequencyCounter.reset();

	cPrecisionClock clock;
	clock.reset();

	double time = 0.0;

	simulationRunning = true;
	simulationFinished = false;

	while (simulationRunning)
	{
		double interval = 0.001; // cMin(0.001, clock.stop());

		clock.start(true);

		time = time + interval;

		cVector3d pos;
		hapticDevice->getPosition(pos);
		pos.mul(workspaceScaleFactor);
		device->setLocalPos(pos);

		cVector3d force(0.0, 0.0, 0.0);

		list<cGELMesh*>::iterator i;

		for (i = defWorld->m_gelMeshes.begin(); i != defWorld->m_gelMeshes.end(); ++i)
		{
			cGELMesh *nextItem = *i;

			if (nextItem->m_useMassParticleModel)
			{
				int numVertices = (int)(nextItem->m_gelVertices.size());
				for (int i = 0; i<numVertices; i++)
				{
					cVector3d nodePos = nextItem->m_gelVertices[i].m_massParticle->m_pos;

					double forceWave = 0.001 * sin(1.0 * (time + nodePos.x() + nodePos.y()));

					cVector3d force = cVector3d(-0.002 * nodePos.x(), -0.002 * nodePos.y(), forceWave);
					if (nodePos.z() < groundLevel)
					{
						double depth = nodePos.z() - groundLevel;
						force.add(cVector3d(0, 0, -100 * depth));
					}
					nextItem->m_gelVertices[i].m_massParticle->setExternalForce(force);
				}
			}

			if (nextItem->m_useSkeletonModel)
			{
				list<cGELSkeletonNode*>::iterator i;
				for (i = nextItem->m_nodes.begin(); i != nextItem->m_nodes.end(); ++i)
				{
					cGELSkeletonNode* node = *i;
					cVector3d nodePos = node->m_pos;
					double radius = node->m_radius;
					cVector3d force = cVector3d(-0.01 * nodePos.x(), -0.01 * (nodePos.y()), 0.0);

					if ((nodePos.z() - radius) < groundLevel)
					{
						double depth = (nodePos.z() - radius) - groundLevel;
						force.add(cVector3d(0, 0, -1.0 * depth));
						node->m_vel.mul(0.95);
					}

					double forceWave = 0.001 * sin(time);
					force.add(0.0, forceWave, 0.0);

					node->setExternalForce(force);
				}
			}
		}

		for (i = defWorld->m_gelMeshes.begin(); i != defWorld->m_gelMeshes.end(); ++i)
		{
			cGELMesh *nextItem = *i;

			if (nextItem->m_useMassParticleModel)
			{
				int numVertices = (int)(nextItem->m_gelVertices.size());
				for (int i = 0; i<numVertices; i++)
				{
					cVector3d nodePos = nextItem->m_gelVertices[i].m_massParticle->m_pos;
					cVector3d f = computeForce(pos, deviceRadius, nodePos, radius, stiffness);
					if (f.lengthsq() > 0)
					{
						cVector3d tmpfrc = cNegate(f);
						nextItem->m_gelVertices[i].m_massParticle->setExternalForce(tmpfrc);
					}
					force.add(cMul(1.0, f));
				}
			}

			if (nextItem->m_useSkeletonModel)
			{
				list<cGELSkeletonNode*>::iterator i;
				for (i = nextItem->m_nodes.begin(); i != nextItem->m_nodes.end(); ++i)
				{
					cGELSkeletonNode* node = *i;
					cVector3d nodePos = node->m_pos;
					double radius = node->m_radius;
					cVector3d f = computeForce(pos, deviceRadius, nodePos, radius, stiffness);
					if (f.lengthsq() > 0)
					{
						cVector3d tmpfrc = cNegate(f);
						node->setExternalForce(tmpfrc);
					}
					force.add(cMul(4.0, f));
				}
			}
		}

		defWorld->updateDynamics(interval);

		force.mul(deviceForceScale);

		if ((pos.z() - deviceRadius) < groundLevel)
		{
			cHapticDeviceInfo info = hapticDevice->getSpecifications();
			double Kv = 0.8 * info.m_maxLinearDamping;

			cVector3d linearVelocity;
			hapticDevice->getLinearVelocity(linearVelocity);

			double val = (groundLevel - (pos.z() - deviceRadius)) / (2.0 * deviceRadius);
			double scale = cClamp(val, 0.1, 1.0);

			cVector3d forceDamping = cMul(-Kv * scale, linearVelocity);
			force.add(forceDamping);
		}

		hapticDevice->setForce(force);

		frequencyCounter.signal(1);
	}

	simulationFinished = true;
}

//---------------------------------------------------------------------------

cVector3d computeForce(const cVector3d& a_cursor,
	double a_cursorRadius,
	const cVector3d& a_spherePos,
	double a_radius,
	double a_stiffness)
{

	cVector3d force;
	force.zero();
	cVector3d vSphereCursor = a_cursor - a_spherePos;

	if (vSphereCursor.length() < 0.0000001)
	{
		return (force);
	}

	if (vSphereCursor.length() > (a_cursorRadius + a_radius))
	{
		return (force);
	}

	double penetrationDistance = (a_cursorRadius + a_radius) - vSphereCursor.length();
	cVector3d forceDirection = cNormalize(vSphereCursor);
	force = cMul(penetrationDistance * a_stiffness, forceDirection);

	return (force);
}
