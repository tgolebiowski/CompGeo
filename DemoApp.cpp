#include "DemoApp.h"
#include <time.h>
#include <random>

DemoApp::DemoApp(PolycodeView *view) {
	core = new Win32Core(view, 1300,700,false, false, 4, 0,30);	  
	CoreServices::getInstance()->getResourceManager()->addArchive("default.pak");
	CoreServices::getInstance()->getResourceManager()->addDirResource("default", false);
	core->getServices()->getRenderer()->setClearColor(Color(1,1,1,1));

	srand(time(NULL));

	scene = new Scene(Scene::SCENE_2D);
	scene->getActiveCamera()->setOrthoSize(1300,700);
	scene->getActiveCamera()->setOrthoMode(true);

	std::vector<Point> pointList;
	Point::PoissonDiscSample(600.0f, 600.0f, 20.0f, &pointList);

	crystal = new QuasiCrystal(7, 25.0f);
	//crystal->generatePointSet(600.0f, 1.0f, 9.0f, &pointList, true, Vector2(-300.0f, -300.0f));

	std::vector<Site*> sites;
	//for (int i = 0; i < pointList.size(); i++) {
	//	ScenePrimitive* dot = new ScenePrimitive(ScenePrimitive::TYPE_CIRCLE, 8, 8, 8);
	//	dot->setPosition(pointList[i].x - 300.0, pointList[i].y - 300.0);

	//	if (dot->getPosition().length() < 300.0f)
	//	{
	//		//scene->addEntity(dot);

	//		Vector2 p = dot->getPosition2D();
	//		Site* newSite = new Site(Pnt(p.x, p.y));
	//		sites.push_back(newSite);
	//	}
	//	else
	//	{
	//		//pointList.erase(pointList.begin() + i);
	//		//i--;
	//	}
	//}

	for (int i = 0; i < pointList.size(); i++) {
		Point p = pointList.at(i);
		p.x -= 300.0f;
		p.y -= 300.0f;
		if (p.length() < 300.0f) {
			sites.push_back(new Site(p));
		}
	}

	triangulator = new Triangulator(&sites);

	for (int i = 0; i < triangulator->edges.size(); i++)
	{
		Point p1 = triangulator->edges.at(i)->site1->p;
		Point p2 = triangulator->edges.at(i)->site2->p;
		SceneLine* newLine = new SceneLine(Vector3(p1.x, p1.y, 0.0f), Vector3(p2.x, p2.y, 0.0f));
		scene->addEntity(newLine);
	}

	//VoronoiGraph voronoiGraph = VoronoiGraph(triangulator);
	//for (int i = 0; i < voronoiGraph.nodeCount; i++)
	//{
	//	VoronoiNode* node = &voronoiGraph.nodes[i];

	//	for (int j = 0; j < 3; j++)
	//	{
	//		VoronoiNode* neighbor = node->neighbors[j];
	//		if (neighbor == NULL) continue;

	//		Vector3 p1 = Vector3(node->position.x, node->position.y, 0.0f);
	//		Vector3 p2 = Vector3(neighbor->position.x, neighbor->position.y, 0.0f);

	//		SceneLine* line = new SceneLine(p1, p2);
	//		scene->addEntity(line);
	//	}
	//}

	click = false;
	cursor = new ScenePrimitive(ScenePrimitive::TYPE_CIRCLE, 12, 12, 8);
	cursor->setColor(0.78, 0.02, 0.91, 1.0);
	scene->addEntity(cursor);
	triMesh = new SceneMesh(Mesh::TRI_MESH);
	
	scene->addEntity(triMesh);

	fpsLabel = new SceneLabel("Fps: ", 14);
	fpsLabel->setPosition(600, 325);

	scene->addEntity(fpsLabel);
}

DemoApp::~DemoApp() {
    
}

bool DemoApp::Update() {

	Vector2 m = core->getInput()->getMousePosition();
	Vector2 p = m - Vector2(1300.0f/2.0f, 700.0f/2.0f);
	cursor->setPosition(p.x, -p.y, 0.0f);

	Vector2 testPoint = Vector2(p.x, -p.y);

	static SceneMesh* tri = NULL;
	if (tri != NULL) {
		scene->removeEntity(tri);
		delete tri;
		tri = NULL;
	}
	Triangle* triangle = triangulator->getTri(testPoint.x, testPoint.y);
	if (triangle != NULL) {
		Point p1 = triangle->sites[0]->p;
		Point p2 = triangle->sites[1]->p;
		Point p3 = triangle->sites[2]->p;
		if ((p3 - p2).crossProd(p1 - p2) < 0.0)	{
			Point holder = p2;
			p2 = p3;
			p3 = holder;
		}
		tri = new SceneMesh(Mesh::TRI_MESH);
		tri->getMesh()->addVertex(p1.x, p1.y, 0.0f);
		tri->getMesh()->addVertex(p2.x, p2.y, 0.0f);
		tri->getMesh()->addVertex(p3.x, p3.y, 0.0f);
		scene->addEntity(tri);
	}

	String fpsText = "FPS: " + String::NumberToString(core->getFPS(), 0);
	fpsLabel->setText(fpsText);

	return core->updateAndRender();
}