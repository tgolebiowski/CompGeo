#include "PolycodeView.h"
#include "Polycode.h"
#include "Triangulator.h"
#include "QuasiCrystal.h"

using namespace Polycode;

class DemoApp {
public:
    DemoApp(PolycodeView *view);
    ~DemoApp();
    
    bool Update();

private:
    Core *core;
	Scene* scene;
	QuasiCrystal* crystal;
	Triangulator* triangulator;

	ScenePrimitive* cursor;
	SceneMesh* triMesh;
	std::vector<SceneMesh*> tris;
	std::vector<ScenePrimitive*> sites;
	std::vector<SceneLine*> debugLines;

	SceneLabel* fpsLabel;
	SceneLabel* timeCheck1;
	SceneLabel* timeCheck2;
	bool click;
};