#include <Polycode.h>
#include "DemoApp.h"
#include "PolycodeView.h"
#include "windows.h"

using namespace Polycode;

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
	PolycodeView *view = new PolycodeView(hInstance, nCmdShow, L"Comp Geo Demo App");
	DemoApp *app = new DemoApp(view);

	MSG Msg;
	do {
		while (PeekMessage(&Msg, NULL, 0, 0, PM_REMOVE)) {
			TranslateMessage(&Msg);
			DispatchMessage(&Msg);
		}
	} while(app->Update());
	return Msg.wParam;
}