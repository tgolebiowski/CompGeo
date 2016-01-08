#pragma once
#include <Polycode.h>
#ifndef QUASICRYSTAL
#define QUASICRYSTAL

using namespace Polycode;
class QuasiCrystal
{
public:
	QuasiCrystal(int waves, float period);
	~QuasiCrystal(void);

	double get(double x, double y, bool contrastOn = true);
	void generatePointSet(double gridSize, double cellSize, float badPointBias, std::vector<Vector2>* pointSetPointer, bool addJitter = true, Vector2 offset = Vector2(0,0));

private:

	double waveAngles[21];
	double period;
	double waveCount;
	double offset;
};
#endif