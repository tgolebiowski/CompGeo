#ifndef TRIANGULATOR_H
#define TRIANGULATOR_H
#include <vector>
#include <algorithm>
#include <assert.h>
#include <unordered_map>
#include <random>
#define PI 3.1415926

struct Triangle;
class Triangulator;

struct Point { 
	float x, y;

	Point() {
		x = 0.0f;
		y = 0.0f;
	}

	Point(float x, float y) {
		this->x = x;
		this->y = y;
	}

	inline Point operator - (Point p) const {
		return Point(this->x - p.x, this->y - p.y);
	}

	inline Point operator + (Point p) const {
		return Point(this->x + p.x, this->y + p.y);
	}

	inline Point operator * (float s) const {
		return Point(this->x * s, this->y * s);
	}

	inline float length() const {
		return std::sqrtf(x * x + y * y);
	}

	inline float crossProd(Point p) const {
		return this->x * p.y - this->y * p.x;
	}

	static void PoissonDiscSample(const float width, const float height, 
		const float minDist, std::vector<Point>* outList);
};

struct Site {
	Point p;

	std::vector<Triangle*> incidentTris;
	std::vector<Site*> neighbors;

	void* userData;

	Site(Point p)	{
		userData = NULL;
		this->p = p;
	}

	~Site()	{

	}
};

struct Edge {
	Site* site1;
	Site* site2;

	Triangle* tri1;
	Triangle* tri2;

	Edge(Site* site1, Site* site2) {
		this->site1 = site1;
		this->site2 = site2;
		tri1 = NULL;
		tri2 = NULL;
	}
};

struct Triangle {
	Site* sites[3];
	//edges[0] : sites[0] & sites[1] /// edges[1] : sites[1] & sites[2] /// edges[2] : sites[2] & sites[0]
	Edge* edges[3];
	Triangle* neighbors[3];

	void* userData;
	void* voronoiHolder;
};

struct Circle {
	Point center;
	float radius;

	Circle() {
		radius = 0.0f;
	};

	Circle(Triangle* tri);
};

struct VoronoiNode
{
	VoronoiNode* neighbors[3];
	Point position;

	VoronoiNode()
	{
		neighbors[0] = NULL;
		neighbors[1] = NULL;
		neighbors[2] = NULL;
	}

	VoronoiNode(Triangle* baseTri)
	{
		position = Circle(baseTri).center;

		neighbors[0] = NULL;
		neighbors[1] = NULL;
		neighbors[2] = NULL;
	}
};

struct VoronoiGraph
{
	int nodeCount;
	VoronoiNode* nodes;

	VoronoiGraph(Triangulator* triangulator);
	~VoronoiGraph();
};

class Triangulator
{
public:
	Triangulator(std::vector<Site*>* siteList);
	~Triangulator();

	void insert(Site* newSite);
	void floodFill(Site* startSite, float cutoff, std::vector<Site*>* outList) const;
	Site* getNearestSite(float x, float y) const;
	Triangle* getTri(const float x, const float y, Triangle* startTri = NULL) const;

	std::vector<Edge*> edges;
	std::vector<Triangle*> tris;
	std::vector<Site*> sites;

private:
	void flip(Triangle* tri1, Triangle* tri2, Edge* edgeInCommon);
	void insertSite(Triangle* enclosingTri, Site* newSite);
	int pointInTri(Triangle* tri, float x, float y) const;

	//validation method
	static bool checkEdgeCorrespondence(Triangle* tri);

	Triangle* centerMostTri;
};
#endif