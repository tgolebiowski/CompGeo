#include "Triangulator.h"

void Point::PoissonDiscSample(const float width, const float height, const float minDist, std::vector<Point>* outList)
{
	//1.4142 = sqrt(2)
	float cellsize = minDist/1.4142;

	//generate grid
	int gridWidth = (int)(width/cellsize) + 1;
	int gridHeight = (int)(height/cellsize) + 1;
	int maxPoints = gridWidth * gridHeight;
	Point* grid = new Point[maxPoints];

	std::vector<Point> processList;
	processList.push_back(Point(rand() % gridWidth, rand() % gridHeight)); //first point

	while(!processList.empty() && outList->size() < maxPoints) {
		Point point = processList.back();
		processList.pop_back();

		for(int i = 0; i < 15; i++) {
			//Generate point around
			float rand1 = ((float)rand() / (RAND_MAX + 1));
			float rand2 = ((float)rand() / (RAND_MAX + 1));
			float randRadius = minDist * (rand1 + 1);
			float randAngle = 2 * PI * rand2;
			float newX = point.x + randRadius * cos(randAngle);
			float newY = point.y + randRadius * sin(randAngle);
			Point newPoint = Point(newX, newY);
			//end Generate point Around

			if(newPoint.x > 0.0f && newPoint.x < width && newPoint.y > 0.0f && newPoint.y < height) {

				auto imageToGrid = [](Point point, float cellsize) -> Point {
					int gridX = point.x / cellsize;
					int gridY = point.y / cellsize;
					return Point(gridX, gridY);
				};

				bool inNeighborhoodTest = false;
				Point gridPoint = imageToGrid(newPoint, cellsize);
				//check in cells around center cell, clamp to inside of total grid
				int areaXMax = std::min((int)gridPoint.x + 2, gridWidth - 1);
				int areaXMin = std::max((int)gridPoint.x - 2, 0);
				int areaYMax = std::min((int)gridPoint.y + 2, gridWidth - 1);
				int areaYMin = std::max((int)gridPoint.y - 2, 0);

				//search loop
				for (int i = areaYMin; i <= areaYMax; i++) {
					for (int j = areaXMin; j <= areaXMax; j++) {
						Point otherPoint = grid[j + i * gridWidth];

						if (otherPoint.x != 0 && otherPoint.y != 0) {
							float dist = (otherPoint - newPoint).length();

							if (dist < minDist) {
								//set flag to true if there is a point too close
								inNeighborhoodTest = true;
								goto breakout;
							}
						}
					}
				}

			breakout:

				if(!inNeighborhoodTest) {
					//update containers
					processList.push_back(newPoint);
					outList->push_back(newPoint);

					Point gridSpot = imageToGrid(newPoint, cellsize);
					grid[((int)gridSpot.x) + ((int)gridSpot.y) * gridWidth] = newPoint;
				}
			}
		}
	}
}

Triangulator::Triangulator(std::vector<Site*>* siteList)
{
	centerMostTri = NULL;

	float minX = siteList->front()->p.x;
	float maxX = siteList->front()->p.x;
	float minY = siteList->front()->p.y;
	float maxY = siteList->front()->p.y;

	for (int i = 0; i < siteList->size(); i++)
	{
		Point p = siteList->at(i)->p;

		maxX = std::max((float)p.x, maxX);
		minX = std::min((float)p.x, minX);
		maxY = std::max((float)p.y, maxY);
		minY = std::min((float)p.y, minY);
	}

	minX *= 1.05f;
	maxX *= 1.05f;
	maxY *= 1.05f;
	minY *= 1.05f;

	std::vector<Site*> startingSites;
	Site* tl = new Site(Point(minX, maxY));
	Site* tr = new Site(Point(maxX, maxY));
	Site* bl = new Site(Point(minX, minY));
	Site* br = new Site(Point(maxX, minY));
	startingSites.push_back(tl); startingSites.push_back(tr); 
	startingSites.push_back(bl); startingSites.push_back(br);

	Edge* tl_tr = new Edge(tl, tr);
	Edge* tr_br = new Edge(tr, br);
	Edge* br_bl = new Edge(br, bl);
	Edge* bl_tl = new Edge(bl, tl);
	Edge* bl_tr = new Edge(bl, tr);
	edges.push_back(tl_tr);
	edges.push_back(tr_br);
	edges.push_back(br_bl);
	edges.push_back(bl_tl);
	edges.push_back(bl_tr);

	Triangle* tri1 = new Triangle(); tri1->edges[0] = tl_tr; tri1->edges[1] = bl_tl; tri1->edges[2] = bl_tr;
	Triangle* tri2 = new Triangle(); tri2->edges[0] = br_bl; tri2->edges[1] = tr_br; tri2->edges[2] = bl_tr;
	tl_tr->tri1 = tri1; tl_tr->tri2 = NULL;
	tr_br->tri1 = tri2; tr_br->tri2 = NULL;
	br_bl->tri1 = tri2; br_bl->tri2 = NULL;
	bl_tl->tri1 = tri1; bl_tl->tri2 = NULL;
	bl_tr->tri1 = tri1; bl_tr->tri2 = tri2;
	tri1->neighbors[0] = NULL; tri1->neighbors[1] = tri2; tri1->neighbors[2] = NULL;
	tri2->neighbors[0] = NULL; tri2->neighbors[1] = NULL; tri2->neighbors[2] = tri1;

	for (int i = 0; i < 2; i++)
	{
		Triangle* tri = NULL;
		if (i == 0) tri = tri1;
		else tri = tri2;

		tri->sites[0] = tri->edges[0]->site1;
		tri->sites[1] = tri->edges[0]->site2;
		if (tri->edges[1]->site1 != tri->sites[0] && tri->edges[1]->site1 != tri->sites[1])
			tri->sites[2] = tri->edges[1]->site1;
		else
			tri->sites[2] = tri->edges[1]->site2;
	}

	tris.push_back(tri1);
	tris.push_back(tri2);

	for (int i = 0; i < siteList->size(); i++)
	{
		insert(siteList->at(i));
	}

	std::vector<Edge*> edgesToDelete;
	std::vector<Triangle*> trisToDelete;
	//collect triangles that are touching the starting points
	for (int i = 0; i < startingSites.size(); i++)
	{
		Site* site = startingSites.at(i);
		for (int j = 0; j < site->incidentTris.size(); j++)
		{
			Triangle* toInsert = site->incidentTris.at(j);
			if (std::find(trisToDelete.begin(), trisToDelete.end(), toInsert) == trisToDelete.end())
			{
				trisToDelete.push_back(toInsert);
				for (int k = 0; k < 3; k++)
				{
					Triangle* neighbor = toInsert->neighbors[k];
					for (int l = 0; neighbor != NULL && l < 3; l++)
					{
						if (neighbor->neighbors[l] == toInsert)
						{
							neighbor->neighbors[l] = NULL;
							break;
						}
					}
				}
			}
		}
	}
	//collect the edges w/n above triangles that have one of the starting sites as one of their two sites
	for (int i = 0; i < trisToDelete.size(); i++)
	{
		Triangle* tri = trisToDelete.at(i);

		for (int j = 0; j < 3; j++)
		{
			Edge* e = trisToDelete.at(i)->edges[j];
			bool connectedToStartSite = false;
			for (int k = 0; k < startingSites.size(); k++)
			{
				connectedToStartSite = e->site1 == startingSites.at(k) || e->site2 == startingSites.at(k);
				if (connectedToStartSite) break;
			}

			if (connectedToStartSite &&
				std::find(edgesToDelete.begin(), edgesToDelete.end(), e) == edgesToDelete.end())
				edgesToDelete.push_back(e);
			else
			{
				if (e->tri1 == tri)
				{
					for (int k = 0; k < 3; k++)
					{
						if (e->tri2->neighbors[k] == tri)
						{
							e->tri2->neighbors[k] = NULL;
							break;
						}
					}
					e->tri1 = NULL;
				}
				else
				{
					for (int k = 0; k < 3; k++)
					{
						if (e->tri1->neighbors[k] == tri)
						{
							e->tri1->neighbors[k] = NULL;
							break;
						}
					}
					e->tri2 = NULL;
				}
			}
		}
	}

	for (int i = 0; i < trisToDelete.size(); i++)
	{
		tris.erase(std::find(tris.begin(), tris.end(), trisToDelete.at(i)));
		for (int j = 0; j < 3; j++)
		{
			Site* site = trisToDelete.at(i)->sites[j];
			for (int k = 0; k < site->incidentTris.size(); k++)
			{
				if (site->incidentTris[k] == trisToDelete[i])
				{
					site->incidentTris.erase(site->incidentTris.begin() + k);
					break;
				}
			}
		}

		delete trisToDelete.at(i);
	}
	trisToDelete.clear();

	for (int i = 0; i < edgesToDelete.size(); i++)
	{
		edges.erase(std::find(edges.begin(), edges.end(), edgesToDelete.at(i)));
		delete edgesToDelete.at(i);
	}
	edgesToDelete.clear();

	for (int i = 0; i < startingSites.size(); i++)
	{
		Site* startSite = startingSites.at(i);

		for (int j = 0; j < startSite->neighbors.size(); j++)
		{
			Site* neighbor = startSite->neighbors.at(j);
			neighbor->neighbors.erase(std::find(neighbor->neighbors.begin(),
				neighbor->neighbors.end(), startSite));
		}
	}

	delete tl;
	delete tr;
	delete bl;
	delete br;
	startingSites.clear();

	float avgX = (maxX + minX) / 2.0f;
	float avgY = (maxY + minY) / 2.0f;
	Point center = Point(avgX, avgY);
	float shortestDistance = (center - Point(maxX, maxY)).length();
	for (int i = 0; i < tris.size(); i++)
	{
		Circle c = Circle(tris.at(i));
		float distance = (c.center - center).length();
		if (distance < shortestDistance)
			centerMostTri = tris.at(i);
	}
}

Triangulator::~Triangulator()
{
	for (int i = 0; i < sites.size(); i++)
		delete sites[i];
	sites.clear();

	for (int i = 0; i < edges.size(); i++)
		delete edges[i];
	edges.clear();

	for (int i = 0; i < tris.size(); i++)
		delete tris[i];
	tris.clear();
}

Triangle* Triangulator::getTri(const float x, const float y, Triangle* startTri) const
{
	if (startTri == NULL)
		startTri = centerMostTri;
	if (startTri == NULL)
		startTri = tris.front();

	int result = -2;
	Triangle* currentTri = startTri;
	do{
		Triangle* tri = currentTri;
		Point p1 = tri->sites[0]->p;
		Point p2 = tri->sites[1]->p;
		Point p3 = tri->sites[2]->p;

		float denominator = ((p2.y - p3.y) * (p1.x - p3.x) + (p3.x - p2.x) * (p1.y - p3.y));

		//a corresponds with edge b/n points 2 & 3
		float a = ((p2.y - p3.y) * (x - p3.x) +	(p3.x - p2.x) * (y - p3.y)) / denominator;
		//b corresponds with edge b/n points 1 & 3
		float b = ((p3.y - p1.y) * (x - p3.x) + (p1.x - p3.x) * (y - p3.y)) / denominator;
		//c corresponds with edge b/n points 1 & 2
		float c = 1 - a - b;

		bool crossEdgeIndex1 = a < 0;
		bool crossEdgeIndex2 = b < 0;
		bool crossEdgeIndex0 = c < 0;

		if (crossEdgeIndex1 && crossEdgeIndex2)
			result = rand() % 2 + 1; //index = 1 or 2
		else if (crossEdgeIndex1 && crossEdgeIndex0)
			result = rand() % 2; //index = 0 or 1
		else if (crossEdgeIndex0 && crossEdgeIndex2)
			result = (rand() % 2 + 2) % 3; //index = 2 or 0
		else if (crossEdgeIndex0)
			result = 0;
		else if (crossEdgeIndex1)
			result = 1;
		else if (crossEdgeIndex2)
			result = 2;
		else
			result = -1; //inside triangle

		//SEARCH FINISHED
		if(result != -1)
			currentTri = currentTri->neighbors[result];

		if(currentTri == NULL)
			result = -1;

	}while(result != -1);

	return currentTri;
}

void Triangulator::insert(Site* site)
{
	Triangle* parentTri = getTri(site->p.x, site->p.y);

	if (parentTri == NULL) return;

	Edge* enclosingEdges[3];
	enclosingEdges[0] = parentTri->edges[0];
	enclosingEdges[1] = parentTri->edges[1];
	enclosingEdges[2] = parentTri->edges[2];
	insertSite(parentTri, site);

	for (int i = 0; i < 3; i++) {
		flip(enclosingEdges[i]->tri1, enclosingEdges[i]->tri2, enclosingEdges[i]);
	}
}

void Triangulator::flip(Triangle* tri1, Triangle* tri2, Edge* edgeInCommon)
{
	if (tri1 == NULL || tri2 == NULL || tri1 == tri2) 
		return;

	//test that they have edge in common
	bool edgeInTri1 = false;
	bool edgeInTri2 = false;
	for(int i = 0; i < 3; i++) {
		if(tri1->edges[i] == edgeInCommon) edgeInTri1 = true;
		if(tri2->edges[i] == edgeInCommon) edgeInTri2 = true;
	}

	assert(edgeInTri1 && edgeInTri2);

	Site* tri1OppositeSite = NULL;
	Site* tri2OppositeSite = NULL;

	for(int i = 0; i < 3; i++) {
		if (tri1->sites[i] != edgeInCommon->site1 && tri1->sites[i] != edgeInCommon->site2)
			tri1OppositeSite = tri1->sites[i];

		if (tri2->sites[i] != edgeInCommon->site1 && tri2->sites[i] != edgeInCommon->site2)
			tri2OppositeSite = tri2->sites[i];
	}

	Circle circleTest1 = Circle(tri1);
	Circle circleTest2 = Circle(tri2);
	float distance1 = (circleTest1.center - tri2OppositeSite->p).length();
	float distance2 = (circleTest2.center - tri1OppositeSite->p).length();
	bool inCircleTest1 = distance1 <= circleTest1.radius;
	bool inCircleTest2 = distance2 <= circleTest2.radius;

	if(!inCircleTest1 || !inCircleTest2)
		return;

	//if we're here, then we need to flip edges

	Site* oldSite1 = edgeInCommon->site1;
	Site* oldSite2 = edgeInCommon->site2;
	//remove neighbor info from old sites
	for (int i = 0; i < oldSite1->neighbors.size(); i++)
		if (oldSite1->neighbors.at(i) == oldSite2)
		{
			oldSite1->neighbors.erase(oldSite1->neighbors.begin() + i); 
			break;
		}

	for (int i = 0; i < oldSite2->neighbors.size(); i++)
		if (oldSite2->neighbors.at(i) == oldSite1)
		{
			oldSite2->neighbors.erase(oldSite2->neighbors.begin() + i); 
			break;
		}

	//add info to now connected sites
	tri1OppositeSite->neighbors.push_back(tri2OppositeSite);
	tri2OppositeSite->neighbors.push_back(tri1OppositeSite);
	//and flip the edge
	edgeInCommon->site1 = tri1OppositeSite;
	edgeInCommon->site2 = tri2OppositeSite;


	//well now the triangles are open, close them back up by switching edges
	auto selectTradeEdge = [edgeInCommon] (Triangle* tri, Point compareVec, 
		int* tradeIndex, Site* oppositeSite) -> Edge* {
		for (int i = 0; i < 3; i++)	{
			if (tri->edges[i] == edgeInCommon) continue;

			Edge* e = tri->edges[i];

			//assert(e->site1 == tri1OppositeSite || e->site2 == tri1OppositeSite);

			Site* nonOppositeSite = NULL;
			Point thisEdgeVec = Point(0, 0);
			if (e->site1 == oppositeSite)
				nonOppositeSite = e->site2;
			else
				nonOppositeSite = e->site1;

			thisEdgeVec = nonOppositeSite->p - oppositeSite->p;

			float crossTest = compareVec.crossProd(thisEdgeVec);

			if (crossTest < 0) {
				//if this test was < 0 then the remaining edge that is not the flipped edge is the 
				//one we want to keep
				*tradeIndex = i;

				for (int j = 0; j < nonOppositeSite->incidentTris.size(); j++) {
					if (nonOppositeSite->incidentTris.at(j) == tri) {
						nonOppositeSite->incidentTris.erase(nonOppositeSite->incidentTris.begin() + j);
						break;
					}
				}
				return e;
			}
		}
	};


	int tri1TradeIndex = -1;
	Point tri1CompareVec = tri2OppositeSite->p - tri1OppositeSite->p;
	Edge* tri1TradeEdge = selectTradeEdge(tri1, tri1CompareVec, &tri1TradeIndex, tri1OppositeSite);

	int tri2TradeIndex = -1;
	Point tri2CompareVec = tri1OppositeSite->p - tri2OppositeSite->p;
	Edge* tri2TradeEdge = selectTradeEdge(tri2, tri2CompareVec, &tri2TradeIndex, tri2OppositeSite);

	//trade info
	tri1->edges[tri1TradeIndex] = tri2TradeEdge;
	tri2->edges[tri2TradeIndex] = tri1TradeEdge;
	tri1OppositeSite->incidentTris.push_back(tri2);
	tri2OppositeSite->incidentTris.push_back(tri1);

	Triangle* triAcrossTri1TradeEdge = NULL;
	if (tri1TradeEdge->tri1 == tri1) {
		tri1TradeEdge->tri1 = tri2;
		triAcrossTri1TradeEdge = tri1TradeEdge->tri2;
	}
	else {
		tri1TradeEdge->tri2 = tri2;
		triAcrossTri1TradeEdge = tri1TradeEdge->tri1;
	}

	Triangle* triAcrossTri2TradeEdge = NULL;
	if (tri2TradeEdge->tri1 == tri2) {
		tri2TradeEdge->tri1 = tri1;
		triAcrossTri2TradeEdge = tri2TradeEdge->tri2;
	}
	else {
		tri2TradeEdge->tri2 = tri1;
		triAcrossTri2TradeEdge = tri2TradeEdge->tri1;
	}

	//rebuild site info
	tri1->sites[0] = tri2OppositeSite;
	tri1->sites[1] = tri1OppositeSite;
	if (tri2TradeEdge->site1 == tri2OppositeSite)
		tri1->sites[2] = tri2TradeEdge->site2;
	else
		tri1->sites[2] = tri2TradeEdge->site1;

	tri2->sites[0] = tri1OppositeSite;
	tri2->sites[1] = tri2OppositeSite;
	if (tri1TradeEdge->site1 == tri1OppositeSite)
		tri2->sites[2] = tri1TradeEdge->site2;
	else
		tri2->sites[2] = tri1TradeEdge->site1;

	//tris hold correct site pointers & in CCW order, their indexes don't nesscessarily correspond with edges
	//the way they are supposed to
	bool toggle = false;
	while (!checkEdgeCorrespondence(tri1)) {
		if (toggle) {
			Site* siteHolder = tri1->sites[1];
			tri1->sites[1] = tri1->sites[0];
			tri1->sites[0] = siteHolder;
		}
		else {
			Site* siteHolder = tri1->sites[1];
			tri1->sites[1] = tri1->sites[2];
			tri1->sites[2] = siteHolder;
		}
		toggle = !toggle;
	}
	//do the same for tri2
	while (!checkEdgeCorrespondence(tri2)) {
		if (toggle) {
			Site* siteHolder = tri2->sites[1];
			tri2->sites[1] = tri2->sites[0];
			tri2->sites[0] = siteHolder;
		}
		else {
			Site* siteHolder = tri2->sites[1];
			tri2->sites[1] = tri2->sites[2];
			tri2->sites[2] = siteHolder;
		}
		toggle = !toggle;
	}

	//and redo neighbor info
	for (int i = 0; i < 3; i++) {
		if (tri1->edges[i]->tri1 == tri1)
			tri1->neighbors[i] = tri1->edges[i]->tri2;
		else
			tri1->neighbors[i] = tri1->edges[i]->tri1;

		if (tri2->edges[i]->tri1 == tri2)
			tri2->neighbors[i] = tri2->edges[i]->tri2;
		else
			tri2->neighbors[i] = tri2->edges[i]->tri1;

		if (triAcrossTri1TradeEdge != NULL && triAcrossTri1TradeEdge->neighbors[i] == tri1)
			triAcrossTri1TradeEdge->neighbors[i] = tri2;

		if (triAcrossTri2TradeEdge != NULL && triAcrossTri2TradeEdge->neighbors[i] == tri2)
			triAcrossTri2TradeEdge->neighbors[i] = tri1;
	}

	//and do recursive calls
	for (int i = 0; i < 3; i++)	{
		Edge* e;
		e = tri1->edges[i];

		if (e != edgeInCommon && e->tri1 != NULL && e->tri2 != NULL)
			flip(e->tri1, e->tri2, e);

		e = tri2->edges[i];

		if (e != edgeInCommon && e->tri1 != NULL && e->tri2 != NULL)
			flip(e->tri1, e->tri2, e);
	}
}

void Triangulator::insertSite(Triangle* enclosingTri, Site* newSite)
{
	Site* p0 = enclosingTri->sites[0];
	Site* p1 = enclosingTri->sites[1];
	Site* p2 = enclosingTri->sites[2];

	sites.push_back(newSite);

	Edge* enclosingEdges[3];

	for (int i = 0; i < 3; i++)	{
		Edge* e = enclosingTri->edges[i];
		if ((e->site1 == p0 && e->site2 == p1) || (e->site1 == p1 && e->site2 == p0))
			enclosingEdges[0] = e;
		else if ((e->site1 == p1 && e->site2 == p2) || (e->site1 == p2 && e->site2 == p1))
			enclosingEdges[1] = e;
		else //assumed that if it didn't match last two criteria, that this edge is p2 & p0
			enclosingEdges[2] = e;
	}
	//now old edges is organized such that: oldEdges[0] = p0 & p1, oldEdges[1] = p1 & p2, oldEdges[2] = p0 & p2

	Edge* newEdge1 = new Edge(newSite, p0);
	Edge* newEdge2 = new Edge(newSite, p1);
	Edge* newEdge3 = new Edge(newSite, p2);
	edges.push_back(newEdge1);
	edges.push_back(newEdge2);
	edges.push_back(newEdge3);

	Triangle* newTri1 = new Triangle();
	Triangle* newTri2 = new Triangle();
	Triangle* newTri3 = new Triangle();

	//register sites within triangles
	newTri1->sites[0] = p0;
	newTri1->sites[1] = p1;
	newTri1->sites[2] = newSite;
	newTri2->sites[0] = p1;
	newTri2->sites[1] = p2;
	newTri2->sites[2] = newSite;
	newTri3->sites[0] = p2;
	newTri3->sites[1] = p0;
	newTri3->sites[2] = newSite;

	//register edges within tri
	//edges are registered such that edges[0] corresponds with sites[0] & sites[1],
	//edges[1] corresponds with sites[1] & sites[2], and edges[2] corresponds with sites[2] & sites[0]
	//this is important in the "getTri" function, which walks across the triangulation
	newTri1->edges[0] = enclosingEdges[0];
	newTri1->edges[1] = newEdge2;
	newTri1->edges[2] = newEdge1;

	newTri2->edges[0] = enclosingEdges[1];
	newTri2->edges[1] = newEdge3;
	newTri2->edges[2] = newEdge2;

	newTri3->edges[0] = enclosingEdges[2];
	newTri3->edges[1] = newEdge1;
	newTri3->edges[2] = newEdge3;

	newEdge1->tri1 = newTri1;
	newEdge1->tri2 = newTri3;
	newEdge2->tri1 = newTri1;
	newEdge2->tri2 = newTri2;
	newEdge3->tri1 = newTri2;
	newEdge3->tri2 = newTri3;

	//register neighbors
	newTri1->neighbors[1] = newTri2;
	newTri2->neighbors[2] = newTri1;
	newTri1->neighbors[2] = newTri3;
	newTri3->neighbors[1] = newTri1;
	newTri2->neighbors[1] = newTri3;
	newTri3->neighbors[2] = newTri2;

	if (enclosingEdges[0]->tri1 == enclosingTri) {
		enclosingEdges[0]->tri1 = newTri1;
		newTri1->neighbors[0] = enclosingEdges[0]->tri2;
	}
	else {
		enclosingEdges[0]->tri2 = newTri1;
		newTri1->neighbors[0] = enclosingEdges[0]->tri1;
	}

	if (enclosingEdges[1]->tri1 == enclosingTri) {
		enclosingEdges[1]->tri1 = newTri2;
		newTri2->neighbors[0] = enclosingEdges[1]->tri2;
	}
	else {
		enclosingEdges[1]->tri2 = newTri2;
		newTri2->neighbors[0] = enclosingEdges[1]->tri1;
	}

	if (enclosingEdges[2]->tri1 == enclosingTri) {
		enclosingEdges[2]->tri1 = newTri3;
		newTri3->neighbors[0] = enclosingEdges[2]->tri2;
	}
	else {
		enclosingEdges[2]->tri2 = newTri3;
		newTri3->neighbors[0] = enclosingEdges[2]->tri1;
	}

	//add reference of newtris to neighbor tris
	Triangle* newTris[3]; newTris[0] = newTri1; newTris[1] = newTri2; newTris[2] = newTri3;
	for (int i = 0; i < 3; i++) {
		Triangle* outsideTri = newTris[i]->neighbors[0];

		if (outsideTri == NULL)
			continue;

		Edge* e = newTris[i]->edges[0];
		for (int j = 0; j < 3; j++)	{
			if (outsideTri->edges[j] == e)	{
				outsideTri->neighbors[j] = newTris[i];
				break;
			}
		}
	}

	for (int i = 0; i < p0->incidentTris.size(); i++) {
		if (p0->incidentTris.at(i) == enclosingTri)	{
			p0->incidentTris.erase(p0->incidentTris.begin() + i);
			break;
		}
	}
	p0->incidentTris.push_back(newTri1);
	p0->incidentTris.push_back(newTri3);

	for (int i = 0; i < p1->incidentTris.size(); i++) {
		if (p1->incidentTris.at(i) == enclosingTri) {
			p1->incidentTris.erase(p1->incidentTris.begin() + i);
			break;
		}
	}
	p1->incidentTris.push_back(newTri1);
	p1->incidentTris.push_back(newTri2);

	for (int i = 0; i < p2->incidentTris.size(); i++) {
		if (p2->incidentTris.at(i) == enclosingTri)	{
			p2->incidentTris.erase(p2->incidentTris.begin() + i);
			break;
		}
	}
	p2->incidentTris.push_back(newTri2);
	p2->incidentTris.push_back(newTri3);

	for (int i = 0; i < tris.size(); i++) {
		if (tris.at(i) == enclosingTri)	{
			tris.erase(tris.begin() + i);
			break;
		}
	}
	delete enclosingTri;

	newSite->incidentTris.push_back(newTri1);
	newSite->incidentTris.push_back(newTri2);
	newSite->incidentTris.push_back(newTri3);
	newSite->neighbors.push_back(p0);
	newSite->neighbors.push_back(p1);
	newSite->neighbors.push_back(p2);
	p0->neighbors.push_back(newSite);
	p1->neighbors.push_back(newSite);
	p2->neighbors.push_back(newSite);

	tris.push_back(newTri1);
	tris.push_back(newTri2);
	tris.push_back(newTri3);
}

bool Triangulator::checkEdgeCorrespondence(Triangle* tri)
{
	bool e1Check = (tri->edges[0]->site1 == tri->sites[0] && tri->edges[0]->site2 == tri->sites[1])
		|| (tri->edges[0]->site1 == tri->sites[1] && tri->edges[0]->site2 == tri->sites[0]);

	bool e2Check = (tri->edges[1]->site1 == tri->sites[1] && tri->edges[1]->site2 == tri->sites[2])
		|| (tri->edges[1]->site1 == tri->sites[2] && tri->edges[1]->site2 == tri->sites[1]);

	bool e3Check = (tri->edges[2]->site1 == tri->sites[2] && tri->edges[2]->site2 == tri->sites[0])
		|| (tri->edges[2]->site1 == tri->sites[0] && tri->edges[2]->site2 == tri->sites[2]);

	return e1Check && e2Check && e3Check;
}

void Triangulator::floodFill(Site* startSite, float cutoff, std::vector<Site*>* outlist) const
{
	if (startSite == NULL) return;

	std::vector<Site*> notExaustedSites;

	notExaustedSites.push_back(startSite);
	outlist->push_back(startSite);

	while (notExaustedSites.size() != 0)
	{
		Site* frontSite = notExaustedSites.front();

		for (int i = 0; i < frontSite->neighbors.size(); i++)
		{
			Site* neighbor = frontSite->neighbors.at(i);
			
			if ((startSite->p - neighbor->p).length() > cutoff)
				continue;

			bool inOutList = std::find(outlist->begin(), outlist->end(), neighbor) != outlist->end();

			if (!inOutList)
				outlist->push_back(neighbor);

			bool inNonExausted = std::find(notExaustedSites.begin(), notExaustedSites.end(), neighbor)
				!= notExaustedSites.end();

			if (!inOutList && !inNonExausted)
				notExaustedSites.push_back(neighbor);

		}

		notExaustedSites.erase(notExaustedSites.begin());
	};
}

Site* Triangulator::getNearestSite(float x, float y) const
{
	Triangle* tri = getTri(x, y);

	if (tri == NULL)
		return NULL;

	Point p = Point(x, y);

	int nearestIndex = 0;
	float closestDist = (p - tri->sites[0]->p).length();

	float dist2 = (p - tri->sites[1]->p).length();
	if (dist2 < closestDist) {
		closestDist = dist2;
		nearestIndex = 1;
	}

	float dist3 = (p - tri->sites[2]->p).length();
	if (dist3 < closestDist)
		nearestIndex = 2;

	return tri->sites[nearestIndex];
}

Circle::Circle(Triangle* triangle)
{	
	Point p0 = triangle->sites[0]->p;
	Point p1 = triangle->sites[1]->p;
	Point p2 = triangle->sites[2]->p;

	Point midPoint01 = (p0 + p1) * 0.5;
	Point midPoint12 = (p1 + p2) * 0.5;

	Point bisector1 = p0 - p1;
	Point bisector2 = p2 - p1;
	bisector1 = Point(-bisector1.y, bisector1.x);
	bisector2 = Point(-bisector2.y, bisector2.x);

	auto Cross = [](Point p1, Point p2) -> float {
		return (p1.x * p2.y) - (p1.y * p2.x);
	};

	float t = Cross((midPoint01 - midPoint12), bisector1) / Cross(bisector2, bisector1);
	Point inter = bisector2;
	inter.x *= t;
	inter.y *= t;

	center = midPoint12 + inter;
	radius = (center - p0).length();
};

VoronoiGraph::VoronoiGraph(Triangulator* triangulator)
{
	nodeCount = triangulator->tris.size();
	nodes = (VoronoiNode*)malloc(sizeof(VoronoiNode) * nodeCount);

	std::unordered_map<Triangle*, int> triangleIndexMap;

	for (int i = 0; i < nodeCount; i++) {
		Triangle* triPointer = triangulator->tris.at(i);
		nodes[i] = VoronoiNode(triPointer);
		//why does this matter depending on project (VS 2010 vs. VS 2012?)
		triangleIndexMap.insert(std::make_pair(triPointer, i));
		//triangleIndexMap.emplace(triPointer, i);
	}

	for (int i = 0; i < nodeCount; i++) {
		Triangle* baseTri = triangulator->tris.at(i);
		VoronoiNode* node = &nodes[i];

		for (int j = 0; j < 3; j++) {
			Triangle* neighbor = baseTri->neighbors[j];
			if (neighbor == NULL) continue;

			int n_index = triangleIndexMap.find(neighbor)->second;
			node->neighbors[j] = &nodes[n_index];
		}
	}
}

VoronoiGraph::~VoronoiGraph()
{
	free(nodes);
}