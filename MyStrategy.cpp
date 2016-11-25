#include "MyStrategy.h"
#include <functional>
#define PI 3.14159265358979323846
#define _USE_MATH_DEFINES

#include <chrono>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <queue>
#include <iostream>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <stack>
#include <list>

#ifdef _zuko3d_pc
//includes
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

//data structures
struct Colour {
	unsigned char r, g, b, a;
};

class TGAImage {

public:

	//Constructor
	TGAImage();

	//Overridden Constructor
	TGAImage(short width, short height);

	//Set all pixels at once
	void setAllPixels(Colour *pixels);

	//set individual pixels
	void setPixel(Colour inputcolor, int xposition, int yposition);

	void WriteImage(string filename);

	//General getters and setters

	void setWidth(short width);
	void setHeight(short height);

	short getWidth();
	short getHeight();

private:

	//store the pixels
	Colour *m_pixels;

	short m_height;
	short m_width;

	//convert 2D to 1D indexing
	int convert2dto1d(int x, int y);



};

//Default Constructor
TGAImage::TGAImage() {

}

//Overridden Constructor
TGAImage::TGAImage(short width, short height) {
	m_width = width;
	m_height = height;
	m_pixels = new Colour[m_width*m_height];
}

//Set all pixels at once
void TGAImage::setAllPixels(Colour *pixels) {
	m_pixels = pixels;
}

//Set indivdual pixels
void TGAImage::setPixel(Colour inputcolor, int x, int y) {
	m_pixels[convert2dto1d(x, y)] = inputcolor;
}

//Convert 2d array indexing to 1d indexing
int TGAImage::convert2dto1d(int x, int y) {
	return m_width * (m_height - y) + x;
}

void TGAImage::WriteImage(string filename) {

	//Error checking
	if (m_width <= 0 || m_height <= 0)
	{
		cout << "Image size is not set properly" << endl;
		return;
	}

	ofstream o(filename.c_str(), ios::out | ios::binary);

	//Write the header
	o.put(0);
	o.put(0);
	o.put(2);                         /* uncompressed RGB */
	o.put(0); 		o.put(0);
	o.put(0); 	o.put(0);
	o.put(0);
	o.put(0); 	o.put(0);           /* X origin */
	o.put(0); 	o.put(0);           /* y origin */
	o.put((m_width & 0x00FF));
	o.put((m_width & 0xFF00) / 256);
	o.put((m_height & 0x00FF));
	o.put((m_height & 0xFF00) / 256);
	o.put(32);                        /* 24 bit bitmap */
	o.put(0);

	//Write the pixel data
	for (int i = 0; i<m_height*m_width; i++) {
		o.put(m_pixels[i].b);
		o.put(m_pixels[i].g);
		o.put(m_pixels[i].r);
		o.put(m_pixels[i].a);
	}

	//close the file
	o.close();

}

#endif

using namespace model;
using namespace std;

model::Faction enFaction, myFaction;

long long time() {
	return chrono::high_resolution_clock::now().time_since_epoch().count();
}

template <class c> void fast_remove(vector<c>& v, const c& op) {
	auto it = &v[0];
	auto it_end = it + v.size();
	bool removed = false;
	while (it != it_end) {
		if (*it == op) {
			*it = v.back();
			v.pop_back();
			break;
		}
		it++;
	}
}

template <class T> inline T sqr(T x) {
	return x * x;
}

namespace std {
	template <>
	struct hash<pair<int, int>>
	{
		std::size_t operator()(const pair<int, int>& k) const
		{
			return ((k.first) << 10) | (k.second); // мы знаем, что этот хэш исопльзуется для клеток, и их значения вполне предсказуемы.
		}
	};
}
struct circle {
	double radius = 0.0;
	double x = 0.0, y = 0.0;
	long long id;

	circle() { }

	circle(const CircularUnit& u, double additionalRadius = 35) : radius(u.getRadius() + additionalRadius), x(u.getX()), y(u.getY()), id(u.getId()) { }

	bool operator<(const circle& c) const {
		return id < c.id;
	}

	bool operator==(const circle& c) const {
		return id == c.id;
	}

	double distanceTo(double px, double py) const {
		return sqrt(sqr(x - px) + sqr(y - py));
	}
};

struct dcircle : circle {
	double sightRange = 0.0;
	int type;
	double attackRange; // максимальное расстояние между центрами, при котором можно нанести урон
	double hpPct;
	double hp;

	dcircle() { }

	dcircle(const LivingUnit& u, double additionalRadius = 38) : circle(u, additionalRadius), hpPct((double) u.getLife() / (double) u.getMaxLife()), hp((double) u.getLife()) { }

	dcircle(const Minion& u, double additionalRadius = 38) : dcircle((LivingUnit) u, additionalRadius) {
		sightRange = u.getVisionRange();
		type = 2; // melee
		attackRange = 50 + 35;
		if (u.getType() == MINION_FETISH_BLOWDART) {
			type = 1;
			attackRange = 305 + 35;
		}
	}
	dcircle(const Wizard& u, double additionalRadius = 38) : dcircle((LivingUnit)u, additionalRadius) {
		sightRange = u.getVisionRange();
		type = 3;
		attackRange = 510 + 35;
	}
	dcircle(const Building& u, double additionalRadius = 38) : dcircle((LivingUnit)u, additionalRadius) {
		sightRange = u.getVisionRange();
		type = 4;
		attackRange = 600;
		if (u.getType() == BUILDING_FACTION_BASE) {
			type = 5;
			attackRange = 800;
		}
	}

	bool operator<(const dcircle& c) const {
		return id < c.id;
	}

	bool operator==(const dcircle& c) const {
		return id == c.id;
	}
};

struct circleSetComarator {
	bool operator() (const circle & c1, const circle& c2) {
		return c1.id < c2.id;
	}
};

double cellSize = 150.0; // не меньше, чем радиус самого большого объекта (85 при данных правилах, не считая базы)
set<circle> prevStaticObjects;

vector< vector< vector< circle > > > staticMap;
vector< vector< vector< dcircle > > > dynamicMap;

class point : public Unit {
public:
	double x, y;
	point(double a, double b) : Unit(0, a, b, 0,0,0,FACTION_ACADEMY) {
		x = a;
		y = b;
	}

	point(const dcircle& c) : Unit(0, c.x, c.y, 0, 0, 0, FACTION_ACADEMY) {
		x = c.x;
		y = c.y;
	}
};

double maxMSX(const Wizard& w) {
	return 3.0;// !!!!!!!!!!!!!!!!!!! НУЖНО УЧИТЫВАТЬ БАФФЫ на скорость
}

double maxMSY(const Wizard& w) {
	return 4.0;// !!!!!!!!!!!!!!!!!!! НУЖНО УЧИТЫВАТЬ БАФФЫ на скорость
}

bool checkCollision(const point& a, const point& b, const circle& c) {
	if (sqr(a.x - c.x) + sqr(a.y - c.y) < sqr(c.radius)) return true;
	if (sqr(b.x - c.x) + sqr(b.y - c.y) < sqr(c.radius)) return true;

	double abx = b.x - a.x;
	double aby = b.y - a.y;

	double bcx = c.x - b.x;
	double bcy = c.y - b.y;

	if ((abx * bcx + aby * bcy) > 0) return false;

	double acx = c.x - a.x;
	double acy = c.y - a.y;

	if ((abx * acx + aby * acy) < 0) return false;

	double d = sqrt(abx * abx + aby * aby);
	double A = -aby / d;
	double B = abx / d;
	double C = -(A * a.x + B * a.y);

	if (fabs(A * c.x + B * c.y + C) < c.radius + 1e-2) return true;
	return false;
}

bool checkCollision(const point& a, const point& b, const CircularUnit& c) {
	if (a.getDistanceTo(c) < c.getRadius()) return true;
	if (b.getDistanceTo(c) < c.getRadius()) return true;

	double abx = b.getX() - a.getX();
	double aby = b.getY() - a.getY();

	double bcx = c.getX() - b.getX();
	double bcy = c.getY() - b.getY();

	if ((abx * bcx + aby * bcy) > 0) return false;

	double acx = c.getX() - a.getX();
	double acy = c.getY() - a.getY();

	if ((abx * acx + aby * acy) < 0) return false;

	double d = sqrt(abx * abx + aby * aby);
	double A = -aby / d;
	double B = abx / d;
	double C = -(A * a.getX() + B * a.getY());

	if (fabs(A * c.getX() + B * c.getY() + C) < c.getRadius() + 1e-2) return true;
	return false;
}

/*
template <class c1, class c2> class default_map : map<c1, c2> {
	c2 default_value;

public:
	default_map(const c2& val) : default_value(val) { }

	c2 operator[](const c1& key) const {
		if(this->find(key) == this->end()) return default_value;
		return ((map<c1,c2>) *this).operator[](key);
	}

	c2& operator[](const c1& key) {
		return ((map<c1, c2>) *this).operator[](key);
	}
};
*/

template <class circle_class> bool checkCollisionsShort(double ax, double ay, double bx, double by, const vector< vector< vector<circle_class> > >& objects) {
	typedef pair<int, int> pt;
	set<pt> cells;
	
	pt a_pt((int)floor(ax / cellSize), (int)floor(ay / cellSize));
	pt b_pt((int)floor(bx / cellSize), (int)floor(by / cellSize));
	if ((a_pt.first >= 0) && (a_pt.first < objects.size())) {
		if ((a_pt.second >= 0) && (a_pt.second < objects.size())) {
			cells.insert(a_pt);
		}
	}
	if (a_pt != b_pt) {
		if ((b_pt.first >= 0) && (b_pt.first < objects.size())) {
			if ((b_pt.second >= 0) && (b_pt.second < objects.size())) {
				cells.insert(b_pt);
			}
		}

		if ((a_pt.first != b_pt.first) && (a_pt.second != b_pt.second)) {
			double A, B, C;
			A = ay - by;
			B = bx - ax;
			C = -(A * ax + B * ay);

			if (a_pt > b_pt) {
				swap(a_pt, b_pt);
				swap(ax, bx);
				swap(ay, by);
			}
			
			double crossx = floor(bx / cellSize) * cellSize;
			double crossy = -(C + A * crossx) / B;
			int sec = (int)floor(crossy / cellSize);
			if ((b_pt.first >= 0) && (b_pt.first < objects.size())) {
				if ((sec >= 0) && (sec < objects.size())) {
					cells.insert(pt(b_pt.first, sec));
				}
			}
		}
	}


	// ----------------------------------------------

	set<circle_class> circles;

	point a(ax, ay), b(bx, by);
	for (auto&c : cells) {
		for (auto& obj : objects[c.first][c.second]) {
			circles.insert(obj);
		}
	}

	for (const auto& obj : circles) {
		if (checkCollision(a, b, obj)) return true;
	}

	return false;
}

template <class circle_class> bool checkCollisions(double ax, double ay, const vector< vector< vector<circle_class> > >& objects) {
	if ((ax < 0.0) || (ax > 4000.0)) return false;
	if ((ay < 0.0) || (ay > 4000.0)) return false;
	
	int xi = (int)floor(ax / cellSize);
	int yi = (int)floor(ay / cellSize);

	for (auto& obj : objects[xi][yi]) {
		if (obj.distanceTo(ax, ay) < obj.radius) return true;
	}
	
	return false;
}


template<class circle_class> vector<circle_class>* cellByXY(double x, double y, vector< vector< vector<circle_class> > > &objects) {
	return &(objects[(int)floor(x / cellSize)][(int)floor(y / cellSize)]);
}

bool checkCollisions(double ax, double ay, double bx, double by) {
	typedef pair<int, int> pt;
	set<pt> cells;
	double A, B, C;
	A = ay - by;
	B = bx - ax;
	C = -(A * ax + B * ay);
	
	/*
	if (fabs(A * bx + B * by + C) > 1e-8) {
		cout << "ERROR in line formula" << endl;
	}
	*/

	double sx = floor(ax / cellSize);
	double sy = floor(ay / cellSize);

	int i, j;
	i = (int)sx;

	cells.insert(pt((int)sx, (int)sy));
	cells.insert(pt((int)floor(bx / cellSize), (int)floor(by / cellSize)));
	sx *= cellSize;
	sy *= cellSize;

	double dx = cellSize * (bx > ax ? 1 : -1);
	int di = (bx > ax ? 1 : -1);
	sx += dx;
	i+=di;
	while (sx * di < bx * di) {
		double y = -(A * sx + C) / B;
		j = (int)floor(y / cellSize);
		cells.insert(pt(i, j));
		cells.insert(pt(i - 1, j));
		i += di;
		sx += dx;
		if (i < 0) {
			cout << "ERROR i < 0!" << endl;
		}
		if (j < 0) {
			cout << "ERROR j < 0!" << endl;
		}
	}

	double dy = cellSize * (by > ay ? 1 : -1);
	int dj = (by > ay ? 1 : -1);
	j = (int) (sy / cellSize);
	sy += dy;
	j += dj;
	while (sy * dj < by * dj) {
		double x = -(B * sy + C) / A;
		i = (int)floor(x / cellSize);
		cells.insert(pt(i, j));
		cells.insert(pt(i, j - 1));
		j += dj;
		sy += dy;
		if (i < 0) {
			cout << "ERROR i < 0!" << endl;
		}
		if (j < 0) {
			cout << "ERROR j < 0!" << endl;
		}
	}

	point a(ax, ay), b(bx, by);
	for (const auto&c : cells) {
		for (const auto& obj : staticMap[c.first][c.second]) {
			if (checkCollision(a, b, obj)) return true;
		}
		for (const auto& obj : dynamicMap[c.first][c.second]) {
			if (checkCollision(a, b, obj)) return true;
		}
	}

	return false;
}

bool checkCollisions(const point& a, const point &b) {
	return checkCollisions(a.getX(), a.getY(), b.getX(), b.getY());
}

void smoothenPath(list<point>& path) {
	list<point>::iterator itl, itr, it_end;
	itl = path.begin();

	int init_size = path.size();

	it_end = path.end();
	auto itl_next = itl;
	itr = it_end;
	itr--;
	itl_next++;
	while (itr != itl_next) {
		if (!checkCollisions(*itl, *itr)) {
			path.erase(itl_next, itr);
			break;
		}
		itr--;
	}
	path.pop_front();
	cout << "path start-optimized for \t" << init_size << " -> \t" << path.size() << endl;
}

void smoothenPathFull(list<point>& path) {
	list<point>::iterator itl, itr, it_end;
	itl = path.begin();
	
	int init_size = path.size();

	it_end = path.end();
	it_end--;
	auto itl_next = itl;
	while (itl != it_end) {
		itr = it_end;
		itl_next = itl;
		itl_next++;
		while (itr != itl_next) {
			if (!checkCollisions(*itl, *itr)) {
				auto ti1 = itl;
				ti1++;
				path.erase(ti1, itr);
				break;
			}

			itr--;
		}

		itl++;
	}

	cout << "path full-optimized for \t" << init_size << " -> \t" << path.size() << endl;
}

long path_ticks = 0;

auto findPath(const point& start, const point& dest, double step = 20.0) {
	typedef pair<int, int> pt;
	typedef pair<double, pt > qtype;

	auto h = [](const pt&a, const pt&b) -> double {return sqrt((double) (sqr(a.first - b.first) + sqr(a.second - b.second))); };

	list<point> path;
		
	unordered_set<pt> closed, openSet;
	priority_queue < qtype, vector<qtype>, std::greater<qtype> > open;
	unordered_map<pt, pt> cameFrom;
	
	unordered_map<pt, double> gscore;

	// !!!!!!!!!!! goal может быть недостижима. Надо как минимум осматривать несколько соседних точек.
	pt pstart(0, 0), goal((int) floor((dest.getX() - start.getX()) / step), (int) floor((dest.getY() - start.getY())/ step)); // !!!!!!!!!!! goal может быть недостижима. Надо как минимум осматривать несколько соседних точек.
	
	gscore[pstart] = 0.0;

	open.push(qtype(h(pstart, goal), pt(0, 0)));
	openSet.insert(pt(0, 0));

	bool good = false;
	qtype cur;
	long cur_path_ticks = 0;
	while (!open.empty()) {
		cur_path_ticks++;
		if (cur_path_ticks > 25000) break;
		cur = open.top();
		while (closed.find(cur.second) != closed.end()) {
			open.pop();
			cur = open.top();
		}
		if (cur.second == goal) {
			good = true;
			break;
		}

		open.pop();
		openSet.erase(cur.second);
		closed.insert(cur.second);
		double curx = step * cur.second.first + start.getX();
		double cury = step * cur.second.second + start.getY();
		for (int i = -1; i < 2; i++) {
			for (int j = -1; j < 2; j++) {
				if ((0 == j) && (i == 0)) continue;
				
				auto neib = pt(cur.second.first + i, cur.second.second + j);

				double nx = step * neib.first + start.getX();
				double ny = step * neib.second + start.getY();
				
				int cellx = (int) floor(nx / cellSize);
				int celly = (int) floor(ny / cellSize);
				if (cellx < 0) {
					closed.insert(neib);
					continue;
				}
				if (celly < 0) {
					closed.insert(neib);
					continue;
				}

				bool inObstalce = false;
				inObstalce = checkCollisionsShort(curx, cury, nx, ny, staticMap);
				if(!inObstalce) inObstalce = checkCollisionsShort(curx, cury, nx, ny, dynamicMap);
				//inObstalce = checkCollisions(curx, cury, nx, ny, staticMap);
				/*
				for (auto &c : staticMap[cellx][celly]) {
					if (sqr(c.radius) > sqr(c.x - nx) + sqr(c.y - ny)) {
						inObstalce = true;
						break;
					}
				}
				*/

				if (inObstalce) {
					closed.insert(neib);
					continue;
				}

				if (closed.find(neib) == closed.end()) {
					auto tScore = gscore[cur.second] + sqrt((double) (abs(i) + abs(j)));
					if (gscore.find(neib) != gscore.end()) {
						if (gscore[neib] <= tScore) continue;
					}
					gscore[neib] = tScore;
					cameFrom[neib] = cur.second;
					open.push(qtype(tScore + h(neib, goal), neib));
					openSet.insert(neib);
				}
			}
		}		
	}

	if (!good) {
		cout << "FAILED to find a path" << endl;
		return path;
	}

	if (cur_path_ticks > path_ticks) {
		path_ticks = cur_path_ticks;
	}

	cout << "path_ticks = \t" << path_ticks << endl;

	pt p_cur = cur.second;
	double x = start.getX();
	double y = start.getY();
	path.clear();
	path.push_front(dest);
	while (p_cur != pstart) {
		path.push_front(point(x + p_cur.first * step, y + p_cur.second * step));
		p_cur = cameFrom[p_cur];
	}
	path.push_front(start);
	if (path.size() < 3) {
		cout << "short path!" << endl;
	}
	smoothenPath(path);
	return path;
}
 
auto findPathToZone(const point& start, const point& dest, const double rad, double step = 20.0) {
	typedef pair<int, int> pt;
	typedef pair<double, pt > qtype;

	auto h = [rad](const pt&a, const pt&b) -> double {return max(0.0, sqrt((double)(sqr(a.first - b.first) + sqr(a.second - b.second))) - rad); };

	list<point> path;

	unordered_set<pt> closed, openSet;
	priority_queue < qtype, vector<qtype>, std::greater<qtype> > open;
	unordered_map<pt, pt> cameFrom;

	unordered_map<pt, double> gscore;

	// !!!!!!!!!!! goal может быть недостижима. Надо как минимум осматривать несколько соседних точек.
	pt pstart(0, 0), goal((int)floor((dest.getX() - start.getX()) / step), (int)floor((dest.getY() - start.getY()) / step)); // !!!!!!!!!!! goal может быть недостижима. Надо как минимум осматривать несколько соседних точек.

	gscore[pstart] = 0.0;

	open.push(qtype(h(pstart, goal), pt(0, 0)));
	openSet.insert(pt(0, 0));

	bool good = false;
	qtype cur;
	long cur_path_ticks = 0;
	while (!open.empty()) {
		cur_path_ticks++;
		if (cur_path_ticks > 25000) break;
		cur = open.top();
		while (closed.find(cur.second) != closed.end()) {
			open.pop();
			cur = open.top();
		}
		if (sqr(cur.second.first - goal.first) + sqr(cur.second.second - goal.second) <= sqr(rad/step)) {
			good = true;
			break;
		}

		open.pop();
		openSet.erase(cur.second);
		closed.insert(cur.second);
		double curx = step * cur.second.first + start.getX();
		double cury = step * cur.second.second + start.getY();
		for (int i = -1; i < 2; i++) {
			for (int j = -1; j < 2; j++) {
				if ((0 == j) && (i == 0)) continue;

				auto neib = pt(cur.second.first + i, cur.second.second + j);

				double nx = step * neib.first + start.getX();
				double ny = step * neib.second + start.getY();

				int cellx = (int)floor(nx / cellSize);
				int celly = (int)floor(ny / cellSize);
				if (cellx < 0) {
					closed.insert(neib);
					continue;
				}
				if (celly < 0) {
					closed.insert(neib);
					continue;
				}

				bool inObstalce = false;
				inObstalce = checkCollisionsShort(curx, cury, nx, ny, staticMap);
				if (!inObstalce) inObstalce = checkCollisionsShort(curx, cury, nx, ny, dynamicMap);
				//inObstalce = checkCollisions(curx, cury, nx, ny, staticMap);
				/*
				for (auto &c : staticMap[cellx][celly]) {
				if (sqr(c.radius) > sqr(c.x - nx) + sqr(c.y - ny)) {
				inObstalce = true;
				break;
				}
				}
				*/

				if (inObstalce) {
					closed.insert(neib);
					continue;
				}

				if (closed.find(neib) == closed.end()) {
					auto tScore = gscore[cur.second] + sqrt((double)(abs(i) + abs(j)));
					if (gscore.find(neib) != gscore.end()) {
						if (gscore[neib] <= tScore) continue;
					}
					gscore[neib] = tScore;
					cameFrom[neib] = cur.second;
					open.push(qtype(tScore + h(neib, goal), neib));
					openSet.insert(neib);
				}
			}
		}
	}

	if (!good) {
		cout << "FAILED to find a path" << endl;
		return path;
	}

	/*
	if (cur_path_ticks > path_ticks) {
		path_ticks = cur_path_ticks;
	}

	cout << "path_ticks = \t" << path_ticks << endl;
	*/

	pt p_cur = cur.second;
	double x = start.getX();
	double y = start.getY();
	path.clear();
	path.push_front(dest);
	while (p_cur != pstart) {
		path.push_front(point(x + p_cur.first * step, y + p_cur.second * step));
		p_cur = cameFrom[p_cur];
	}
	path.push_front(start);
	if (path.size() < 3) {
		cout << "short path!" << endl;
	}
	smoothenPath(path);
	return path;
}


list <point> curWay;
double posX = 0, posY = 0;
void setMoveToPoint(const Wizard& self, Move& move, const Unit& pt) {
	double angle = self.getAngleTo(pt);
	//cout << "angle = " << angle << endl;
	//cout << "dx = " << self.getX() - posX << "\t dy = " << self.getY() - posY << endl;
	posX = self.getX();
	posY = self.getY();
	//cout << self.getX() << ", " << self.getY() << endl;
	
	double maxSpeedX = maxMSX(self); // !!!!!!!!!!!!!!!!!!! НУЖНО УЧИТЫВАТЬ БАФФЫ на скорость
	double maxSpeedY = maxMSY(self); // !!!!!!!!!!!!!!!!!!! НУЖНО УЧИТЫВАТЬ БАФФЫ на скорость

	double dist = self.getDistanceTo(pt);
	if (dist <= maxSpeedX) {
		move.setStrafeSpeed(sin(angle) * dist);
		move.setSpeed(cos(angle) * dist);
	} else {
		if (fabs(angle) < 1e-17) {
			move.setSpeed(maxSpeedY);
			move.setStrafeSpeed(0.0);
		} else {
			if (fabs(angle - PI / 2.0) < 1e-17) {
				move.setSpeed(-maxSpeedX);
				move.setStrafeSpeed(0.0);
			}
			else {
				if (fabs(angle) <= PI / 2.0) {
					double dx = sqrt(1.0 / (sqr(1.0 / maxSpeedX) + sqr(1.0 / tan(angle) / maxSpeedY)));
					move.setSpeed(fabs(dx / tan(angle)));
					if (angle < 0) dx *= -1;
					move.setStrafeSpeed(dx);
				}
				else {
					move.setStrafeSpeed(sin(angle) * maxSpeedX);
					move.setSpeed(cos(angle) * maxSpeedX);
				}
			}
		}
	}
}

void setMoveToPoint(const Wizard& self, Move& move, double x, double y) {
	setMoveToPoint(self, move, point(x,y));
}

template <class circle_class> void addObjectsToMap(vector< vector< vector< circle_class > > > &mp, const set<circle_class>& objects, double map_width, double map_height ) {
	for (const auto& c : objects) {
		//cout << "new obj: " << c.x << "\t" << c.y << "\t" << c.radius << endl;

		bool top = false;
		bool bot = false;
		bool left = false;
		bool right = false;

		int i = (int) floor(c.x / map_width);
		int j = (int) floor(c.y / map_height);

		mp[i][j].push_back(c);

		if (j > (int)floor((c.y - c.radius) / map_height)) {
			top = true;
			if (j > 0) mp[i][j - 1].push_back(c);
		}

		if (j < (int)floor((c.y + c.radius) / map_height)) {
			bot = true;
			mp[i][j + 1].push_back(c);
		}

		if (i >(int)floor((c.x - c.radius) / map_width)) {
			left = true;
			if (i > 0) mp[i - 1][j].push_back(c);
		}

		if (i < (int)floor((c.x + c.radius) / map_width)) {
			right = true;
			mp[i + 1][j].push_back(c);
		}

		if (top && left) {
			double x = i * map_width;
			double y = j * map_height;
			if (i > 0) if (j > 0) if (sqr(x - c.x) + sqr(y - c.y) < c.radius) mp[i - 1][j - 1].push_back(c);
		}

		if (top && right) {
			double x = (i + 1) * map_width;
			double y = j * map_height;
			if (j > 0) if (sqr(x - c.x) + sqr(y - c.y) < c.radius) mp[i + 1][j - 1].push_back(c);
		}

		if (bot && left) {
			double x = i * map_width;
			double y = (j + 1) * map_height;
			if (i > 0) if (sqr(x - c.x) + sqr(y - c.y) < c.radius) mp[i - 1][j + 1].push_back(c);
		}

		if (bot && right) {
			double x = (i + 1) * map_width;
			double y = (j + 1) * map_height;
			if (sqr(x - c.x) + sqr(y - c.y) < c.radius) mp[i + 1][j + 1].push_back(c);
		}

	}
}

void removeObjectsFromMap(vector< vector< vector< circle > > > &mp, const set<circle>& objects, double map_width, double map_height) {
	for (const auto& c : objects) {
		//cout << "new obj: " << c.x << "\t" << c.y << "\t" << c.radius << endl;

		bool top = false;
		bool bot = false;
		bool left = false;
		bool right = false;

		int i = (int) floor(c.x / map_width);
		int j = (int) floor(c.y / map_height);

		fast_remove(mp[i][j], c);

		if (j > (int)floor((c.y - c.radius) / map_height)) {
			top = true;
			if (j > 0) fast_remove(mp[i][j - 1], c);
		}

		if (j < (int)floor((c.y + c.radius) / map_height)) {
			bot = true;
			fast_remove(mp[i][j + 1], c);
		}

		if (i >(int)floor((c.x - c.radius) / map_width)) {
			left = true;
			if (i > 0) fast_remove(mp[i - 1][j], c);
		}

		if (i < (int)floor((c.x + c.radius) / map_width)) {
			right = true;
			fast_remove(mp[i + 1][j], c);
		}

		if (top && left) {
			double x = i * map_width;
			double y = j * map_height;
			if (i > 0) if (j > 0) if (sqr(x - c.x) + sqr(y - c.y) < c.radius) fast_remove(mp[i - 1][j - 1], c);
		}

		if (top && right) {
			double x = (i + 1) * map_width;
			double y = j * map_height;
			if (j > 0) if (sqr(x - c.x) + sqr(y - c.y) < c.radius) fast_remove(mp[i + 1][j - 1], c);
		}

		if (bot && left) {
			double x = i * map_width;
			double y = (j + 1) * map_height;
			if (i > 0) if (sqr(x - c.x) + sqr(y - c.y) < c.radius) fast_remove(mp[i - 1][j + 1], c);
		}

		if (bot && right) {
			double x = (i + 1) * map_width;
			double y = (j + 1) * map_height;
			if (sqr(x - c.x) + sqr(y - c.y) < c.radius) fast_remove(mp[i + 1][j + 1], c);
		}
	}
}

#define TOP 1
#define MID 2
#define BOT 3
#define OTHER 4
int getUnitLane(double x, double y) {
	if ((x < 800) || (y < 800)) return TOP;
	if ((x > 3200) || (y > 3200)) return BOT;
	if (fabs(x - y) < 400) return MID;

	return OTHER;
}

point getTopFront(const World& world) {
	double maxx = 0, miny = 4000.0;

	auto& minions = world.getMinions();
	for (auto& m : minions) {
		if (m.getFaction() == myFaction) {
			double mx = m.getX();
			double my = m.getY();
			if (getUnitLane(mx, my) == TOP) {
				mx += m.getVisionRange();
				my -= m.getVisionRange();
				if (mx > maxx) maxx = mx;
				if (my < miny) miny = my;
			}
		}
	}

	auto& buildings = world.getBuildings();
	for (auto& m : buildings) {
		if (m.getFaction() == myFaction) {
			double mx = m.getX();
			double my = m.getY();
			if (getUnitLane(mx, my) == TOP) {
				mx += m.getVisionRange();
				my -= m.getVisionRange();
				if (mx > maxx) maxx = mx;
				if (my < miny) miny = my;
			}
		}
	}

	auto& wizards = world.getWizards();
	for (auto& m : wizards) {
		if (m.getFaction() == myFaction) {
			double mx = m.getX();
			double my = m.getY();
			if (getUnitLane(mx, my) == TOP) {
				mx += m.getVisionRange();
				my -= m.getVisionRange();
				if (mx > maxx) maxx = mx;
				if (my < miny) miny = my;
			}
		}
	}

	if (miny < 200) return point(maxx, 200);
	return point(200, miny);
}

inline double sigmoid(const double& x, const double& dx, const double& dy) {
	return dy * (1.0 - (x - dx) / (1.0 + fabs(x - dx)));
}

const double dyingXpCoef = 10.0;
const double damageDealCoef = 2.0;
const double staffDamageDealCoef = 3.0;
const double damageTakeFromMinionCoef = 2.5;
const double damageTakeFromWizardCoef = 1.5;
const double scareFactor = 0.03;
const double obstacleCoef = 10.0;
const double projectileCoef = 10.0;

vector<double> attackRangeByType = { 0, 305, 50, 565, 765 };

double getPotential(double x, double y, const World& world, const set<dcircle>& enemiesNear, double myHPpct) {
	double positive = 0;
	double negative = 0;

	if ((x < 0.0) || (x > 4000.0)) return -1e10;
	if ((y < 0.0) || (y > 4000.0)) return -1e10;

	int xi = (int)floor(x / cellSize);
	int yi = (int)floor(y / cellSize);

	double obstacles = 0;
	obstacles += sigmoid(x, 0, obstacleCoef);
	obstacles += sigmoid(y, 0, obstacleCoef);

	for (auto& obj : staticMap[xi][yi]) {
		obstacles +=  sigmoid(obj.distanceTo(x, y), obj.radius, obstacleCoef);
	}

	for (auto& obj : dynamicMap[xi][yi]) {
		obstacles += sigmoid(obj.distanceTo(x, y), obj.radius, obstacleCoef);
	}

	set<double, std::greater<double> > enemiesAttackable, enemiesAggroed;

	for (auto&e : enemiesNear) {
		double dist = e.distanceTo(x,y);
		positive += sigmoid(dist, 600, dyingXpCoef / (0.1 + floor(e.hp / 12.0))); //!!!!!!!!!!!!!! 600 - радиус полчения экспы. Возможно, нужно считать не расстояние между центрами, а расстояние от своего центра до любой точки умирающего
		
		//мы хотим атаковать (стрелять)
		enemiesAttackable.insert(sigmoid(dist, (500 + 10 - 38) + e.radius, damageDealCoef)); // 500 - дальность моей атаки, 10 - радиус снаряда, 38 - дополнительный радиус дин. объекта
		
		//мы хотим бить руками
		enemiesAttackable.insert(sigmoid(dist, 70 - 38 + e.radius, staffDamageDealCoef)); // 70 - дальность моей атаки
			
		if (e.type != 3) {
			enemiesAggroed.insert(sigmoid(dist, e.attackRange, damageTakeFromMinionCoef / (scareFactor + myHPpct)));
		}
		else {
			enemiesAggroed.insert(sigmoid(dist, e.attackRange + 20.0, damageTakeFromWizardCoef / (scareFactor + myHPpct)));
			enemiesAggroed.insert(sigmoid(dist, 105.0, damageTakeFromWizardCoef / (scareFactor + myHPpct)));
		}
	}
	double mod = 1.0;
	for (auto& d : enemiesAttackable) {
		positive += mod * d;
		mod /= 3.0;
	}
	mod = 1.0;
	for (auto& d : enemiesAggroed) {
		negative += mod * d;
		mod *= 2.5;
	}
	//!!!!!!!!!!!!!!!!! отходить от летящих снарядов
	auto& proj = world.getProjectiles();

	for (auto&p : proj) {
		if (p.getFaction() != myFaction) {
			double A = p.getSpeedY();
			double B = -p.getSpeedX();

			if (-(x - p.getX()) * B + (y - p.getY()) * A > 0){
				double C = -(A * p.getX() + B * p.getY());

				negative += sigmoid(fabs(A * x + B * y + C) / sqrt(A * A + B * B), 45, projectileCoef / (scareFactor + myHPpct));
			}
		}
	}

	return positive - negative - obstacles;
}

long long cum_move_time = 0;
int tick = 0;
set<long long> treesSet;
point destination(0,0), top_rune(1200, 1200), bot_rune(2800, 2800);
#define PUSH_TOP 1
int cur_strategy = PUSH_TOP;

#define GET_RUNE 1
#define PUSH 2
int cur_aim = PUSH;

const double nearTreshold = 650; 

void MyStrategy::move(const Wizard& self, const World& world, const Game& game, Move& move) {
	bool inBattle = false;
	auto move_start_time = time();
	set<circle> staticObjects, newStaticObjects, killedStaticObjects;
	set<dcircle> dynamicObjects;
	auto buildings = world.getBuildings();
	auto trees = world.getTrees();
	auto minions = world.getMinions();
	auto players = world.getWizards();
	set<circle> toRemove;	//!!!!!!!!!!!!!!!!!! здесь скорее всего лучше вектора использовать, ане сеты
	set <dcircle> enemies, enemiesNear;

	vector<CircularUnit> enemiesNearCircular;
	enemiesNearCircular.reserve(20);
	
	bool newTrees = false;
	double selfx = self.getX();
	double selfy = self.getY();

	tick++;

	myFaction = self.getFaction();
	if (self.getFaction() == FACTION_ACADEMY) {
		enFaction = FACTION_RENEGADES;
	}
	else {
		enFaction = FACTION_ACADEMY;
	}

	// ========================================================= Prepare map ===================================================================================================================================================================
	{
		for (auto &row : dynamicMap) {
			for (auto& cell : row) {
				cell.resize(0);
			}
		}


		for (auto& b : buildings) {
			staticObjects.insert(b);
			if (b.getFaction() == enFaction) {
				enemiesNearCircular.push_back(b);
				enemiesNear.insert(b);
			}
		}
		int sizeWas = treesSet.size();
		for (auto& t : trees) {
			staticObjects.insert(t);
			treesSet.insert(t.getId());
		}
		if (treesSet.size() > sizeWas) newTrees = true;

		for (auto & m : minions) {
			dcircle d(m);
			if (sqr(selfx - d.x) + sqr(selfy - d.y) < sqr(d.radius + 1e-2)) {
				d.radius -= 3;
			}
			dynamicObjects.insert(d);

			if (m.getFaction() == enFaction) {
				enemies.insert(d);
				if (sqr(selfx - d.x) + sqr(selfy - d.y) < sqr(d.radius - 30 + nearTreshold)) {
					enemiesNear.insert(d);
					enemiesNearCircular.push_back(m);
				}
			}
		}

		for (auto & m : players) {
			if (!m.isMe()) {
				dcircle d(m);
				if (sqr(selfx - d.x) + sqr(selfy - d.y) < sqr(d.radius + 1e-2)) {
					d.radius -= 3;
				}
				dynamicObjects.insert(d);

				if (m.getFaction() == enFaction) {
					enemies.insert(d);
					if (sqr(selfx - d.x) + sqr(selfy - d.y) < sqr(d.radius - 30 + nearTreshold)) {
						enemiesNear.insert(d);
						enemiesNearCircular.push_back(m);
					}
				}
			}
		}


		set_difference(staticObjects.begin(), staticObjects.end(), prevStaticObjects.begin(), prevStaticObjects.end(), inserter(newStaticObjects, newStaticObjects.end()));
		set_difference(prevStaticObjects.begin(), prevStaticObjects.end(), staticObjects.begin(), staticObjects.end(), inserter(killedStaticObjects, killedStaticObjects.end()));

		prevStaticObjects = staticObjects;

		addObjectsToMap(staticMap, newStaticObjects, cellSize, cellSize);
		addObjectsToMap(dynamicMap, dynamicObjects, cellSize, cellSize);

		// ----------- remove killed objects

		for (auto& c : killedStaticObjects) {
			//!!!!!!!!!!!!!!!!!!!!!!!!!! check if any dynamic obj. or *building* should see that object
			auto cell = cellByXY(c.x, c.y, dynamicMap);
			bool insight = false;
			for (auto& d : (*cell)) {
				if (sqr(d.sightRange) < sqr(d.x - c.x) + sqr(d.y - c.y)) {
					toRemove.insert(c);
					insight = true;
					break;
				}
			}
			if (insight) break;
			for (auto& b : buildings) {
				if (sqr(b.getVisionRange()) < sqr(b.getX() - c.x) + sqr(b.getY() - c.y)) {
					toRemove.insert(c);
					break;
				}
			}
		}

		removeObjectsFromMap(staticMap, toRemove, cellSize, cellSize);
	}

	// ========================================================= Combat? ===================================================================================================================================================================
#pragma region "Combat"
	inBattle = false;
	if (enemiesNear.size() > 0) {
		inBattle = true;
	}

#pragma endregion


	// ========================================================= Global Movement ===================================================================================================================================================================
	if(!inBattle) {
		//!!!!!!!!!!!!!!!! CHECK RUNE

		//assume we shouldn't go for the rune
		if (cur_strategy == PUSH_TOP) {
			destination = getTopFront(world);
			cout << "DEST: " << floor(destination.x) << " \t" << destination.y << endl;
		}
	}

	// ========================================================= Micro ===================================================================================================================================================================
#pragma region "Micro"
	if (inBattle) {
		double ts = time();

		double myHPpct = ((double)self.getLife()) / ((double)self.getMaxLife());
		auto best = getPotential(selfx, selfy, world, enemiesNear, myHPpct);
		cout << "cur potential: " << best << endl;
		point best_pt(selfx, selfy);

		double pot_step = 30.0;
		double search_rad = 10;
		double pot;

#ifdef _zuko3d_pc
		vector<double> pre_bmp;
		int wd = 2 * ((int)search_rad) + 1;
		TGAImage bmp(wd,wd);
		pre_bmp.resize(sqr(wd));
		double xstart = selfx - pot_step * search_rad;
		double ystart = selfy - pot_step * search_rad;
		
#endif
		for (double y(selfy - pot_step * search_rad), y_end(selfy + pot_step * search_rad); y <= y_end; y += pot_step) {
			for (double x(selfx - pot_step * search_rad), x_end(selfx + pot_step * search_rad); x <= x_end; x += pot_step) {
 				pot = getPotential(x, y, world, enemiesNear, myHPpct);
#ifdef _zuko3d_pc
				int xi = (int) ((x - xstart) / pot_step);
				int yi = (int) ((y - ystart) / pot_step);
				//cout << "\t" << floor(pot * 10.0);
				pre_bmp[xi + yi * wd] = pot;
				
#endif
				if (pot > best) {
					best = pot;
					best_pt = point(x, y);
				}
			}
			//cout << endl;
		}
		cout << "best: " << best << " \t coord: " << floor(best_pt.x) << " \t" << floor(best_pt.y) << endl;
		
		curWay = findPathToZone(point(selfx, selfy), best_pt, 10.0, 5.0);
		if (!curWay.empty()) {
			setMoveToPoint(self, move, curWay.front());
		}
		else {
			setMoveToPoint(self, move, best_pt);
		}

		cout << "pot. calc time: " << time() - ts << endl;
#ifdef _zuko3d_pc
		//saturate
		double mx = 0;
		for (auto&p : pre_bmp) {
			if (mx < fabs(p)) mx = fabs(p);
		}
		mx = 254 / mx;
		Colour clr;
		clr.a = 255;
		clr.r = 0;
		clr.g = 0;
		clr.b = 0;
		
		for(unsigned int i = 0; i < pre_bmp.size(); i++){
			if (pre_bmp[i] < 0) {
				//bmp[i * 3] = (BYTE)floor(-pre_bmp[i] * mx);
				clr.r = (unsigned char) floor(-pre_bmp[i] * mx);
				clr.g = 0;
				bmp.setPixel(clr, i % wd, i / wd);
			}
			else {
				clr.g = (unsigned char)floor(pre_bmp[i] * mx);
				clr.r = 0;
				bmp.setPixel(clr, i % wd, i / wd);
			}
		}
		int xi = (int)((best_pt.x - xstart) / pot_step);
		int yi = (int)((best_pt.y - ystart) / pot_step);

		cout << "xi/yi: " << xi << " \t" << yi << endl;

		clr.r = 255;
		clr.g = 255;
		clr.b = 255;
		bmp.setPixel(clr, xi, yi);

		bmp.WriteImage("potential.tga");
		//system("pause");
#endif
	}
#pragma endregion

	// ========================================================= Attack ===================================================================================================================================================================
	{
		if (inBattle) {
			if (!self.getRemainingActionCooldownTicks()) {
				if (self.getRemainingCooldownTicksByAction()[ACTION_MAGIC_MISSILE] < 1) {
					point aim(0,0);
					double best = 0.0;
					for (auto&w : enemiesNear) {
						point pt(w);
						if (fabs(self.getAngleTo(pt)) < PI / 12.0) {
							if (self.getDistanceTo(pt) < 500) {
								double pts = 1000.0 / w.hp;
								if (w.type == 3) pts += 20.0;
								if (w.type == 4) pts += 50.0;
								if (w.type == 5) pts += 60.0;
								if (pts > best) {
									aim = pt;
									best = pts;
								}
							}
						}
					}

					if (best > 0.0) {
						move.setAction(ACTION_MAGIC_MISSILE);
						move.setCastAngle(self.getAngleTo(aim));
						move.setMinCastDistance(self.getDistanceTo(aim) - 8 - 30);
						move.setMaxCastDistance(550);
					}
				}
				else {
					if (self.getRemainingCooldownTicksByAction()[ACTION_STAFF] < 1) move.setAction(ACTION_STAFF);
				}
				double cd = (double) self.getRemainingCooldownTicksByAction()[ACTION_MAGIC_MISSILE];

				point aim(0,0);
				double best = 0;
				for (auto&w : enemiesNear) {
					point pt(w);
					if (self.getDistanceTo(pt) < 530) {
						double pts = (5.0 + max(0.0, self.getAngleTo(pt) / PI * 30.0 - cd)) / w.hp;
						if (pts > best) {
							aim = pt;
							best = pts;
						}
					}
				}

				if (best > 0) {
					move.setTurn(self.getAngleTo(aim));
				}
			}
		}
	}

	// ========================================================= Out-of-combat Movement ===================================================================================================================================================================
	{
		if (!inBattle) {
			if (!curWay.empty()) {
				bool needToRecountPath = false;
				auto dobjects = cellByXY(self.getX(), self.getY(), dynamicMap);
				if (dobjects->size() > 0) {
					bool tooClose = false;
					for (auto& d : (*dobjects)) {
						if (sqrt(sqr(selfx - d.x) + sqr(selfy - d.y)) < d.radius + 10.0) {
							tooClose = true;
							break;
						}
					}
					if (tooClose) {
						cout << "too close dyn objects here: " << dobjects->size() << endl;
						needToRecountPath = true;
					}
				}
				if (self.getDistanceTo(curWay.front()) < 1e-2) {
					cout << "recalc route, too close to WP: " << self.getX() << " \t" << self.getY() << endl;
					needToRecountPath = true;
				}
				if ((tick > 1) && (fabs(self.getSpeedX()) + fabs(self.getSpeedY()) < 1e-2)) {
					needToRecountPath = true;
					cout << "We are struck!" << endl;
				}
				if (newTrees) {
					needToRecountPath = true;
					cout << "New tree spawned!" << endl;
				}
				if (needToRecountPath) {
					auto ts = time();
					static long long cumulative = 0;
					if (cur_aim == PUSH) {
						curWay = findPathToZone(point(self.getX(), self.getY()), destination, 400.0);
					}
					else {
						curWay = findPath(point(self.getX(), self.getY()), destination);
					}
					cumulative += time() - ts;
					cout << "time spent for path: " << (time() - ts) / 1000000 << " \t " << cumulative / 1000000 << endl;
					cout << "next WP: " << curWay.front().x << " \t" << curWay.front().y << endl;

#ifndef _DEBUG
					//system("pause");
#endif
				}
				if (!curWay.empty()) {
					setMoveToPoint(self, move, curWay.front());

					if (!inBattle) {
						move.setTurn(self.getAngleTo(curWay.front()));
					}
				}
			}
			else {
				curWay = findPath(point(selfx, selfy), destination);
			}
		}
	}

	cum_move_time += time() - move_start_time;
	//cout << tick << ":\t Time per tick(us): " << cum_move_time / tick / 1000 << endl;
}

MyStrategy::MyStrategy() { 
	staticMap.resize(5 + (int)ceil(4000.0 / cellSize));
	for (auto &row : staticMap) {
		row.resize(5 + (int)ceil(4000.0 / cellSize));
		for (auto& cell : row) {
			cell.reserve(100);
		}
	}

	dynamicMap.resize(5 + (int)ceil(4000.0 / cellSize));
	for (auto &row : dynamicMap) {
		row.resize(5 + (int)ceil(4000.0 / cellSize));
		for (auto& cell : row) {
			cell.reserve(100);
		}
	}
}
