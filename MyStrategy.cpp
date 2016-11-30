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
	return m_width * (m_height - y - 1) + x;
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
double maxSpeedX, maxSpeedY;
double myAttackRange;
vector<model::Bonus> runes;
vector<model::Projectile> projNear;

#ifdef _zuko3d_output_path_stats
#include <iostream>
ofstream fout_path("path_stats.txt", ios_base::app);
#endif

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
	struct circle;
	struct dcircle;
}

class point : public Unit {
public:
	double x, y;

	point();
	point(double a, double b);
	point(const dcircle& c);
};

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
struct std::circle {
	double radius = 0.0;
	double x = 0.0, y = 0.0;
	long long id;
	bool isTree = false;

	circle() { }

	circle(const CircularUnit& u, double additionalRadius = 35, bool tree = false) : radius(u.getRadius() + additionalRadius), x(u.getX()), y(u.getY()), id(u.getId()), isTree(tree) { }

	bool operator<(const circle& c) const {
		return id < c.id;
	}

	bool operator==(const circle& c) const {
		return id == c.id;
	}
	
	double distanceTo(double px, double py) const {
		return sqrt(sqr(x - px) + sqr(y - py));
	}

	double distanceTo(const point& p) const {
		return distanceTo(p.x, p.y);
	}
};

list<point> pathToBonus;
double pathLenghtToBonus = 1e10;

struct std::dcircle : circle {
	double sightRange = 0.0;
	int type;
	double attackRange; // максимальное расстояние между центрами, при котором можно нанести урон
	double angle;
	double hpPct;
	double hp;
	double cd;
	int faction;
	double refugeDistance = 0.0;
	double predicted_x, predicted_y;

	dcircle() { }

	dcircle(const LivingUnit& u, double additionalRadius = 38) : circle(u, additionalRadius), predicted_x(u.getX() + 3.0 * u.getSpeedX()), predicted_y(u.getY() + 3.0 * u.getSpeedY()), hpPct((double) u.getLife() / (double) u.getMaxLife()), hp((double) u.getLife()), angle(u.getAngle()), faction(u.getFaction()) { }

	dcircle(const Minion& u, double additionalRadius = 38) : dcircle((LivingUnit) u, additionalRadius) {
		sightRange = u.getVisionRange();
		cd = u.getRemainingActionCooldownTicks();
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
		attackRange = 510 + 35; // !!!!!!!!!!!!!!!!!!!!!! дальность атаки увеличивается с уровнем!
		cd = u.getRemainingCooldownTicksByAction()[ACTION_MAGIC_MISSILE];
	}
	dcircle(const Building& u, double additionalRadius = 38) : dcircle((LivingUnit)u, additionalRadius) {
		sightRange = u.getVisionRange();
		cd = u.getRemainingActionCooldownTicks();
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

	double getAngleTo(const dcircle& c) const {
		return getAngleTo(c.x, c.y);
	}

	double getAngleTo(double px, double py) const {
		double dx = px - x;
		double dy = py - y;
		
		double ret = atan2(dx, dy) - angle;
		if (ret > PI) ret -= 2.0 * PI;

		return ret;
	}
};

point::point() : Unit(0, 0, 0, 0, 0, 0, FACTION_ACADEMY) { }

point::point(double a, double b) : Unit(0, a, b, 0, 0, 0, FACTION_ACADEMY) {
	x = a;
	y = b;
}

point::point(const dcircle& c) : Unit(0, c.x, c.y, 0, 0, 0, FACTION_ACADEMY) {
	x = c.x;
	y = c.y;
}


struct circleSetComarator {
	bool operator() (const circle & c1, const circle& c2) {
		return c1.id < c2.id;
	}
};

double cellSize = 150.0; // не меньше, чем радиус самого большого объекта (85 при данных правилах, не считая базы)
set<circle> prevStaticObjects;

vector< vector< vector< circle > > > staticMap;
vector< vector< vector< dcircle > > > dynamicMap;

vector<vector<point> > refugePoints;
point refugePoint;

double maxMSX(const Wizard& w) {
	return 3.0;// !!!!!!!!!!!!!!!!!!! НУЖНО УЧИТЫВАТЬ БАФФЫ на скорость
}

double maxMSY(const Wizard& w) {
	return 4.0;// !!!!!!!!!!!!!!!!!!! НУЖНО УЧИТЫВАТЬ БАФФЫ на скорость
}

bool checkCollision(const point& a, const point& b, const circle& c, double d_rad = 0.0) {
	if (sqr(a.x - c.x) + sqr(a.y - c.y) < sqr(c.radius + d_rad)) return true;
	if (sqr(b.x - c.x) + sqr(b.y - c.y) < sqr(c.radius + d_rad)) return true;

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

	if (fabs(A * c.x + B * c.y + C) < c.radius + d_rad + 1e-2) return true;
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

bool checkCollisions(double ax, double ay, double d_rad = 0.0) {
	if ((ax < 0.0) || (ax > 4000.0)) return false;
	if ((ay < 0.0) || (ay > 4000.0)) return false;

	int xi = (int)floor(ax / cellSize);
	int yi = (int)floor(ay / cellSize);

	for (auto& obj : staticMap[xi][yi]) {
		if (obj.distanceTo(ax, ay) < obj.radius - d_rad) return true;
	}

	for (auto& obj : dynamicMap[xi][yi]) {
		if (obj.distanceTo(ax, ay) < obj.radius - d_rad) return true;
	}

	return false;
}

bool checkCollisions(const point& pt, double d_rad = 0.0) {
	return checkCollisions(pt.x, pt.y, d_rad);
}

template<class circle_class> vector<circle_class>* cellByXY(double x, double y, vector< vector< vector<circle_class> > > &objects) {
	return &(objects[(int)floor(x / cellSize)][(int)floor(y / cellSize)]);
}

bool checkCollisions(double ax, double ay, double bx, double by, bool treesOnly = false, double d_rad = 0.0) {
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
		if(i > 0) cells.insert(pt(i - 1, j));
		i += di;
		sx += dx;
		if (i < -1) {
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
		if(j > 0) cells.insert(pt(i, j - 1));
		j += dj;
		sy += dy;
		if (i < 0) {
			cout << "ERROR i < 0!" << endl;
		}
		if (j < -1) {
			cout << "ERROR j < 0!" << endl;
		}
	}

	point a(ax, ay), b(bx, by);
	for (const auto&c : cells) {
		for (const auto& obj : staticMap[c.first][c.second]) {
			if ((!treesOnly) || obj.isTree) {
				if (checkCollision(a, b, obj, d_rad)) return true;
			}
		}
		if (!treesOnly) {
			for (const auto& obj : dynamicMap[c.first][c.second]) {
				if (checkCollision(a, b, obj, d_rad)) return true;
			}
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

/*
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
*/

list<point> findPathToZonePrecise(const point& start, const point& dest, double rad, double step = 20.0, int max_Ticks = 50000) {
	typedef pair<int, int> pt;
	typedef pair<double, pt > qtype;

	auto h = [rad](const pt&a, const pt&b) -> double {return max(0.0, sqrt((double)(sqr(a.first - b.first) + sqr(a.second - b.second))) - rad); };

	list<point> path;
	if (checkCollisions(dest, rad)) return path;
	if (!checkCollisions(start, dest)) {
		path.push_back(dest);
		return path;
	}

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
		if (cur_path_ticks > max_Ticks) break;
		cur = open.top();
		while (closed.find(cur.second) != closed.end()) {
			open.pop();
			if (open.empty()) {
				break;
			}
			cur = open.top();
		}

		if (sqr(cur.second.first * step + start.x - dest.x) + sqr(cur.second.second * step + start.y - dest.y) <= sqr(rad)) {
			good = true;
			break;
		}

		if (open.empty()) break;

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

	{
		cout << "Dist: " << start.getDistanceTo(dest) << " (rad = " << rad << ") , ticks: \t" << cur_path_ticks << endl;
#ifdef _zuko3d_output_path_stats
		fout_path << start.getDistanceTo(dest) << " \t" << rad << " \t" << cur_path_ticks << endl;
		fout_path.flush();
#endif // _z
	}

	if (!good) {
		cout << "!!!!!!!!! FAILED to find a path" << endl;

#ifdef _zuko3d_pc
		int minx = 4000, maxx = 0, miny = 4000, maxy = 0;
		for (auto& p : openSet) {
			if (minx > p.first) minx = p.first;
			if (maxx < p.first) maxx = p.first;
			if (miny > p.second) miny = p.second;
			if (maxy < p.second) maxy = p.second;
		}

		for (auto& p : closed) {
			if (minx > p.first) minx = p.first;
			if (maxx < p.first) maxx = p.first;
			if (miny > p.second) miny = p.second;
			if (maxy < p.second) maxy = p.second;
		}


		TGAImage path_pic;

		//path_pic.WriteImage("path.tga");

		//system("pause");
#endif
		return path;
	}

	
	if (cur_path_ticks > path_ticks) {
		path_ticks = cur_path_ticks;
	}

	cout << "path_ticks = \t" << path_ticks << "\t" << cur_path_ticks << endl;
	

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

list<point> findPathToZone(const point& start, const point& dest, double rad, double step = 20.0, int max_Ticks = 50000) {
	list<point> ret;
	ret.clear();

	while (ret.empty()) {
		ret = findPathToZonePrecise(start, dest, rad, step, max_Ticks);
		rad = (rad + 5.0) * 2.0;
	}

	return ret;
}

list<point> findPath(const point& start, const point& dest, double step = 20.0) {
	if (checkCollisions(dest)) return list<point>();
	return findPathToZone(start, dest, 0.05, step);
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
				mx += m.getVisionRange() / 4.0;
				my -= m.getVisionRange() / 4.0;
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
				mx += m.getVisionRange() / 4.0;
				my -= m.getVisionRange() / 4.0;
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
				mx += m.getVisionRange() / 4.0;
				my -= m.getVisionRange() / 4.0;
				if (mx > maxx) maxx = mx;
				if (my < miny) miny = my;
			}
		}
	}

	if (miny < 200) return point(maxx, 200);
	return point(200, miny);
}

inline double sigmoid(const double& x, const double& dx, const double& dy) { // 0 ... 2
	return dy * (1.0 - (x - dx) / (1.0 + fabs(x - dx)));
}


inline double sigmoidStrict(const double& x, const double& dx, const double& dy) { // 0 ... 2
	return dy * (1 - tanh(x - dx));
}


int top_rune_cond = 0;
int bot_rune_cond = 0;
point destination(0, 0), top_rune(1200, 1200), bot_rune(2800, 2800);

const double dyingXpCoef = 10.0;
const double damageDealCoef = 2.0;
const double staffDamageDealCoef = 7.5;
const double damageTakeFromMinionCoef = 2.5;
const double damageTakeFromWizardCoef = 1.5;
const double scareFactor = 0.03;
const double obstacleCoef = 500.0;
const double projectileCoef = 20.0;
const double catchFearCoef = 0.01;


vector<double> attackRangeByType = { 0, 305, 50, 565, 765 };

double getPotential(double x, double y, const World& world, const vector<dcircle>& enemiesNear, double myHPpct, const Wizard& self) {
	double positive = 0;
	double negative = 0;
	double cought = 0.0;

	auto cd = self.getRemainingCooldownTicksByAction()[ACTION_MAGIC_MISSILE];

	if ((x < 0.0) || (x > 4000.0)) return -10.0 * obstacleCoef;
	if ((y < 0.0) || (y > 4000.0)) return -10.0 * obstacleCoef;

	int xi = (int)floor(x / cellSize);
	int yi = (int)floor(y / cellSize);

	double obstacles = 0;
	obstacles += sigmoidStrict(x, 35, obstacleCoef);
	obstacles += sigmoid(x, 0, obstacleCoef);
	obstacles += sigmoidStrict(y, 0, obstacleCoef);

	for (auto& obj : staticMap[xi][yi]) {
		obstacles += sigmoidStrict(obj.distanceTo(x, y), obj.radius, obstacleCoef);
	}

	for (auto& obj : dynamicMap[xi][yi]) {
		obstacles += sigmoidStrict(obj.distanceTo(x, y), obj.radius, obstacleCoef);
	}

	set<double, std::greater<double> > enemiesAttackable, enemiesAggroed;

	double timeToCome = ceil(sqrt(sqr(x - self.getX()) + sqr(y - self.getY())) / maxSpeedX);
	double myRefugeDistance = refugePoint.getDistanceTo(x, y) + 50.0 / sqr(myHPpct);

	for (auto&e : enemiesNear) {
		double dist = e.distanceTo(x,y);
		double predicted_dist = sqrt(sqr(e.predicted_x - x) + sqr(e.predicted_y - y));

		positive += sigmoid(dist, 600, dyingXpCoef / (0.3 + floor(e.hp / 12.0))); //!!!!!!!!!!!!!! 600 - радиус полчения экспы. Возможно, нужно считать не расстояние между центрами, а расстояние от своего центра до любой точки умирающего
		
		if (e.faction == enFaction) {
			cought += max(-30.0, myRefugeDistance - e.refugeDistance) / (scareFactor + myHPpct) * (e.type == 3 ? 4.0: 1.0); // !!!!!!!!! только для топ-линии! для остальных нужно свои формулы делать
		}
		double angleCoef = sigmoid(fabs(e.getAngleTo(x, y)), PI / 12.0 * (1.0 + timeToCome * 0.4), 0.5);
		//double myAngleCoef = sigmoid(fabs(e.getAngleTo(x, y)), PI / 12.0 * (1.0 + timeToCome * 0.4), 0.5);
		//хотим бить руками
		enemiesAttackable.insert(sigmoid(dist, 70 - 38 - 5 + e.radius, staffDamageDealCoef * myHPpct)); // !!!!!!!!!!!!!!!!!! не учитываем угол для удара!

		if (!checkCollisions(x, y, e.x, e.y, true, -25.0)) {
			if (cd <= timeToCome) {
				//мы хотим атаковать (стрелять)
				enemiesAttackable.insert(sigmoid(predicted_dist, myAttackRange - 51 + e.radius, damageDealCoef / (1.0 - sigmoid(e.hp, 20.0, 0.5)))); // 500 - дальность моей атаки, 10 - радиус снаряда, 38 - дополнительный радиус дин. объекта
			}

			if (angleCoef > 0.1) {
				double tmp1 = e.attackRange - max(0.0, e.cd - 4.0 / myHPpct) * maxSpeedX + 5.0 / myHPpct + timeToCome; // добавочное расстояние за счёт того, что к нам могут подойдти
				if (e.type != 3) {
					//!!!!!!!!!!!!!!!!!!! можно предсказывать, кого моб будет бить
					enemiesAggroed.insert(sigmoid(dist, tmp1, angleCoef * damageTakeFromMinionCoef / (scareFactor + myHPpct)));
				}
				else {
					enemiesAggroed.insert(sigmoid(dist, tmp1, angleCoef * damageTakeFromWizardCoef / (scareFactor + myHPpct)));
					enemiesAggroed.insert(sigmoid(dist, 105.0, angleCoef * damageTakeFromWizardCoef / (scareFactor + myHPpct)));
				}
			}
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

	for (auto&p : projNear) {
		double A = p.getSpeedY();
		double B = -p.getSpeedX();

		if (-(x - p.getX()) * B + (y - p.getY()) * A > 0){
			double C = -(A * p.getX() + B * p.getY());

			negative += sigmoid(fabs(A * x + B * y + C) / sqrt(A * A + B * B), 45, projectileCoef / (scareFactor + myHPpct));
		}
	}

	if(top_rune_cond > 0) { //!!!!!!!!!!!!!!!! для бот-руны тоже надо!
		//double dist = sqrt(sqr(x - top_rune.getX()) + sqr(y - top_rune.getY()));
		//double dist = fabs(x - top_rune.getX()) + fabs(y - top_rune.getY());
		double dist = sqrt(sqr(pathToBonus.front().x - x) + sqr(pathToBonus.front().y - y));
		positive += max(0.0, 3000.0 - pathLenghtToBonus - dist ) * 3.0;
	}

	return positive - negative - obstacles - cought * catchFearCoef;
}

double pathLength(const list<point> & path, const Wizard& self) {
	list<point>::const_iterator it, it_next, it_end;
	double ret = self.getDistanceTo(path.front());

	it = path.begin();
	it_end = path.end();
	it_next = it;
	it_next++;

	while (it_next != it_end) {
		ret += it->getDistanceTo(*it_next);

		it++;
		it_next++;
	}
	return ret;
}

long long cum_move_time = 0;
int tick = 0;
set<long long> treesSet;
#define PUSH_TOP 1
int cur_strategy = PUSH_TOP;

#define GET_RUNE 1
#define PUSH 2
int cur_aim = PUSH;

const double nearTreshold = 650; 
int prev_tick = 0;

vector<model::SkillType> skillBuild;

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
	vector <dcircle> enemies, enemiesNear;
	enemies.reserve(20);
	enemiesNear.reserve(20);

	vector<CircularUnit> enemiesNearCircular;
	enemiesNearCircular.reserve(20);
	auto gtick = world.getTickIndex();
	// ======================================================== Update Variables===================================================================================================================================================================
	
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

	maxSpeedX = maxMSX(self); // !!!!!!!!!!!!!!!!!!! НУЖНО УЧИТЫВАТЬ БАФФЫ на скорость
	maxSpeedY = maxMSY(self); // !!!!!!!!!!!!!!!!!!! НУЖНО УЧИТЫВАТЬ БАФФЫ на скорость

	bool haste = false;
	auto myStatuses = self.getStatuses();
	for (auto&s : myStatuses) {
		if (s.getType() == STATUS_HASTENED) {
			haste = true;
			break;
		}
	}
	if (haste) {
		maxSpeedX *= 1.0 + game.getHastenedMovementBonusFactor();
		maxSpeedY *= 1.0 + game.getHastenedMovementBonusFactor();
	}

	myAttackRange = 500 + 10 - 30 + 25.0 * max(0.0, (double) self.getLevel() - 4.0); 

	auto myLane = getUnitLane(selfx, selfy);
	if (myLane == TOP) { // !!!!!!!!!! другие лэйны тоже надо учитывать
		if (selfy > 600) {
			for (auto&p : refugePoints[TOP]) {
				if (selfy > p.y) break;
				refugePoint = p;
			}
		}
		else {
			refugePoint = refugePoints[TOP].front();
			for (auto&p : refugePoints[TOP]) {
				if (selfx < p.x) break;
				refugePoint = p;
			}
		}
	}

	auto proj = world.getProjectiles();
	projNear.resize(0);
	projNear.reserve(10);
	for (auto&p : proj) {
		if (p.getFaction() != myFaction) {
			if (self.getDistanceTo(p) < 550) {
				projNear.push_back(p);
			}
		}
	}

	// ========================================================= Level up? ===================================================================================================================================================================
	move.setSkillToLearn(skillBuild[self.getLevel()]);

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
				if (sqr(selfx - b.getX()) + sqr(selfy - b.getY()) < sqr(b.getRadius() + nearTreshold)) {
					enemiesNearCircular.push_back(b);
					enemiesNear.push_back(b);
				}
			}
		}
		auto sizeWas = treesSet.size();
		for (auto& t : trees) {
			circle c(t, 35, true);
			staticObjects.insert(c);
			treesSet.insert(t.getId());
		}
		if (treesSet.size() > sizeWas) newTrees = true;

		for (auto & m : minions) {
			dcircle d(m);
			if (sqr(selfx - d.x) + sqr(selfy - d.y) < sqr(d.radius + 1e-2)) {
				d.radius -= 3;
			}
			dynamicObjects.insert(d);

			if (m.getFaction() != myFaction) {
				enemies.push_back(d);
				if (sqr(selfx - d.x) + sqr(selfy - d.y) < sqr(d.radius - 30 + nearTreshold)) {
					enemiesNear.push_back(d);
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
					enemies.push_back(d);
					if (sqr(selfx - d.x) + sqr(selfy - d.y) < sqr(d.radius - 30 + nearTreshold)) {
						enemiesNear.push_back(d);
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

		for (auto&p : enemiesNear) {
			p.refugeDistance = p.distanceTo(refugePoint);
		}
	}

	// ========================================================= Combat? ===================================================================================================================================================================
#pragma region "Combat"
	inBattle = false;
	for (auto& e : enemiesNear) { 
		if (e.faction == enFaction) {// !!!!!!!!!!!!!!! нужно ещё проверять нейтралов, которые зааггрились
			inBattle = true;
			break;
		}
	}

#pragma endregion
	// ========================================================= Rune ===================================================================================================================================================================
#pragma region "Runes"
	runes = world.getBonuses();
	
	if (prev_tick % 2500 > world.getTickIndex() % 2500) {
		bot_rune_cond = top_rune_cond = 2;
	}

	bool seen_top_rune = false, seen_bot_rune = false;
	for (auto& r : runes) {
		if (r.getX() > 2000) {
			seen_bot_rune = true;
		}
		else {
			seen_top_rune = true;
		}
	}

	for (auto&w : players) {
		if (w.getFaction() == myFaction) {
			if (w.getDistanceTo(top_rune) < w.getVisionRange()) {
				if (!seen_top_rune) top_rune_cond = 0;
			}

			if (w.getDistanceTo(bot_rune) < w.getVisionRange()) {
				if (!seen_bot_rune) bot_rune_cond = 0;
			}
		}
	}

	if (top_rune_cond > 0) {
		if ((world.getTickIndex() % 10) == 0) {
			pathToBonus = findPathToZone((point)self, top_rune, 200.0, 20.0, 5000);
			pathLenghtToBonus = pathLength(pathToBonus, self);
			cout << "pathLenghtToBonus = " << pathLenghtToBonus << endl;
		}
	}
	else {
		pathToBonus.clear();
	}
#pragma endregion	
	// ========================================================= Global Movement ===================================================================================================================================================================
	if(!inBattle) {
		if (top_rune_cond > 0) {
			destination = top_rune;
		}
		else {
			if (cur_strategy == PUSH_TOP) {
				destination = getTopFront(world);
			}
		}
	}

	// ========================================================= Micro ===================================================================================================================================================================
#pragma region "Micro"
	if (inBattle) {
		auto ts = time();

		double myHPpct = ((double)self.getLife()) / ((double)self.getMaxLife());
		auto best = getPotential(selfx, selfy, world, enemiesNear, myHPpct, self);
		cout << "cur potential: " << best << endl;
		point best_pt(selfx, selfy);

		double pot_step = 40.0;
		double search_rad = 7;
		double pot;

#ifdef _zuko3d_pc
		vector<double> pre_bmp;
		int wd = 2 * ((int)search_rad) + 1;
		TGAImage bmp(wd,wd);
		pre_bmp.resize(sqr(wd));
#endif

		double xstart = selfx - pot_step * search_rad;
		double ystart = selfy - pot_step * search_rad;
		for (double y(ystart), y_end(selfy + pot_step * search_rad); y <= y_end; y += pot_step) {
			for (double x(xstart), x_end(selfx + pot_step * search_rad); x <= x_end; x += pot_step) {
#ifdef _zuko3d_pc
				int xi = (int)floor(0.5 + (x - xstart) / pot_step);
				int yi = (int)floor(0.5 + (y - ystart) / pot_step);
#endif
 				pot = getPotential(x, y, world, enemiesNear, myHPpct, self);
#ifdef _zuko3d_pc
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
		
		curWay = findPathToZone(point(selfx, selfy), best_pt, pot_step, 10.0, 3000);
		if (!curWay.empty()) {
			setMoveToPoint(self, move, curWay.front());
			destination = curWay.front();
		}
		else {
			setMoveToPoint(self, move, best_pt);
			destination = best_pt;
		}

		cout << "pot. calc time: " << time() - ts << endl;
#ifdef _zuko3d_pc
		//saturate
		double mx = 0;
		for (auto&p : pre_bmp) {
			if (mx < fabs(p)) mx = fabs(p);
		}
		mx = 254 / log(mx + 1.0);
		Colour clr;
		clr.a = 255;
		clr.r = 0;
		clr.g = 0;
		clr.b = 0;
		
		for(unsigned int i = 0; i < pre_bmp.size(); i++){
			if (pre_bmp[i] < 0) {
				//bmp[i * 3] = (BYTE)floor(-pre_bmp[i] * mx);
				clr.r = (unsigned char) floor(log(-pre_bmp[i] + 1.0) * mx);
				clr.g = 0;
				bmp.setPixel(clr, i % wd, i / wd);
			}
			else {
				clr.g = (unsigned char)floor(log(pre_bmp[i] + 1.0) * mx);
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
					point aim(0, 0);
					double best = 0.0;
					for (auto&w : enemiesNear) {
						point pt(w.predicted_x, w.predicted_y);
						if (fabs(self.getAngleTo(pt)) < PI / 12.0) {
							if (self.getDistanceTo(pt) < myAttackRange + w.radius - 50) { // 38 - доп. радиус, -12 чтобы было тяжелее увернуться
								if (!checkCollisions(selfx, selfy, pt.x, pt.y, true, -25.0)) {
									double pts = 1000.0 / w.hp;
									if (w.type == 3) pts += 10.0;
									if (w.type == 4) pts += 15.0;
									if (w.type == 5) pts += 20.0;
									if (pts > best) {
										aim = pt;
										best = pts;
									}
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
					if (self.getRemainingCooldownTicksByAction()[ACTION_STAFF] < 1) {
						for (auto&e : enemiesNear) {
							if (e.distanceTo((point)self) < 30 + e.radius) {
								point pt(e);
								if (fabs(self.getAngleTo(pt)) <= PI / 12.0) {
									move.setAction(ACTION_STAFF);
									break;
								}
							}
						}
					}
				}
			}
			
			double cd = max((double) self.getRemainingActionCooldownTicks(), (double)self.getRemainingCooldownTicksByAction()[ACTION_MAGIC_MISSILE]);
			double staff_cd = max((double)self.getRemainingActionCooldownTicks(), (double)self.getRemainingCooldownTicksByAction()[ACTION_STAFF]);
			double all_cd = min(staff_cd, cd);

			point aim(0, 0);
			double best = 0;
			for (auto&w : enemiesNear) {
				point pt(w.predicted_x, w.predicted_y);
				double dist = self.getDistanceTo(pt);
				if (dist < 540) {
					if (!checkCollisions(selfx, selfy, pt.x, pt.y, true, -25.0)) {
						double pts = (5.0 + max(0.0, self.getAngleTo(pt) / PI * 30.0 - cd)) / w.hp;
						if (dist < 30 + w.radius) pts *= 1.5 + max(0.0, cd - staff_cd) / 3.0; // we can strike it with STAFF!
						if (w.type == 3) pts += 10.0;
						if (w.type == 4) pts += 15.0;
						if (w.type == 5) pts += 20.0;
						if (pts > best) {
							aim = pt;
							best = pts;
						}
					}
				}
			}

			if (best > 0) {
				double aim_angle = self.getAngleTo(aim);
				if (fabs(aim_angle / PI * 12.0) < all_cd) {
					move.setTurn(self.getAngleTo(destination));
				}
				else {
					move.setTurn(aim_angle);
				}
			}
			
		}
	}

	// ========================================================= Out-of-combat Movement ===================================================================================================================================================================
	{
		if (!inBattle) {
			bool needToRecountPath = false;
			if (!curWay.empty()) {

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
						if (checkCollisions((point) self, curWay.front())) {
							cout << "too close dyn objects here: " << dobjects->size() << endl;
							needToRecountPath = true;
						}
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
			}
			else {
				needToRecountPath = true;
			}
			if (needToRecountPath) {
				auto ts = time();
				static long long cumulative = 0;
				if (cur_aim == PUSH) {
					curWay = findPathToZone(point(self.getX(), self.getY()), destination, 300.0);
				}
				else {
					curWay = findPath(point(self.getX(), self.getY()), destination);
				}
				cumulative += time() - ts;
				//cout << "time spent for path: " << (time() - ts) / 1000000 << " \t " << cumulative / 1000000 << endl;
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
	}

	cum_move_time += time() - move_start_time;
	prev_tick = world.getTickIndex();
	//cout << prev_tick << ":\t Time per tick(us): " << cum_move_time / tick / 1000 << endl;
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

	refugePoints.resize(5);
	vector<point> tmp;
	tmp.reserve(10);
	tmp.push_back(point(100, 3700));
	tmp.push_back(point(200, 3100));
	tmp.push_back(point(200, 1800));
	tmp.push_back(point(200, 900));
	tmp.push_back(point(750, 350));

	refugePoints[TOP] = tmp;

	tmp.resize(0);
	skillBuild.push_back(SKILL_MAGICAL_DAMAGE_BONUS_PASSIVE_1);
	skillBuild.push_back(SKILL_MAGICAL_DAMAGE_BONUS_PASSIVE_1);
	skillBuild.push_back(SKILL_MAGICAL_DAMAGE_BONUS_AURA_1);
	skillBuild.push_back(SKILL_MAGICAL_DAMAGE_BONUS_PASSIVE_2);
	skillBuild.push_back(SKILL_MAGICAL_DAMAGE_BONUS_AURA_2);
	skillBuild.push_back(SKILL_RANGE_BONUS_PASSIVE_1);
	skillBuild.push_back(SKILL_RANGE_BONUS_AURA_1);
	skillBuild.push_back(SKILL_RANGE_BONUS_PASSIVE_2);
	skillBuild.push_back(SKILL_RANGE_BONUS_AURA_2);
	skillBuild.push_back(SKILL_ADVANCED_MAGIC_MISSILE);
	skillBuild.push_back(SKILL_STAFF_DAMAGE_BONUS_PASSIVE_1);
	skillBuild.push_back(SKILL_STAFF_DAMAGE_BONUS_AURA_1);
	skillBuild.push_back(SKILL_STAFF_DAMAGE_BONUS_PASSIVE_2);
	skillBuild.push_back(SKILL_STAFF_DAMAGE_BONUS_AURA_2);
	skillBuild.push_back(SKILL_MOVEMENT_BONUS_FACTOR_PASSIVE_1);
	skillBuild.push_back(SKILL_MOVEMENT_BONUS_FACTOR_AURA_1);
	skillBuild.push_back(SKILL_MOVEMENT_BONUS_FACTOR_PASSIVE_2);
	skillBuild.push_back(SKILL_MOVEMENT_BONUS_FACTOR_AURA_2);
}

#ifdef _zuko3d_output_path_stats
MyStrategy::~MyStrategy() {
	fout_path.close();
}
#endif // 