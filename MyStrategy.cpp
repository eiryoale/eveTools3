#include "MyStrategy.h"
#include <functional>
#define PI 3.14159265358979323846
#define _USE_MATH_DEFINES

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

using namespace model;
using namespace std;

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

template <class T> inline double sqr(T x) {
	return x * x;
}

template<>
struct hash<pair<int, int> >
	: public _Bitwise_hash<pair<int, int> >
{	// hash functor for bool
};

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
};

struct circleSetComarator {
	bool operator() (const circle & c1, const circle& c2) {
		return c1.id < c2.id;
	}
};

double staticMapWidth = 150.0; // �� ������, ��� ������ ������ �������� ������� (85 ��� ������ ��������, �� ������ ����)
double staticMapHeight = 150.0;
set<circle, circleSetComarator> prevStaticObjects;

vector< vector< vector< circle > > > staticMap;

class point : public Unit {
public:
	double x, y;
	point(double a, double b) : Unit(0, a, b, 0,0,0,FACTION_ACADEMY) {
		x = a;
		y = b;
	}
};

double maxMSX(const Wizard& w) {
	return 3.0;// !!!!!!!!!!!!!!!!!!! ����� ��������� ����� �� ��������
}

double maxMSY(const Wizard& w) {
	return 4.0;// !!!!!!!!!!!!!!!!!!! ����� ��������� ����� �� ��������
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

	if (A * c.x + B * c.y + C < c.radius) return true;
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

	if (A * c.getX() + B * c.getY() + C < c.getRadius()) return true;
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

bool checkCollisions(double ax, double ay, double bx, double by, const vector< vector< vector<circle> > >& objects) {
	typedef pair<int, int> pt;
	set<pt> cells;
	double A, B, C;
	A = ay - by;
	B = bx - ax;
	C = -(A * ax + B * ay);
	
	double sx = floor(ax / staticMapWidth);
	double sy = floor(ay / staticMapHeight);

	int i, j;
	i = (int)sx;

	cells.insert(pt((int)sx, (int)sy));
	sx *= staticMapWidth;
	sy *= staticMapHeight;

	double dx = staticMapWidth * (bx > ax ? 1 : -1);
	int di = (bx > ax ? 1 : -1);
	sx += dx;
	i+=di;
	while (sx * di < bx * di) {
		static double y = -(A * sx + C) / B;
		j = (int)floor(y / staticMapHeight);
		cells.insert(pt(i, j));
		cells.insert(pt(i - 1, j));
		i += di;
		sx += dx;
	}

	double dy = staticMapHeight * (by > ay ? 1 : -1);
	int dj = (by > ay ? 1 : -1);
	j = (int)sy / staticMapHeight;
	while (sy * dj < by * dj) {
		static double x = -(B * sy + C) / A;
		i = (int)floor(x / staticMapHeight);
		cells.insert(pt(i, j));
		cells.insert(pt(i, j - 1));
		j += dj;
		sy += dy;
	}

	point a(ax, ay), b(bx, by);
	for (const auto&c : cells) {
		for (const auto& obj : objects[c.first][c.second]) {
			if (checkCollision(a, b, obj)) return true;
		}
	}

	return false;
}

auto findPath(const point& start, const point& dest, double step = 20.0) {
	typedef pair<int, int> pt;
	typedef pair<double, pt > qtype;

	auto h = [](const pt&a, const pt&b) -> double {return sqrt((double) (sqr(a.first - b.first) + sqr(a.second - b.second))); };

	list<point> path;
		
	unordered_set<pt> closed, openSet;
	priority_queue < qtype, vector<qtype>, std::greater<qtype> > open;
	unordered_map<pt, pt> cameFrom;
	
	unordered_map<pt, double> gscore;

	pt pstart(0, 0), goal((int) floor((dest.getX() - start.getX()) / step), (int) floor((dest.getY() - start.getY())/ step)); // !!!!!!!!!!! goal ����� ���� �����������. ���� ��� ������� ����������� ��������� �������� �����.
	
	gscore[pstart] = 0.0;

	open.push(qtype(h(pstart, goal), pt(0, 0)));
	openSet.insert(pt(0, 0));

	bool good = false;
	qtype cur;
	while (!open.empty()) {
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
				
				int cellx = (int) floor(nx / staticMapWidth);
				int celly = (int) floor(ny / staticMapHeight);
				if (cellx < 0) {
					closed.insert(neib);
					continue;
				}
				if (celly < 0) {
					closed.insert(neib);
					continue;
				}

				bool inObstalce = checkCollisions(curx, cury, nx,ny, staticMap);
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

	if (!good) return path;
	pt p_cur = cur.second;
	double x = start.getX();
	double y = start.getY();
	path.clear();
	path.push_front(dest);
	while (p_cur != pstart) {
		path.push_front(point(x + p_cur.first * step, y + p_cur.second * step));
		p_cur = cameFrom[p_cur];
	}
	return path;
}

list <point> curWay;
double posX = 0, posY = 0;
void setMoveToPoint(const Wizard& self, Move& move, const Unit& pt) {
	double angle = self.getAngleTo(pt);
	cout << "angle = " << angle << endl;
	cout << "dx = " << self.getX() - posX << "\t dy = " << self.getY() - posY << endl;
	posX = self.getX();
	posY = self.getY();
	cout << self.getX() << ", " << self.getY() << endl;
	
	double maxSpeedX = maxMSX(self); // !!!!!!!!!!!!!!!!!!! ����� ��������� ����� �� ��������
	double maxSpeedY = maxMSY(self); // !!!!!!!!!!!!!!!!!!! ����� ��������� ����� �� ��������

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


void addObjectsToMap(vector< vector< vector< circle > > > &mp, const set<circle>& objects, double map_width, double map_height ) {
	for (const auto& c : objects) {
		//cout << "new obj: " << c.x << "\t" << c.y << "\t" << c.radius << endl;

		bool top = false;
		bool bot = false;
		bool left = false;
		bool right = false;

		int i = floor(c.x / map_width);
		int j = floor(c.y / map_height);

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

		int i = floor(c.x / map_width);
		int j = floor(c.y / map_height);

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

void MyStrategy::move(const Wizard& self, const World& world, const Game& game, Move& move) {
	bool inBattle = false;
	
	// ========================================================= Prepare map ===================================================================================================================================================================
	set<circle> staticObjects, newStaticObjects, killedStaticObjects;

	auto buildings = world.getBuildings();
	for (auto& b : buildings) staticObjects.insert(b);
	auto trees = world.getTrees();
	for (auto& t : trees) staticObjects.insert(t);

	
	set_difference(staticObjects.begin(), staticObjects.end(), prevStaticObjects.begin(), prevStaticObjects.end(), inserter(newStaticObjects, newStaticObjects.end()));
	set_difference(prevStaticObjects.begin(), prevStaticObjects.end(), staticObjects.begin(), staticObjects.end(), inserter(killedStaticObjects, killedStaticObjects.end()));

	addObjectsToMap(staticMap, newStaticObjects, staticMapWidth, staticMapHeight);
	removeObjectsFromMap(staticMap, killedStaticObjects, staticMapWidth, staticMapHeight);

	/*
	double r_2 = c.radius * c.radius;
	double leftx = c.x - c.radius;
	double rightx = c.x + c.radius;
	int s = (int) floor(leftx / staticMapWidth);
	int e = (int) floor(rightx / staticMapWidth);

	for (int i = s; i <= e; i++) {
	double x0 = i * staticMapWidth;

	double ytop = c.y - sqrt(r_2 - sqr(x0 - c.x));
	double ybot = c.y + sqrt(r_2 - sqr(x0 - c.x));

	int sy = (int)floor(ytop / staticMapHeight);
	int ey = (int)floor(ybot / staticMapHeight);

	for (int j = sy; j <= ey; j++) {
	staticMap[i][j].push_back(c);
	}

	//������� �����!
	}
	*/

	// ----------- remove killed objects
	// .............................................

	// ----------- deal with trees
	// ���� � ����� ���� ��������� ��� ������, ������� ������ ���� - ������ ��� �� ������ ��������.
	// .............................................

	// ========================================================= Movement ===================================================================================================================================================================
	{
		auto ts = clock();
		
		if (!curWay.empty()) {
			if (self.getDistanceTo(curWay.front()) < 1e-2) {
				curWay = findPath(point(self.getX(), self.getY()), point(1200, 1200));
				cout << "time spent for path: " << clock() - ts << endl;
			}
			if (!curWay.empty()) {
				setMoveToPoint(self, move, curWay.front());

				if (!inBattle) {
					move.setTurn(self.getAngleTo(curWay.front()));
				}
			}
		}
	}
}

MyStrategy::MyStrategy() { 
	staticMap.resize(5 + ceil(4000.0 / staticMapHeight));
	for (auto &row : staticMap) {
		row.resize(5 + ceil(4000.0 / staticMapWidth));
		for (auto& cell : row) {
			cell.reserve(100);
		}
	}

	curWay = findPath(point(100, 3700), point(1200, 1200));
}
