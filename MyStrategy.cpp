#include "MyStrategy.h"
#include <functional>
#define PI 3.14159265358979323846
#define _USE_MATH_DEFINES

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

	circle(const CircularUnit& u) : radius(u.getRadius()), x(u.getX()), y(u.getY()), id(u.getId()) { }

	bool operator<(const circle& c) const {
		return id < c.id;
	}
};

struct circleSetComarator {
	bool operator() (const circle & c1, const circle& c2) {
		return c1.id < c2.id;
	}
};

double staticMapWidth = 150.0; // не меньше, чем радиус самого большого объекта (85 при данных правилах, не счита€ базы)
double staticMapHeight = 150.0;
set<circle, circleSetComarator> prevStaticObjects;

vector< vector< vector< circle > > > staticMap;

class point : public Unit {
public:
	point(double a, double b) : Unit(0, a, b, 0,0,0,FACTION_ACADEMY) {
	}
};

double maxMSX(const Wizard& w) {
	return 3.0;// !!!!!!!!!!!!!!!!!!! Ќ”∆Ќќ ”„»“џ¬ј“№ Ѕј‘‘џ на скорость
}

double maxMSY(const Wizard& w) {
	return 4.0;// !!!!!!!!!!!!!!!!!!! Ќ”∆Ќќ ”„»“џ¬ј“№ Ѕј‘‘џ на скорость
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

auto findPath(const point& start, const point& dest) {
	typedef pair<int, int> pt;
	typedef pair<double, pt > qtype;

	auto h = [](const pt&a, const pt&b) -> double {return sqrt((double) (sqr(a.first - b.first) + sqr(a.second - b.second))); };

	stack<point> path;

	double step = 20.0;
	
	unordered_set<pt> closed, openSet;
	priority_queue < qtype, vector<qtype>, std::greater<qtype> > open;
	unordered_map<pt, pt> cameFrom;
	
	unordered_map<pt, double> gscore;

	pt pstart(0, 0), goal((int) floor((dest.getX() - start.getX()) / step), (int) floor((dest.getY() - start.getY())/ step)); // !!!!!!!!!!! goal может быть недостижима. Ќадо как минимум осматривать несколько соседних точек.
	
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
		for (int i = -1; i < 2; i++) {
			for (int j = -1; j < 2; j++) {
				if ((0 == j) && (i == 0)) continue;
				
				auto neib = pt(cur.second.first + i, cur.second.second + j);

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
	path.push(dest);
	while (p_cur != pstart) {
		path.push(point(x + p_cur.first * step, y + p_cur.second * step));
		p_cur = cameFrom[p_cur];
	}
	return path;
}

stack <point> curWay;
double posX = 0, posY = 0;
void setMoveToPoint(const Wizard& self, Move& move, const Unit& pt) {
	double angle = self.getAngleTo(pt);
	cout << "angle = " << angle << endl;
	cout << "dx = " << self.getX() - posX << "\t dy = " << self.getY() - posY << endl;
	posX = self.getX();
	posY = self.getY();
	cout << self.getX() << ", " << self.getY() << endl;
	
	double maxSpeedX = maxMSX(self); // !!!!!!!!!!!!!!!!!!! Ќ”∆Ќќ ”„»“џ¬ј“№ Ѕј‘‘џ на скорость
	double maxSpeedY = maxMSY(self); // !!!!!!!!!!!!!!!!!!! Ќ”∆Ќќ ”„»“џ¬ј“№ Ѕј‘‘џ на скорость

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


void MyStrategy::move(const Wizard& self, const World& world, const Game& game, Move& move) {
	bool inBattle = false;
	
	// ========================================================= Prepare map ===================================================================================================================================================================
	set<circle, circleSetComarator> staticObjects, newStaticObjects, killedStaticObjects;

	auto buildings = world.getBuildings();
	for (auto& b : buildings) staticObjects.insert(b);
	auto trees = world.getTrees();
	for (auto& t : trees) staticObjects.insert(t);

	
	set_difference(staticObjects.begin(), staticObjects.end(), prevStaticObjects.begin(), prevStaticObjects.end(), inserter(newStaticObjects, newStaticObjects.end()));
	set_difference(prevStaticObjects.begin(), prevStaticObjects.end(), staticObjects.begin(), staticObjects.end(), inserter(killedStaticObjects, killedStaticObjects.end()));

	// ----------- add new objects
	for (auto& c : newStaticObjects) {
		//cout << "new obj: " << c.x << "\t" << c.y << "\t" << c.radius << endl;

		bool top = false;
		bool bot = false;
		bool left = false;
		bool right = false;

		int i = floor(c.x / staticMapWidth);
		int j = floor(c.y / staticMapHeight);

		staticMap[i][j].push_back(c);

		if (j > (int)floor((c.y - c.radius) / staticMapHeight)) {
			top = true;
			staticMap[i][j - 1].push_back(c);
		}

		if (j < (int)floor((c.y + c.radius) / staticMapHeight)) {
			bot = true;
			staticMap[i][j + 1].push_back(c);
		}

		if (i > (int)floor((c.x - c.radius) / staticMapWidth)) {
			left = true;
			staticMap[i - 1][j].push_back(c);
		}

		if (i < (int)floor((c.x + c.radius) / staticMapWidth)) {
			right = true;
			staticMap[i + 1][j].push_back(c);
		}

		if (top && left) {
			double x = i * staticMapWidth;
			double y = j * staticMapHeight;
			if(sqr(x - c.x) + sqr(y - c.y) < c.radius) staticMap[i - 1][j - 1].push_back(c);
		}

		if (top && right) {
			double x = (i + 1) * staticMapWidth;
			double y = j * staticMapHeight;
			if (sqr(x - c.x) + sqr(y - c.y) < c.radius) staticMap[i + 1][j - 1].push_back(c);
		}

		if (bot && left) {
			double x = i * staticMapWidth;
			double y = (j + 1) * staticMapHeight;
			if (sqr(x - c.x) + sqr(y - c.y) < c.radius) staticMap[i - 1][j + 1].push_back(c);
		}

		if (bot && right) {
			double x = (i + 1) * staticMapWidth;
			double y = (j + 1) * staticMapHeight;
			if (sqr(x - c.x) + sqr(y - c.y) < c.radius) staticMap[i + 1][j + 1].push_back(c);
		}

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

			//крайние точки!
		}
		*/
	}
	
	// ----------- remove killed objects
	// .............................................

	// ----------- deal with trees
	// если в €вной зоне видимости нет дерева, которое должно быть - убрать его из списка деревьев.
	// .............................................

	// ========================================================= Movement ===================================================================================================================================================================
	{
		if (!curWay.empty()) {
			if (self.getDistanceTo(curWay.top()) < 1e-2) {
				curWay.pop();
			}
			if (!curWay.empty()) {
				setMoveToPoint(self, move, curWay.top());

				if (!inBattle) {
					move.setTurn(self.getAngleTo(curWay.top()));
				}
			}
		}
	}
}

MyStrategy::MyStrategy() { 
	staticMap.resize(ceil(4000.0 / staticMapHeight));
	for (auto &row : staticMap) {
		row.resize(ceil(4000.0 / staticMapWidth));
		for (auto& cell : row) {
			cell.reserve(100);
		}
	}

	curWay = findPath(point(100, 3700), point(2000, 2000));
}
