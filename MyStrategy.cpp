#include "MyStrategy.h"
#include <functional>
#define PI 3.14159265358979323846
#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdlib>
#include <queue>
#include <iostream>
#include <set>
#include <map>
#include <stack>

using namespace model;
using namespace std;

template <class T> inline double sqr(T x) {
	return x * x;
}

class point : public Unit {
public:
	point(double a, double b) : Unit(0, a, b, 0,0,0,FACTION_ACADEMY) {
	}
};

double maxMSX(const Wizard& w) {
	return 3.0;// !!!!!!!!!!!!!!!!!!! НУЖНО УЧИТЫВАТЬ БАФФЫ на скорость
}

double maxMSY(const Wizard& w) {
	return 4.0;// !!!!!!!!!!!!!!!!!!! НУЖНО УЧИТЫВАТЬ БАФФЫ на скорость
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

template <class c1, class c2> class default_map : map<c1, c2> {
	c2 default_value;

public:
	default_map(const c2& val) : default_map(val) { }

	c2& operator[](const c1& key) {
		if (this->find(key) == this->end()) return default_value;
		return this->operator[](key);
	}
};

auto findPath(const point& start, const point& dest) {
	auto h = [](const point&a, const point&b) -> double {return a.getDistanceTo(b); };

	stack<point> path;

	double step = 10.0;
	set<pair<int, int> > closed;
	
	typedef pair<int, int> pt;
	typedef pair<double, pt > qtype;

	priority_queue < qtype, vector<qtype>, std::greater<qtype> > open;
	map<pt, pt> cameFrom;
	
	default_map<pt, double> gscore(1e100), fscore(1e100);

	pt pstart(0, 0), goal((int) floor(dest.getX() - start.getX()), (int) floor(dest.getY() - start.getY())); // !!!!!!!!!!! goal может быть недостижима. Надо как минимум осматривать несколько соседних точек.
	
	gscore[pstart] = 0.0;
	fscore[pstart] = h(start, dest);

	open.push(qtype(h(start, dest), pt(0, 0)));

	bool good = false;
	qtype cur;
	while (!open.empty()) {
		cur = open.top();
		if (cur.second == goal) {
			good = true;
			break;
		}
		
	}

	if (!good) return path;
	pt p_cur = cur.second;
	double x = start.getX();
	double y = start.getX();
	path.push(dest);
	while (p_cur != pstart) {
		path.push(point(x + p_cur.first * step, y + p_cur.second * step));
		p_cur = cameFrom[p_cur];
	}
	return path;
}

queue <Unit> curWay;
double posX = 0, posY = 0;
void setMoveToPoint(const Wizard& self, Move& move, const Unit& pt) {
	double angle = self.getAngleTo(pt);
	cout << "angle = " << angle << endl;
	cout << "dx = " << self.getX() - posX << "\t dy = " << self.getY() - posY << endl;
	posX = self.getX();
	posY = self.getY();
	cout << self.getX() << ", " << self.getY() << endl;
	
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


void MyStrategy::move(const Wizard& self, const World& world, const Game& game, Move& move) {
	bool inBattle = false;

	// ========================================================= Movement ===================================================================================================================================================================
	{
		if (!curWay.empty()) {
			if (self.getDistanceTo(curWay.front()) < 1e-2) {
				curWay.pop();
			}
			if (!curWay.empty()) {
				setMoveToPoint(self, move, curWay.back());

				if (!inBattle) {
					move.setTurn(self.getAngleTo(curWay.front()));
				}
			}
		}
	}
}

MyStrategy::MyStrategy() { 
	curWay.push(point(90, 3950));
}
