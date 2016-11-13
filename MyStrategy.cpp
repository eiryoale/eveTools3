#include "MyStrategy.h"

#define PI 3.14159265358979323846
#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdlib>
#include <queue>
#include <iostream>

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
	return 3.0;// !!!!!!!!!!!!!!!!!!! ÍÓÆÍÎ Ó×ÈÒÛÂÀÒÜ ÁÀÔÔÛ íà ñêîğîñòü
}

double maxMSY(const Wizard& w) {
	return 4.0;// !!!!!!!!!!!!!!!!!!! ÍÓÆÍÎ Ó×ÈÒÛÂÀÒÜ ÁÀÔÔÛ íà ñêîğîñòü
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
	
	double maxSpeedX = maxMSX(self); // !!!!!!!!!!!!!!!!!!! ÍÓÆÍÎ Ó×ÈÒÛÂÀÒÜ ÁÀÔÔÛ íà ñêîğîñòü
	double maxSpeedY = maxMSY(self); // !!!!!!!!!!!!!!!!!!! ÍÓÆÍÎ Ó×ÈÒÛÂÀÒÜ ÁÀÔÔÛ íà ñêîğîñòü

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
