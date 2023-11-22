#include <iostream>
#include <cmath>
#include <fstream>
#include "C16_CONTINENTAL_Tire_Data.h"
#include "constants.h"
using namespace C16_CONTINENTAL_Tire_Data;
using namespace CONSTANTS_JBS;
using namespace std;


struct state {
	double X, Y, yaw, vx, vy, w;

	state operator*=(float a) {
		return state{
			this->X * a,
			this->Y * a,
			this->yaw * a,
			this->vx * a,
			this->vy * a,
			this->w * a };
	}

	state operator+(state a) {
		return state{
			this->X + a.X,
			this->Y + a.Y,
			this->yaw + a.yaw,
			this->vx + a.vx,
			this->vy + a.vy,
			this->w + a.w };
	}
};

state operator*(state a, float b) {
	return a *= b;
}

state operator*(float a, state b) {
	return b *= a;
}

struct controlInfluence {
	double throttle, steeringAngle, brakes;
};


class DynamicBycicleModel {
private:
	state carState = { 0, 0, 0, 0, 0, 0 };
	double t = 0;

	//Forces of the longitudinal dynamics of the car
	double Fdrv(controlInfluence input) {
		return input.throttle * Cm;
	}
	double Frrr(state actual) {
		return Crr * tanh(actual.vx);
	}
	double Frrf(state actual) {
		return Crr * tanh(actual.vx);
	}
	double Fdrag(state actual) {
		return Cd * actual.vx * actual.vx;
	}
	double Fbf(controlInfluence input, state actual) {
		return input.brakes * Cbf * tanh(actual.vx);
	}
	double Fbr(controlInfluence input, state actual) {
		return input.brakes * Cbr * tanh(actual.vx);
	}
	//Forces of the lateral dynamics of the car
	double Fry(double ar) {
		return 2 * Cx * ar;
	}
	double Ffy(double af) {
		return 2 * Cx * af;
	}

	double Ftransversal(controlInfluence input, state actual, double af, double ar) {
		return Fdrv(input)
			- Frrr(actual)
			- Frrf(actual) * cos(input.steeringAngle)
			- Fdrag(actual)
			- Fbf(input, actual) * cos(input.steeringAngle)
			- Fbr(input, actual)
			- Ffy(af) * sin(input.steeringAngle);
	}

	double Flateral(controlInfluence input, state actual, double af, double ar) {
		return -Frrf(actual) * sin(input.steeringAngle)
			- Fbf(input, actual) * sin(input.steeringAngle)
			+ Fry(ar)
			+ Ffy(af) * cos(input.steeringAngle);
	}

	double L(controlInfluence input, state actual, double af, double ar) {
		return -Frrf(actual) * sin(input.steeringAngle) * lf
			- Fbf(input, actual) * sin(input.steeringAngle) * lf
			- Fry(ar) * lr
			+ Ffy(af) * cos(input.steeringAngle) * lf;
	}

	state Derivatives(controlInfluence input, state actual, double af, double ar) {
		double dxdt = actual.vx * cos(actual.yaw) - actual.vy * sin(actual.yaw);
		double dydt = actual.vx * sin(actual.yaw) + actual.vy * cos(actual.yaw);
		double dyawdt = actual.w;
		double dvxdt = (Ftransversal(input, actual, af, ar) / m + actual.vy * actual.w);
		double dvydt = (Flateral(input, actual, af, ar) / m - actual.vx * actual.w);
		//cout << L(input, actual, af, ar) << endl;
		double dwdt = L(input, actual, af, ar) / Iz;
		//cout << dwdt << "\n";

		return state{
			dxdt,
			dydt,
			dyawdt,
			dvxdt,
			dvydt,
			dwdt
		};
	}

public:
	float getX() const { return carState.X; }
	float getY() const { return carState.Y; }
	float getyaw() const { return carState.yaw; }
	float getvx() const { return carState.vx; }
	float getvy() const { return carState.vy; }
	float getr() const { return carState.w; }
	float gett() const { return t; }

	void updateEuler(controlInfluence input) {
		double af = 0;
		double ar = 0;
		float vrx = carState.vx;
		float vry = carState.vy - carState.w * lr;
		if (vrx == 0) {
			ar = 0;
		}
		else {
			ar = atan2((carState.vy - lr * carState.w), carState.vx);
		}
		float vfx = carState.vx;
		float vfy = carState.vy + carState.w * lf;
		float ve = vfx * cos(input.steeringAngle) + vfy * sin(input.steeringAngle);
		float vn = -vfx * sin(input.steeringAngle) + vfy * cos(input.steeringAngle);
		if (ve == 0) {
			af = 0;
		}
		else {
			af = atan2((carState.vy + lf * carState.w), carState.vx) - input.steeringAngle;;
		}
		t += dt;
		carState = carState + dt * Derivatives(input, carState, af, ar);
	}
};

ostream& operator<<(ostream& s, const DynamicBycicleModel& sdbm) {
	return (s << "x(" << sdbm.gett() << ") = {" << sdbm.getX() <<
		" " << sdbm.getY() << " " << sdbm.getyaw() <<
		" " << sdbm.getvx() << " " << sdbm.getvy() << " " << sdbm.getr() << " }");
}


int main()
{
	int iterations_by_one_step = 50;
	DynamicBycicleModel A;
	std::ifstream in("input.txt");
	if (in.is_open())
	{
		int n;
		in >> n;
		for (int i = 0; i < n; i++) {
			float a, sa, br;
			in >> a >> sa >> br;

			for (int j = 0; j < iterations_by_one_step; j++) {
				A.updateEuler(controlInfluence{ a, sa, br });
			}
			cout << A << endl;
		}
	}
	in.close();
}
