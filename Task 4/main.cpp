#include <iostream>
#include <cmath>
#include <fstream>
#include "C16_CONTINENTAL_Tire_Data.h"
#include "constants.h"
using namespace C16_CONTINENTAL_Tire_Data;
using namespace CONSTANTS_JBS;
using namespace std;

float sgn(float v) {
	return (v > 0) - (v < 0);
}

struct state {
	double X, Y, yaw, vx, vy, w, omegaf, omegar;

	state operator*=(float a) {
		return state{
			this->X * a,
			this->Y * a,
			this->yaw * a,
			this->vx * a,
			this->vy * a,
			this->w * a,
			this->omegaf * a,
			this->omegar * a };
	}

	state operator+(state a) {
		return state{
			this->X + a.X,
			this->Y + a.Y,
			this->yaw + a.yaw,
			this->vx + a.vx,
			this->vy + a.vy,
			this->w + a.w,
			this->omegaf + a.omegaf,
			this->omegar + a.omegar };
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
	state carState = { 0, 0, 0, 0, 0, 0, 0, 0 };
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
	double FPacejkalateral(float alpha, float gamma, float Fz) {
		float gammay = gamma * LGAY;
		float Fz0 = FNOMIN;

		float dfz = (Fz - Fz0 * LFZO) / (Fz0 * LFZO);
		float muy = (PDY1 + PDY2 * dfz) * (1 - PDY3 * gammay * gammay) * LMUY;;
		float Cy = PCY1 * LCY;

		float Dy = muy * Fz;

		float Ky = PKY1 * Fz0 * sin(2 * atan2(Fz, (PKY2 * Fz0 * LFZO))) * (1 - PKY3 * abs(gammay)) * LFZO * LKY;
		float By = Ky / (Cy * Dy);

		float SHy = (PHY1 + PHY2 * dfz) * LHY + PHY3 * gammay;

		float alphay = alpha + SHy;
		float Ey = (PEY1 + PEY2 * dfz) * (1 - (PEY3 + PEY4 * gammay) * sgn(alphay)) * LEY;

		float SVy = Fz * ((PVY1 + PVY2 * dfz) * LVY + (PVY3 + PVY4 * dfz) * gammay) * LMUY;

		float Fy0 = Dy * sin(Cy * atan(By * alphay - Ey * (By * alphay - atan(By * alphay)))) + SVy;
		return Fy0;
	}
	double Fry(double ar) {
		return 2 * FPacejkalateral(ar, 0, m * 9.81 / 4);
	}
	double Ffy(double af) {
		return 2 * FPacejkalateral(af, 0, m * 9.81 / 4);
	}
	double FPacejkatransversal(float kappa, float gamma, float Fz) {
		float gammax = gamma * LGAX;
		float Fz0 = FNOMIN;

		float dfz = (Fz - Fz0 * LFZO) / (Fz0 * LFZO);
		float mux = (PDX1 + PDX2 * dfz) * (1 - PDX3 * gammax * gammax) * LMUX;;
		float Cx = PCX1 * LCX;

		float Dx = mux * Fz;

		float Kx = Fz0 * (PKX1 + PKX2 * dfz) * exp(PKX3 * dfz) * LKX;
		float Bx = Kx / (Cx * Dx);

		float SHx = (PHX1 + PHX2 * dfz) * LHX;

		float kappax = kappa + SHx;
		float Ex = (PEX1 + PEX2 * dfz + PEX3 * dfz * dfz) * (1 - PEX4 * sgn(kappax)) * LEX;

		float Svx = Fz * (PVX1 + PVX2 * dfz) * LVX * LMUX;

		float Fx0 = Dx * sin(Cx * atan(Bx * kappax - Ex * (Bx * kappax - atan(Bx * kappax)))) + Svx;
		return Fx0;
	}
	double Frx(double kappar) {
		return 2 * FPacejkatransversal(kappar, 0, m * 9.81 / 4);
	}
	double Ffx(double kappaf) {
		return 2 * FPacejkatransversal(kappaf, 0, m * 9.81 / 4);
	}

	double Ftransversal(controlInfluence input, state actual, double af, double ar, double kappaf, double kappar) {
		return Fdrv(input)
			- Frrr(actual)
			- Frrf(actual) * cos(input.steeringAngle)
			- Fdrag(actual)
			- Fbf(input, actual) * cos(input.steeringAngle)
			- Fbr(input, actual)
			- Ffy(af) * sin(input.steeringAngle)
			+ Ffx(kappaf) * cos(input.steeringAngle)
			+ Frx(kappar);
	}

	double Flateral(controlInfluence input, state actual, double af, double ar, double kappaf) {
		return -Frrf(actual) * sin(input.steeringAngle)
			- Fbf(input, actual) * sin(input.steeringAngle)
			+ Fry(ar)
			+ Ffy(af) * cos(input.steeringAngle)
			+ Ffx(kappaf) * sin(input.steeringAngle);
	}

	double L(controlInfluence input, state actual, double af, double ar, double kappaf) {
		return -Frrf(actual) * sin(input.steeringAngle) * lf
			- Fbf(input, actual) * sin(input.steeringAngle) * lf
			- Fry(ar) * lr
			+ Ffy(af) * cos(input.steeringAngle) * lf
			+ Ffx(kappaf) * sin(input.steeringAngle) * lf;
	}

	double frontWheelAngleAcceleration(controlInfluence input, state actual, double kappaf) {
		double frontWheelMomentum = UNLOADED_RADIUS * (0.5 * Ffx(kappaf) - Frrf(actual) - Fbf(input, actual));
		double epsilonwr = frontWheelMomentum / Iw;
		return epsilonwr;
	}
	double rearWheelAngleAcceleration(controlInfluence input, state actual, double kappar) {
		double rearWheelMomentum = UNLOADED_RADIUS * (Fdrv(input) - 0.5 * Frx(kappar) - Frrr(actual) - Fbr(input, actual));
		double epsilonwr = rearWheelMomentum / Iw;
		return epsilonwr;
	}

	state Derivatives(controlInfluence input, state actual) {
		double af = 0;
		double ar = 0;
		double vrx = actual.vx;
		double vry = actual.vy - actual.w * lr;
		if (vrx == 0) {
			ar = 0;
		}
		else {
			ar = atan2((actual.vy - lr * actual.w), actual.vx);
		}
		double vfx = actual.vx;
		double vfy = actual.vy + actual.w * lf;
		double ve = vfx * cos(input.steeringAngle) + vfy * sin(input.steeringAngle);
		double vn = -vfx * sin(input.steeringAngle) + vfy * cos(input.steeringAngle);
		if (ve == 0) {
			af = 0;
		}
		else {
			af = atan2((actual.vy + lf * actual.w), actual.vx) - input.steeringAngle;
		}

		double kappaf = (actual.omegaf * UNLOADED_RADIUS - ve) / max(ve, vxmin);
		double kappar = (actual.omegar * UNLOADED_RADIUS - vrx) / max(vrx, vxmin);


		double dxdt = actual.vx * cos(actual.yaw) - actual.vy * sin(actual.yaw);
		double dydt = actual.vx * sin(actual.yaw) + actual.vy * cos(actual.yaw);
		double dyawdt = actual.w;
		//cout << Ftransversal(input, actual, af, ar, kappaf, kappar) << endl;
		double dvxdt = (Ftransversal(input, actual, af, ar, kappaf, kappar) / m + actual.vy * actual.w);
		double dvydt = (Flateral(input, actual, af, ar, kappaf) / m - actual.vx * actual.w);
		//cout << L(input, actual, af, ar) << endl;
		double dwdt = L(input, actual, af, ar, kappaf) / Iz;
		//cout << dwdt << "\n";
		double domegafdt = frontWheelAngleAcceleration(input, actual, kappaf);
		double domegardt = rearWheelAngleAcceleration(input, actual, kappar);

		return state{
			dxdt,
			dydt,
			dyawdt,
			dvxdt,
			dvydt,
			dwdt,
			domegafdt,
			domegardt
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

	void updateRK4(controlInfluence input) {
		float h = dt;
		state k1 = Derivatives(input, carState);
		state k2 = Derivatives(input, carState + k1 * (h / 2));
		state k3 = Derivatives(input, carState + k2 * (h / 2));
		state k4 = Derivatives(input, carState + k3 * h);

		carState = carState + (k1 + k2 * 2 + k3 * 2 + k4) * (h / 6);
		//carState = carState + Derivatives(input, carState, af, ar, kappaf, kappar) * dt;
		t += dt;
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
				A.updateRK4(controlInfluence{ a, sa, br });
			}
			cout << A << endl;
		}
	}
	in.close();
}
