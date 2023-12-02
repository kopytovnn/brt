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


class WheelModel{
public:
	double gamma, gammay, Fz0, dfz, muy, Fz, Cy, Dy, Ky, By, SHy, SVy, Cyk, Eyk, SHyk;
	double gammax, mux, Cx, Dx, Kx, Bx, SHx, Svx, Exa, Cxa, SHxa;
	double alpha, kappa;
	WheelModel(){
		gamma = 0.;
		gammay = gamma * LGAY;
		Fz0 = FNOMIN;
		Fz = m * 9.81 / 4;
		dfz = (Fz - Fz0 * LFZO) / (Fz0 * LFZO);
		muy = (PDY1 + PDY2 * dfz) * (1 - PDY3 * gammay * gammay) * LMUY;
		Cy = PCY1 * LCY;
		Dy = muy * Fz;
		Ky = PKY1 * Fz0 * sin(2 * atan2(Fz, (PKY2 * Fz0 * LFZO))) * (1 - PKY3 * abs(gammay)) * LFZO * LKY;
		By = Ky / (Cy * Dy);
		SHy = (PHY1 + PHY2 * dfz) * LHY + PHY3 * gammay;
		SVy = Fz * ((PVY1 + PVY2 * dfz) * LVY + (PVY3 + PVY4 * dfz) * gammay) * LMUY;
		Cyk = RCY1;
		Eyk = REY1 + REY2 * dfz;
		SHyk = RHY1 + RHY2 * dfz;

		gammax = gamma * LGAX;
		mux = (PDX1 + PDX2 * dfz) * (1 - PDX3 * gammax * gammax) * LMUX;
		Cx = PCX1 * LCX;
		Dx = mux * Fz;
		Kx = Fz * (PKX1 + PKX2 * dfz) * exp(PKX3 * dfz) * LKX;
		Bx = Kx / (Cx * Dx);
		SHx = (PHX1 + PHX2 * dfz) * LHX;
		Svx = Fz * (PVX1 + PVX2 * dfz) * LVX * LMUX;
	    Exa = REX1 + REX2 * REX2 * dfz;
		Cxa = RCX1;
		SHxa = RHX1;
	}
	double angleVelocity;
	WheelModel operator*(float a){
		angleVelocity *= a;
	}
	WheelModel operator+(WheelModel a){
		angleVelocity += a.angleVelocity;
	}
	double pureLateralSlip(double alpha){
		float alphay = alpha + SHy;
		float Ey = (PEY1 + PEY2 * dfz) * (1 - (PEY3 + PEY4 * gammay) * sgn(alphay)) * LEY;
		float Fy0 = Dy * sin(Cy * atan(By * alphay - Ey * (By * alphay - atan(By * alphay)))) + SVy;
		return Fy0;
	}
	double combinedLateralSlip(double alpha, double kappa){
		double Fy0 = pureLateralSlip(alpha);
		double kappas = kappa + SHyk;
		double Byk = RBY1 * cos(atan(RBY2 * (alpha - RBY3))) * LYKA;
		double Dyk = Fy0 / cos(Cyk * atan(Byk * SHyk - Eyk * (Byk * SHyk - atan(Byk * SHyk))));
		double DVyk = muy * Fz * (RVY1 + RVY2 * dfz + RVY3 * gamma) * cos(atan(RVY4 * alpha));
		double SVyk = DVyk * sin(RVY5 * atan(RVY6 * kappa)) * LVYKA;
		double Fy = Dyk * cos(Cyk * atan(Byk * kappas - Eyk * (Byk * kappas - atan(Byk * kappas)))) + SVyk;
		return Fy;
	}
	double pureLongitudinalSlip(double kappa){
		double kappax = kappa + SHx;
		double Ex = (PEX1 + PEX2 * dfz + PEX3 * dfz * dfz) * (1 - PEX4 * sgn(kappax)) * LEX;
		double Fx0 = Dx * sin(Cx * atan(Bx * kappax - Ex * (Bx * kappax - atan(Bx * kappax)))) + Svx;
		return Fx0;
	}
	double combinedLongitudinalSlip(double alpha, double kappa){
		double Fx0 = pureLongitudinalSlip(kappa);
		double Bxa = RBX1 * cos(atan(RBX2 * kappa)) * LXAL;
		double alphas = alpha + SHxa;
		double Dxa = Fx0 / cos(Cxa * atan(Bxa * SHxa - Exa * (Bxa * SHxa - atan(Bxa * SHxa))));
		double Fx = Dxa * cos(Cxa * atan(Bxa * alphas - Exa * (Bxa * alphas - atan(Bxa * alphas))));
		return Fx;	
	}
private:

};

class FrontWheel : WheelModel{
public:
	FrontWheel operator*=(float a){
		angleVelocity *= a;
	}
	FrontWheel operator+(FrontWheel a){
		angleVelocity += a.angleVelocity;
	}

	void coeffs(controlInfluence input, state actual){
		double vfx = actual.vx;
		double vfy = actual.vy + actual.w * lf;
		double a = atan2((actual.vy + lf * actual.w), actual.vx) - input.steeringAngle;
		alpha = a;

		double ve = vfx * cos(input.steeringAngle) + vfy * sin(input.steeringAngle);
		double k = (angleVelocity * UNLOADED_RADIUS - ve) / max(ve, vxmin);
		kappa = k;
	}

	double Flongitudinal(state parentCar, controlInfluence ci){
		coeffs(ci, parentCar);
		return combinedLongitudinalSlip(alpha, kappa);
	}
	double Flateral(state parentCar, controlInfluence ci){
		coeffs(ci, parentCar);
		return combinedLateralSlip(alpha, kappa);
	}

	double wheelAngleAcceleration(state parentCar, controlInfluence ci){
		double Fx = Flongitudinal(parentCar, ci);
		double Frrf = Crr * tanh(parentCar.vx);
		double Fbf = ci.brakes * Cbf * tanh(parentCar.vx);
		double frontWheelMomentum = UNLOADED_RADIUS * (-Fx - Frrf - Fbf) / 2.;
		return frontWheelMomentum / Iw;
	}

	FrontWheel wheelDerivate(state parentCar, controlInfluence ci){
		FrontWheel d = FrontWheel();
		d.angleVelocity = wheelAngleAcceleration(parentCar, ci);
		return d;
	}
};


class RearWheel : WheelModel{
public:
	RearWheel operator*=(float a){
		angleVelocity *= a;
	}
	RearWheel operator+(RearWheel a){
		angleVelocity += a.angleVelocity;
	}

	void coeffs(controlInfluence input, state actual){
		double vrx = actual.vx;
		double vry = actual.vy - actual.w * lr;
		double a = atan2((actual.vy - lr * actual.w), actual.vx);
		alpha = a;

		double vfx = actual.vx;
		double vfy = actual.vy + actual.w * lf;
		double ve = vfx * cos(input.steeringAngle) + vfy * sin(input.steeringAngle);
		double k = (angleVelocity * UNLOADED_RADIUS - vrx) / max(vrx, vxmin);
		kappa = k;
	}

	double Flongitudinal(state parentCar, controlInfluence ci){
		coeffs(ci, parentCar);
		return combinedLongitudinalSlip(alpha, kappa);
	}
	double Flateral(state parentCar, controlInfluence ci){
		coeffs(ci, parentCar);
		return combinedLateralSlip(alpha, kappa);
	}

	double wheelAngleAcceleration(state parentCar, controlInfluence ci){
		double Fx = Flongitudinal(parentCar, ci);
		double Fdrv = ci.throttle  * Cm;
		double Frrr = Crr * tanh(parentCar.vx);
		double Fbr = ci.brakes * Cbr * tanh(parentCar.vx);
		double frontWheelMomentum = UNLOADED_RADIUS * (Fdrv - Fx - Frrr - Fbr) / 2.;
		return frontWheelMomentum / Iw;
	}

	RearWheel wheelDerivate(state parentCar, controlInfluence ci){
		RearWheel d = RearWheel();
		d.angleVelocity = wheelAngleAcceleration(parentCar, ci);
		return d;
	}
};

FrontWheel operator*(FrontWheel a, double b){
	return a *= b;
}

FrontWheel operator*(double a, FrontWheel b){
	return b *= a;
}

RearWheel operator*(RearWheel a, double b){
	return a *= b;
}

RearWheel operator*(double a, RearWheel b){
	return b *= a;
}

struct state
{
	double X, Y, yaw, vx, vy, w; 
	FrontWheel front_left, front_right;
	RearWheel rear_left, rear_right;

	state operator*=(double a) {
		return state{
			X * a,
			Y * a,
			yaw * a,
			vx * a,
			vy * a,
			w * a,
			front_left * a,
			front_right * a, 
			rear_left * a,
			rear_right * a};
		}

	state operator+(state a) {
		return state{
			X + a.X,
			Y + a.Y,
			yaw + a.yaw,
			vx + a.vx,
			vy + a.vy,
			w + a.w,
			front_left + a.front_left,
			front_right + a.front_right, 
			rear_left + a.rear_left,
			rear_right + a.rear_left};
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

class Dynamic4WheelsModel{
public:
	state carState;
	double t = 0;
	Dynamic4WheelsModel(){
		FrontWheel front_left = FrontWheel();
		FrontWheel front_right = FrontWheel();
		RearWheel rear_left = RearWheel();
		RearWheel rear_right = RearWheel();
		carState = {0, 0, 0, 0, 0, 0, front_left, front_right, rear_left, rear_right};
	}


	double Flongitudinal(controlInfluence input, state actual){
		double steeringAngle = input.steeringAngle;

		double flWheel = actual.front_left.Flongitudinal(actual, input) * cos(steeringAngle)
						-actual.front_left.Flateral(actual, input) * sin(steeringAngle);
		
		double frWheel = actual.front_right.Flongitudinal(actual, input) * cos(steeringAngle)
				-actual.front_right.Flateral(actual, input) * sin(steeringAngle);

		double rlWheel = actual.rear_left.Flongitudinal(actual, input);
		
		double rrWheel = actual.rear_right.Flongitudinal(actual, input);

		return flWheel + frWheel + rlWheel + rrWheel;
	}

	double Flateral(controlInfluence input, state actual){
		double steeringAngle = input.steeringAngle;

		double flWheel = actual.front_left.Flongitudinal(actual, input) * sin(steeringAngle)
						+ actual.front_left.Flateral(actual, input) * cos(steeringAngle);
		
		double frWheel = actual.front_right.Flongitudinal(actual, input) * sin(steeringAngle)
				+ actual.front_right.Flateral(actual, input) * cos(steeringAngle);

		double rlWheel = actual.rear_left.Flateral(actual, input);
		
		double rrWheel = actual.rear_right.Flateral(actual, input);

		return flWheel + frWheel + rlWheel + rrWheel;
	}

	double L(controlInfluence input, state actual){
		double steeringAngle = input.steeringAngle;

		double flWheelT = actual.front_left.Flongitudinal(actual, input) * cos(steeringAngle)
						-actual.front_left.Flateral(actual, input) * sin(steeringAngle);
		double flTM = flWheelT * (- b);
		double flWheelL = actual.front_left.Flongitudinal(actual, input) * sin(steeringAngle)
						+ actual.front_left.Flateral(actual, input) * cos(steeringAngle);
		double flLM = flWheelL * lf;
		
		double frWheelT = actual.front_right.Flongitudinal(actual, input) * cos(steeringAngle)
				-actual.front_right.Flateral(actual, input) * sin(steeringAngle);
		double frTM = frWheelT * b;
		double frWheelL = actual.front_right.Flongitudinal(actual, input) * sin(steeringAngle)
				+ actual.front_right.Flateral(actual, input) * cos(steeringAngle);
		double frLM = frWheelL * lf;


		double rlWheelT = actual.rear_left.Flongitudinal(actual, input);
		double rlTM = rlWheelT * (- b);
		double rlWheelL = actual.rear_left.Flateral(actual, input);
		double rlLM = rlWheelL * (- lr);

		double rrWheelT = actual.rear_right.Flongitudinal(actual, input);
		double rrTM = rrWheelT * b;
		double rrWheelL = actual.rear_right.Flateral(actual, input);
		double rrLM = rrWheelL * (- lr);

		return flTM + flLM + frLM + frTM + rlTM + rlLM + rrTM + rrLM;
	}

	state Derivatives(controlInfluence input, state actual){
		double dxdt = actual.vx * cos(actual.yaw) - actual.vy * sin(actual.yaw);
		double dydt = actual.vx * sin(actual.yaw) + actual.vy * cos(actual.yaw);
		double dyawdt = actual.w;
		double dvxdt = (Flongitudinal(input, actual) / m + actual.vy * actual.w);
		double dvydt = (Flateral(input, actual) / m - actual.vx * actual.w);
		double dwdt = L(input, actual) / Iz;
		
		FrontWheel dFLdt = actual.front_left.wheelDerivate(actual, input);
		FrontWheel dFRdt = actual.front_right.wheelDerivate(actual, input);
		RearWheel dRLdt = actual.rear_left.wheelDerivate(actual, input);
		RearWheel dRRdt = actual.rear_right.wheelDerivate(actual, input);

		return state{
			dxdt,
			dydt,
			dyawdt,
			dvxdt,
			dvydt,
			dwdt,
			dFLdt,
			dFRdt,
			dRLdt,
			dRRdt
		};
	}

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

ostream& operator<<(ostream& s, const Dynamic4WheelsModel& sdbm) {
	return (s << "x(" << sdbm.gett() << ") = {" << sdbm.getX() <<
		" " << sdbm.getY() << " " << sdbm.getyaw() <<
		" " << sdbm.getvx() << " " << sdbm.getvy() << " " << sdbm.getr() << " }");
}


int main()
{
	int iterations_by_one_step = 50;
	Dynamic4WheelsModel A;
	Dynamic4WheelsModel aas = Dynamic4WheelsModel();
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
