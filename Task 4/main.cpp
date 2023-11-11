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

struct y {
    float X;
    float Y;
    float yaw;
    float vx;
    float vy;
    float r; 

    y operator*(float a) {
        return y{ 
            this->X * a, 
            this->Y * a, 
            this->yaw * a, 
            this->vx * a, 
            this->vy * a,
            this->r * a};
    }

    y operator+(float a) {
        return y{
            this->X + a,
            this->Y + a,
            this->yaw + a,
            this->vx + a,
            this->vy + a,
            this->r + a };
    }

    y operator+(y a) {
        return y{
            this->X + a.X,
            this->Y + a.Y,
            this->yaw + a.yaw,
            this->vx + a.vx,
            this->vy + a.vy,
            this->r + a.r };
    }
};

struct u {
    float throttle;
    float steeringAngle;
    float brakes;

    u operator+(u a) {
        return u{
            this->throttle + a.throttle,
            this->steeringAngle + a.steeringAngle,
            this->brakes + a.brakes
        };
    }

    u operator+(float a) {
        return u{
            this->throttle + a,
            this->steeringAngle + a,
            this->brakes + a 
        };
    }

    u operator*(float a) {
        return u{
            this->throttle + a,
            this->steeringAngle + a,
            this->brakes + a
        };
    }
};

class DynamicBycicleModel {
private:
    y old = { 0, 0, 0, 0, 0, 0 };
    float ar = 0.0f;
    float af = 0.0f;
    float t = 0.0f;
    u instant;

    float Fdrv() {
        return instant.throttle * Cm * instant.throttle;
    }

    float Frrr() {
        return Crr * tanh(old.vx);
    }

    float Frrf() {
        return Frrr();
    }

    float Fdrag() {
        return Cd * old.vx * old.vx;
    }

    float Fbf() {
        return instant.brakes * Cbf * tanh(old.vx);
    }

    float Fbr() {
        return instant.brakes * Cbr * tanh(old.vx);
    }

    float Fry() {
        return magicFy(ar, 0, m / 4);
    }

    float Ffy() {
        return magicFy(af, 0, m / 4);
    }

    float magicFy(float alpha, float gamma, float Fz) {
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

    float L() {  // Angular momentum
        //cout << "\tFfy() = " << Ffy() << endl;
        //cout << "\tFbf() = " << Fbf() << endl;
        //cout << "\tFrrf() = " << Frrf() << endl;
        //cout << "\tFry() = " << Fry() << endl;
        //cout << "\tFrrf() = " << Frrf() << endl;
        return -lr * Fry() + lf * (Ffy() * cos(instant.steeringAngle) - Fbf() * sin(instant.steeringAngle) - Frrf() * sin(instant.steeringAngle));
    }

    float Flateral() {
        float k = cos(atan(old.vy / old.vx));
        float Faero = 0.0f;
        if (k == k) {
            Faero = Fdrag() * k;
        }
        return -Frrr() - Fbr() + Fdrv()
            - Faero
            - Fbf() * cos(instant.steeringAngle) - Frrf() * cos(instant.steeringAngle) - Ffy() * sin(instant.steeringAngle);
    }

    float Ftransversal() {
        float k = sin(atan(old.vy / old.vx));
        float Faero = 0.0f;
        if (k == k) {
            Faero = Fdrag() * k;
        }
        return Fry() - Faero - (Fbf() + Frrf()) * sin(instant.steeringAngle) + Ffy() * cos(instant.steeringAngle);
    }

    float Ftotal() {
        cout << "\tFlateral() = " << Flateral() << endl;
        cout << "\tFtransversal() = " << Ftransversal() << endl;
        return sqrt(pow(Flateral(), 2) + pow(Ftransversal(), 2));
    }

    y f(y data) {
        float X_intermidiate = old.vx * cos(old.yaw) + old.vy * sin(old.yaw);
        float Y_intermidiate = old.vx * sin(old.yaw) + old.vy * cos(old.yaw);
        float yaw_intermidiate = old.yaw + old.r;
        float vx_intermidiate = old.vx + (Flateral() / m + old.vy * old.r);
        float vy_intermidiate = (Ftransversal() / m - old.vx * old.r);
        float r_intermidiate = L() / Iz;

        return y{
        X_intermidiate,
        Y_intermidiate,
        yaw_intermidiate,
        vx_intermidiate,
        vy_intermidiate,
        r_intermidiate };
    }


public:
    float getX() const { return old.X; }
    float getY() const { return old.Y; }
    float getyaw() const { return old.yaw; }
    float getvx() const { return old.vx; }
    float getvy() const { return old.vy; }
    float getr() const { return old.r; }
    float gett() const { return t; }


    void update(u un) {
        instant.throttle = un.throttle;
        instant.steeringAngle = un.steeringAngle;
        instant.brakes = un.brakes;

        float h = dt;
        y k1 = f(old);
        y k2 = f(old + k1 * (h / 2));
        y k3 = f(old + k2 * (h / 2));
        y k4 = f(old + k3 * h);

        y n1 = old + (k1 + k2 * 2 + k3 * 2 + k4) * (h / 6);

        old = n1;

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
                A.update(u{a, sa, br});
            }
            cout << A << endl;
        }
    }
    in.close();
}