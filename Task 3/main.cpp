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
    float X;
    float Y;
    float yaw;
    float vx;
    float vy;
    float r;

    state operator*=(float a) {
        return state{
            this->X * a,
            this->Y * a,
            this->yaw * a,
            this->vx * a,
            this->vy * a,
            this->r * a };
    }

    state operator+(state a) {
        return state{
            this->X + a.X,
            this->Y + a.Y,
            this->yaw + a.yaw,
            this->vx + a.vx,
            this->vy + a.vy,
            this->r + a.r };
    }
};

state operator*(state a, float b) {
    return a *= b;
}

state operator*(float a, state b) {
    return b *= a;
}

struct input {
    float throttle;
    float steeringAngle;
    float brakes;
};


class DynamicBycicleModel3 {
private:
    state old = { 0, 0, 0, 0, 0, 0 };
    float ar = 0.0f;
    float af = 0.0f;
    float kappa = 0;
    float t = 0.0f;
    input instant;

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

    float magicFnu(float alpha, float gamma, float Fz) {
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

    float Ffy() {
        return magicFnu(af, 0, m / 4);
        //return - Cx * alpha;
    }

    float Fry() {
        return magicFnu(ar, 0, m / 4);
        //return - Cx * alpha;
    }

    float L(int toOutput=0) {
        if (toOutput == 1) {
            cout << Ffy() << '\t' << Frrf() << '\t' << Fbf() << endl;
        }
        return Ffy() * cos(instant.steeringAngle) * lf
            - Frrf() * sin(instant.steeringAngle) * lf
            - Fbf() * sin(instant.steeringAngle) * lf
            - Fry() * lr;
    }

    float Flateral() {
        return -Ffy() * cos(instant.steeringAngle)
            - Fbf() * sin(instant.steeringAngle)
            - Frrf() * sin(instant.steeringAngle)
            + Fry();
    }

    float Ftransversal() {
        return -Ffy() * sin(instant.steeringAngle)
            - Fbf() * cos(instant.steeringAngle)
            - Frrf() * cos(instant.steeringAngle)
            - Fdrag()
            + Fdrv()
            - Fbr();

    }

    state f(state data) {
        float X_intermidiate = data.vx * cos(data.yaw) - data.vy * sin(data.yaw);
        float Y_intermidiate = data.vx * sin(data.yaw) + data.vy * cos(data.yaw);
        float yaw_intermidiate = data.r;
        float vx_intermidiate = (Ftransversal() / m + data.vy * data.r);
        float vy_intermidiate = (Flateral() / m - data.vx * data.r);
        float r_intermidiate = L() / Iz;

        return state{
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

    void update(input un) {
        instant.throttle = un.throttle;
        instant.steeringAngle = un.steeringAngle;
        instant.brakes = un.brakes;

        float vrx = old.vx;
        float vry = old.vy - old.r * lr;
        if (vrx == 0) {
            ar = 0;
        }
        else {
            ar = vry / vrx;
        }

        float vfx = old.vx;
        float vfy = old.vy + old.r * lf;
        float ve = vfx * cos(instant.steeringAngle) + vfy * sin(instant.steeringAngle);
        float vn = -vfx * sin(instant.steeringAngle) + vfy * cos(instant.steeringAngle);
        if (ve == 0) {
            af = 0;
        }
        else {
            af = vn / ve;
        }

        cout << "\tL: " << L() << "\taf: " << af << endl;
        cout << "\t\tvx: " << old.vx << "\tvy: " << old.vy << endl;

        float h = dt;
        state k1 = f(old);
        state k2 = f(old + k1 * (h / 2));
        state k3 = f(old + k2 * (h / 2));
        state k4 = f(old + k3 * h);

        state n1 = old + (k1 + k2 * 2 + k3 * 2 + k4) * (h / 6);

        //state n1 = old + f(old) * dt;
        old = n1;

        t += dt;
    }
};

ostream& operator<<(ostream& s, const DynamicBycicleModel3& sdbm) {
    return (s << "x(" << sdbm.gett() << ") = {" << sdbm.getX() <<
        " " << sdbm.getY() << " " << sdbm.getyaw() <<
        " " << sdbm.getvx() << " " << sdbm.getvy() << " " << sdbm.getr() << " }");
}

int main()
{
    int iterations_by_one_step = 50;
    DynamicBycicleModel3 A;
    std::ifstream in("input.txt");
    if (in.is_open())
    {
        int n;
        in >> n;
        for (int i = 0; i < n; i++) {
            float a, sa, br;
            in >> a >> sa >> br;

            for (int j = 0; j < iterations_by_one_step; j++) {
                A.update(input{ a, sa, br });
            }
            cout << A << endl;
        }
    }
    in.close();
}
