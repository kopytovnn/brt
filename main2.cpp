#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;


class SimpleDynamicBycicleModel {
private:
    struct x {
        float X;
        float Y;
        float yaw;
        float vx;
        float vy;
        float r;
    };

    struct u {
        float throttle;
        float steeringAngle;
        float brakes;
    };

    x old = { 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0 };
    u instant;
    float t = 0.0f;

    const float dt = 0.5f;
    const float m = 300.0f;
    const float lf = 0.721f;
    const float lr = 0.823f;
    const float Cm = 3600.0f;
    const float Crr = 200.0f;
    const float Cd = 1.53f;
    const float Cbf = 5411.0f;
    const float Cbr = 2650.0f;
    const float Cx = 20000.0f;
    const float Iz = 134.0f;
    float ar = 0.0f;
    float af = 0.0f;

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
        if (ar != ar) {
            return 0.0f;
        }
        return -Cx * ar;
    }

    float Ffy() {
        if (af != af) {
            return 0.0f;
        }
        return -Cx * af;
    }

    float Fx() {
        return Fdrv() - Frrr() - Frrf() * cos(instant.steeringAngle)
            - Fdrag() * cos(instant.steeringAngle) - Fbr()
            - Fbf() * cos(instant.steeringAngle);
    }

    float Fy() {
        return -Frrf() * sin(instant.steeringAngle) - Fdrag() * sin(instant.steeringAngle)
            - Fbf() * sin(instant.steeringAngle) + Fry() + Ffy() * cos(instant.steeringAngle);
    }

    float F() {
        return sqrt(Fx() * Fx() + Fy() * Fy());
    }

    float Rr() {  // Distance between ICR & rear wheel
        return (lr + lf) / tan(instant.steeringAngle);
    }

    float Rс() {  // Distance between ICR & center of mass
        return sqrt(Rr() * Rr() + lr * lr);
    }

    float Rf() {
        return (lr + lf) / sin(instant.steeringAngle);
    }

    float L() {  // Angular momentum
        return -lr * Fry() + lf * (Ffy() - Fbf() - Frrf()) * cos(instant.steeringAngle);
    }

public:
    float getX() const { return old.X; }
    float getY() const { return old.Y; }
    float getyaw() const { return old.yaw; }
    float getvx() const { return old.vx; }
    float getvy() const { return old.vy; }
    float getr() const { return old.r; }
    float gett() const { return t; }

    void update(float throttle, float steeringAngle, float brakes) {
        instant.throttle = throttle;
        instant.steeringAngle = steeringAngle;
        instant.brakes = brakes;

        float vrx = old.vx * Rr() / Rс();
        float vry = old.vy * Rr() / Rс();
        ar = atan(vry / vrx);

        float vfx = old.vx * Rf() / Rс();
        float vfy = old.vy * Rf() / Rс();
        float ve = vfx * cos(instant.steeringAngle) + vfy * sin(instant.steeringAngle);
        float vn = vfy * cos(instant.steeringAngle) - vfx * sin(instant.steeringAngle);
        af = atan(vn / ve);

        float angularAcceleration = L() / Iz;

        float vx_p1 = old.vx + instant.throttle * cos(instant.steeringAngle) * dt;
        float vy_p1 = old.vy + instant.throttle * sin(instant.steeringAngle) * dt;

        float X_p1 = old.X + (old.vx * cos(old.yaw) - old.vy * sin(old.yaw)) * dt;
        float Y_p1 = old.Y + (old.vx * sin(old.yaw) + old.vy * cos(old.yaw)) * dt;



        float yaw_p1 = old.yaw + old.r * dt;
        float r_p1 = old.r + angularAcceleration * dt;

        t += dt;
        old = { X_p1, Y_p1, yaw_p1, vx_p1, vy_p1, r_p1 };
    }
};


ostream& operator<<(ostream& s, const SimpleDynamicBycicleModel& sdbm) {
    return (s << "x(" << sdbm.gett() << ") = {" << sdbm.getX() <<
        " " << sdbm.getY() << " " << sdbm.getyaw() <<
        " " << sdbm.getvx() << " " << sdbm.getvy() << " " << sdbm.getr() << " }");
}


int main()
{
    //SimpleKinematicBycicleModel A;
    SimpleDynamicBycicleModel A;
    std::ifstream in("input.txt");
    if (in.is_open())
    {
        int n;
        in >> n;
        for (int i = 0; i < n; i++) {
            float a, sa, br;
            in >> a >> sa >> br;
            A.update(a, sa, br);
            cout << A << endl;
        }
    }
    in.close();
}
