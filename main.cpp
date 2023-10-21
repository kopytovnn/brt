#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;


class SimpleKinematicBycicleModel {
private:
    struct x {
        float X;
        float Y;
        float yaw;
        float v;
    };

    x old = {0.0, 0.0, 0.0 , 0.0};
    float t = 0.0f;

    const float L = 1.5f;
    const float dt = 0.5f;

public:
    float getX() const { return old.X; }
    float getY() const { return old.Y; }
    float getyaw() const { return old.yaw; }
    float getv() const { return old.v; }
    float gett() const { return t; }

    void upd(float acceleration, float steeringAngle) {
        float v_p1 = old.v + dt * acceleration;
        float yaw_p1 = old.yaw + dt * (old.v * tan(steeringAngle) / L);
        float X_p1 = old.X + dt * (old.v * cos(old.yaw));
        float Y_p1 = old.Y + dt * (old.v * sin(old.yaw));
        
        t += dt;
        old = { X_p1, Y_p1, yaw_p1, v_p1 };
    }
};

ostream& operator<<(ostream& s, const SimpleKinematicBycicleModel& skbm) {
    return (s << "x(" << skbm.gett() << ") = {" << skbm.getX() << 
        " " << skbm.getY() << " " << skbm.getyaw() << 
        " " << skbm.getv() << " }");
}


int main()
{
    SimpleKinematicBycicleModel A;
    std::ifstream in("input.txt");
    if (in.is_open())
    {
        int n;
        in >> n;
        for (int i = 0; i < n; i++) {
            float a, sa;
            in >> a >> sa;
            A.upd(a, sa);
            cout << A << endl;
        }
    }
    in.close();
}
