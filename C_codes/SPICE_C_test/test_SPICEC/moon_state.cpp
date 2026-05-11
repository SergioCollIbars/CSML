#include <iostream>
#include <iomanip>

extern "C" {
#include "SpiceUsr.h"
}

using namespace std;

int main()
{
    SpiceDouble et;
    SpiceDouble state[6];
    SpiceDouble lt;

    furnsh_c("/Users/sergiocollibars/Documents/MATLAB/kernels/naif0012.tls");
    furnsh_c("/Users/sergiocollibars/Documents/MATLAB/kernels/de421.bsp");

    str2et_c("2015 APR 29 12:00:00 UTC", &et);

    spkezr_c("MOON", et, "J2000", "NONE", "EARTH", state, &lt);

    cout << fixed << setprecision(6);
    cout << "Moon state relative to Earth in J2000\n";
    cout << "Position [km]: "
         << state[0] << " "
         << state[1] << " "
         << state[2] << "\n";

    cout << setprecision(9);
    cout << "Velocity [km/s]: "
         << state[3] << " "
         << state[4] << " "
         << state[5] << "\n";

    kclear_c();

    return 0;
}