#include "grid.h"


inline double VelocityEnv::operator()(const int i) const {
    assert(-2 < i && i < 2);
    return u[1 + i];
}

inline double PressureEnv::operator()(const int i, const int j) const {
    assert(-2 < i && i < 2 && -2 < j && j < 2 && std::abs(i) + std::abs(j) < 2);  // TODO: test all cases
    return p[(j + 1) * 2 + i];
}
