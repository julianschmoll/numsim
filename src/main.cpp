#include <iostream>

#ifdef VTK_AVAILABLE
#include <vtkMath.h>
#endif

int main() {
    std::cout << "Hi" << std::endl;

#ifdef VTK_AVAILABLE
    double p0[3] = {0.0, 0.0, 0.0};
    double p1[3] = {1.0, 1.0, 1.0};
    double squaredDistance = vtkMath::Distance2BetweenPoints(p0, p1);
    std::cout << "vtkMath::Distance2BetweenPoints(p0, p1) = " << squaredDistance << std::endl;
#endif
}