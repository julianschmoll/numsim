#include "redBlackSolver.h"
#include "grid/dataField.h"
#include "macros.h"
#include "simulation/partitioning.h"
#include <limits>
#include <vector>

RedBlackSolver::RedBlackSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning,
        double epsilon, int maximumNumberOfIterations, double omega) : 
    grid_(std::move(grid)),
    partitioning_(std::move(partitioning)), 
    epsilon_(epsilon), 
    maxNumberOfIterations_(maximumNumberOfIterations), 
    omega_(omega) 
{}

void RedBlackSolver::solve() {
    int it = 0;

    DataField &p = grid_->p();
    DataField &rhs = grid_->rhs();

    const double dx2 = grid_->dx() * grid_->dx();
    const double dy2 = grid_->dy() * grid_->dy();
    const double invDx2 = 1 / dx2;
    const double invDy2 = 1 / dy2;
    const double scalingFactor = 0.5 * dx2 * dy2 / (dx2 + dy2);

    globalPressureError_ = std::numeric_limits<double>::max();

    const int beginI = p.beginI() + 1;
    const int endI = p.endI() - 1;

    const int beginJ = p.beginJ() + 1;
    const int endJ = p.endJ() - 1;

    while (it < maxNumberOfIterations_) {
        localPressureError_ = 0.0;
        // rb = 0 is red, rb=1 is black pass
        #pragma omp parallel for reduction(+ : localPressureError_)
        for (int rb = 0; rb < 2; ++rb) {
            #pragma omp simd reduction(+ : localPressureError_)
            for (int j = beginJ; j < endJ; ++j) {
                // ((beginI + j + rb) & 1) calculates offset 0 or 1
                for (int i = beginI + ((beginI + j + rb) & 1); i < endI; i += 2) {
                    const double pDxx = (p(i - 1, j) + p(i + 1, j)) * invDx2;
                    const double pDyy = (p(i, j - 1) + p(i, j + 1)) * invDy2;
                    const double pNew = (1 - omega_) * p(i, j) + omega_ * scalingFactor * (pDxx + pDyy - rhs(i, j));
                    localPressureError_ += (p(i, j) - pNew) * (p(i, j) - pNew);
                    p(i, j) = pNew;
                }
            }
            setBoundaryValues();
            partitioning_->exchange(std::vector{&p});
        }

        MPI_Allreduce(&localPressureError_, &globalPressureError_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // DEBUG(std::cout << "Local pressure Error: " << localPressureError_ << ", global pressure error: " << globalPressureError_ << std::endl;)
        if (globalPressureError_ < epsilon_ * epsilon_) {
            break;
        }
        it++;
    }
}

void RedBlackSolver::setBoundaryValues() {
    DataField &p = grid_->p();

    if (partitioning_->ownContainsBoundary<Direction::Bottom>()) {
        for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
            p(i, p.beginJ()) = p(i, p.beginJ() + 1);
        }
    }
    if (partitioning_->ownContainsBoundary<Direction::Top>()) {
        for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
            p(i, p.endJ() - 1) = p(i, p.endJ() - 2);
        }
    }
    if (partitioning_->ownContainsBoundary<Direction::Left>()) {
        for (int j = p.beginJ(); j < p.endJ(); j++) {
            p(p.beginI(), j) = p(p.beginI() + 1, j);
        }
    }
    if (partitioning_->ownContainsBoundary<Direction::Right>()) {
        for (int j = p.beginJ(); j < p.endJ(); j++) {
            p(p.endI() - 1, j) = p(p.endI() - 2, j);
        }
    }
}
