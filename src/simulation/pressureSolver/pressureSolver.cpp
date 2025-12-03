
#include "simulation/pressureSolver/pressureSolver.h"

PressureSolver::PressureSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning, double epsilon, int maximumNumberOfIterations) : 
    grid_(std::move(grid)),
    partitioning_(std::move(partitioning)),
    epsilon_(epsilon), 
    maxNumberOfIterations_(maximumNumberOfIterations)
{}

void PressureSolver::setBoundaryValues() {
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
