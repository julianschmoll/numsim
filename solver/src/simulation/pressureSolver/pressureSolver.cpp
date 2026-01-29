
#include "simulation/pressureSolver/pressureSolver.h"
#include <cassert>

PressureSolver::PressureSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning, const Settings &settings)
    : partitioning_(std::move(partitioning)), settings_(settings), grid_(std::move(grid)) {}

void PressureSolver::setBoundaryValues(DataField &p) {

    if (partitioning_->ownContainsBoundary<Direction::Bottom>()) {
        switch (settings_.boundaryBottom) {
            case BoundaryType::InflowNoSlip: {
                for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
                    p(i, p.beginJ()) = p(i, p.beginJ() + 1);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
                    p(i, p.beginJ()) = 0;
                }
                break;
            }
            case BoundaryType::Elastic: break; // TODO: setStructureBoundaries() macht das eigentlich
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Top>()) {
        switch (settings_.boundaryTop) {
            case BoundaryType::InflowNoSlip: {
                for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
                    p(i, p.endJ() - 1) = p(i, p.endJ() - 2);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int i = p.beginI() + 1; i < p.endI() - 1; i++) {
                    p(i, p.endJ() - 1) = 0;
                }
                break;
            }
            case BoundaryType::Elastic: break; // TODO: setStructureBoundaries() macht das eigentlich
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Left>()) {
        switch (settings_.boundaryLeft) {
            case BoundaryType::InflowNoSlip: {
                for (int j = p.beginJ(); j < p.endJ(); j++) {
                    p(p.beginI(), j) = p(p.beginI() + 1, j);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int j = p.beginJ(); j < p.endJ(); j++) {
                    p(p.beginI(), j) = 0;
                }
                break;
            }
            
            case BoundaryType::Elastic: assert(false);
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Right>()) {
        switch (settings_.boundaryRight) {
            case BoundaryType::InflowNoSlip: {
                for (int j = p.beginJ(); j < p.endJ(); j++) {
                    p(p.endI() - 1, j) = p(p.endI() - 2, j);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int j = p.beginJ(); j < p.endJ(); j++) {
                    p(p.endI() - 1, j) = 0;
                }
                break;
            }

            case BoundaryType::Elastic: assert(false);
        }
    }
}

void PressureSolver::setStructureBoundaries(DataField &p) { // Randwerte im inneren falls struktur rand berürt werden die anderen verwendet?
    for (int j = p.beginJ() + 1; j < p.endJ() - 1; ++j) { // TODO: p.endJ() - 1 oder p.endJ() - 2
        for (int i = p.beginI() + 1; i < p.endI() - 1; ++i) {
            if (grid_->isSolid(i, j)) {
                const bool leftIsFluid = grid_->isFluid(i - 1, j);
                const bool rightIsFluid = grid_->isFluid(i + 1, j);
                const bool bottomIsFluid = grid_->isFluid(i, j - 1);
                const bool topIsFluid = grid_->isFluid(i, j + 1);
                const int fluidCells = (int(bottomIsFluid) + int(topIsFluid) + (leftIsFluid) + int(rightIsFluid));
                if (fluidCells == 0) {
                    p(i, j) = 0;
                } else {
                    p(i, j) = int(bottomIsFluid) * p(i, j - 1) + int(topIsFluid) * p(i, j + 1) + int(leftIsFluid) * p(i - 1, j) + int(rightIsFluid) * p(i + 1, j);
                    p(i, j) /= fluidCells;
                }
            }
        }
    }
}