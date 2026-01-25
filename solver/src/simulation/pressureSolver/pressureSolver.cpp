
#include "simulation/pressureSolver/pressureSolver.h"

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
        }
    }
}
