
#include "simulation/pressureSolver/pressureSolver.h"

PressureSolver::PressureSolver(std::shared_ptr<StaggeredGrid> grid, std::shared_ptr<Partitioning> partitioning, const Settings &settings)
    : partitioning_(std::move(partitioning)), settings_(settings), grid_(std::move(grid)) {}

void PressureSolver::setBoundaryValues() {
    DataField &p = grid_->p();

    const int beginI = grid_->beginI(p);
    const int beginJ = grid_->beginJ(p);
    const int endI = grid_->endI(p);
    const int endJ = grid_->endJ(p);

    if (partitioning_->ownContainsBoundary<Direction::Bottom>()) {
        switch (settings_.boundaryBottom) {
            case BoundaryType::InflowNoSlip: {
                for (int i = beginI + 1; i < endI - 1; i++) {
                    p(i, beginJ) = p(i, beginJ + 1);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int i = beginI + 1; i < endI - 1; i++) {
                    p(i, beginJ) = 0;
                }
                break;
            }
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Top>()) {
        switch (settings_.boundaryTop) {
            case BoundaryType::InflowNoSlip: {
                for (int i = beginI + 1; i < endI - 1; i++) {
                    p(i, endJ - 1) = p(i, endJ - 2);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int i = beginI + 1; i < endI - 1; i++) {
                    p(i, endJ - 1) = 0;
                }
                break;
            }
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Left>()) {
        switch (settings_.boundaryLeft) {
            case BoundaryType::InflowNoSlip: {
                for (int j = beginJ; j < endJ; j++) {
                    p(beginI, j) = p(beginI + 1, j);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int j = beginJ; j < endJ; j++) {
                    p(beginI, j) = 0;
                }
                break;
            }
        }
    }

    if (partitioning_->ownContainsBoundary<Direction::Right>()) {
        switch (settings_.boundaryRight) {
            case BoundaryType::InflowNoSlip: {
                for (int j = beginJ; j < endJ; j++) {
                    p(endI - 1, j) = p(endI - 2, j);
                }
                break;
            }

            case BoundaryType::Outflow: {
                for (int j = beginJ; j < endJ; j++) {
                    p(endI - 1, j) = 0;
                }
                break;
            }
        }
    }
}
