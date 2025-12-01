#pragma once

#include "grid/array2d.h"
#include "mpi.h"

/**
 * @class DataField
 * @brief 2D Data representation containing method for interpolation.
 */
class DataField : public Array2d {
public:
    /**
     * Constructs DataField with grid size, offset, and index ranges.
     *
     * @param size Dimensions of the data field (number of cells in x and y).
     * @param meshWidth Physical size of the grid.
     * @param offset Offset from the lower-left cell corner to the data point location.
     */
     //TODO: should parameters be const references?
    explicit DataField(std::array<int, 2> size, std::array<double, 2> meshWidth, std::array<double, 2> offset, int fieldID = 0);

    DataField();

    /**
     * Destructor for DataField.
     */
    ~DataField() override;

    /**
     * Interpolates the value at a given coordinate (x, y).
     *
     * Uses bilinear interpolation based on mesh width and offset.
     *
     * @param x x-coordinate.
     * @param y y-coordinate.
     * @return Interpolated value at (x, y).
     */
    [[nodiscard]] double interpolateAt(double x, double y) const;

    /**
     * Returns the starting index in the j-direction.
     *
     * @return First valid j-index.
     */
    int beginJ() const;

    /**
     * Returns the ending index in the j-direction.
     *
     * @return Last valid j-index.
     */
    int endJ() const;

    /**
     * Returns the starting index in the i-direction.
     *
     * @return First valid i-index.
     */
    int beginI() const;

    /**
     * Returns the ending index in the i-direction.
     *
     * @return Last valid i-index.
     */
    int endI() const;

    void setToZero();

    int getID() const;

    MPI_Datatype getMPIColType() const;
    MPI_Datatype getMPIRowType() const;

    int rows() const;
    int cols() const;

private:
    // Grid Spacing
    std::array<double, 2> meshWidth_;

    /**
     * Describes the offset of the position of the data point from the lower left cell edge.
     *
     * (0,1)---(1,1)
     *   |       |
     *   |       |
     * (0,0)---(1,0)
     *
     * - Pressure:                0.5, 0.5
     * - Velocity in x direction: 0,   0.5
     * - Velocity in y direction: 0.5, 0
     */
    std::array<double, 2> offset_;

    MPI_Datatype mpiColType_;

    int fieldID_;
};
