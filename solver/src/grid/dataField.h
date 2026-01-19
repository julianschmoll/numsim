#pragma once

#include "grid/array2d.h"
#include "mpi.h"
#include <iostream>

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
    explicit DataField(std::array<int, 2> size, std::array<double, 2> meshWidth, std::array<double, 2> offset = {0, 0}, int fieldID = 0);

    /**
     * Move constructor for DataField.
     *
     * Needed to handle MPI Datatype.
     *
     * @param other The DataField object to move resources from.
     */
    DataField(DataField &&other) noexcept;

    /**
     * Default constructor.
     */
    DataField();

    /**
     * Move assignment operator for DataField.
     *
     * Needed for MPI Datatype.
     *
     * @param other The DataField object to move resources from.
     * @return Reference to the current object.
     */
    DataField &operator=(DataField &&other) noexcept;

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

    /**
     * Sets all values on DataField to 0.
     */
    void setToZero();

    /**
     * Gets ID of the DataField.
     *
     * @return ID of the DataField.
     */
    int getID() const;

    /**
     * Gets number of rows.
     *
     * @return Number of rows of the DataField.
     */
    int rows() const;

    /**
     * Gets number of columns.
     *
     * @return Number of columns of the DataField.
     */
    int cols() const;

    /// MPI Datatype to exchange columns.
    MPI_Datatype mpiColType() const;

private:
    /// Grid Spacing
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

    /// ID of the DataField.
    int fieldID_;

    /// DataType to correctly exchange columns of DataField.
    MPI_Datatype mpiColType_ = MPI_DATATYPE_NULL;
};
