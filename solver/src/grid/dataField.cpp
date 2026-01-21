#include "grid/dataField.h"

#include "macros.h"

#include <cassert>
#include <cmath>
#include <mpi.h>

DataField::DataField()
    : meshWidth_({0, 0}), offset_({0, 0}), fieldID_(EMPTY_ID) {}

DataField::DataField(const std::array<int, 2> size, const std::array<double, 2> meshWidth, const std::array<double, 2> offset, int fieldID)
    : Array2d(size), meshWidth_(meshWidth), offset_(offset), fieldID_(fieldID) {
    const int numberX = endI() - beginI();
    const int numberY = endJ() - beginJ();

    MPI_Type_vector(numberY, 1, numberX, MPI_DOUBLE, &mpiColType_);
    MPI_Type_commit(&mpiColType_);
}

DataField::DataField(DataField &&other) noexcept
    : Array2d(std::move(other)), meshWidth_(other.meshWidth_), offset_(other.offset_), fieldID_(other.fieldID_), mpiColType_(other.mpiColType_) {
    other.mpiColType_ = MPI_DATATYPE_NULL;
}

DataField &DataField::operator=(DataField &&other) noexcept {
    if (this == &other) {
        return *this;
    }
    if (mpiColType_ != MPI_DATATYPE_NULL) {
        MPI_Type_free(&mpiColType_);
    }
    Array2d::operator=(std::move(other));
    meshWidth_ = other.meshWidth_;
    offset_ = other.offset_;
    fieldID_ = other.fieldID_;
    mpiColType_ = other.mpiColType_;
    other.mpiColType_ = MPI_DATATYPE_NULL;
    return *this;
}

DataField::DataField(const DataField &other)
    : Array2d(other),
      meshWidth_(other.meshWidth_),
      offset_(other.offset_),
      fieldID_(other.fieldID_) {
    if (other.mpiColType_ != MPI_DATATYPE_NULL) {
        MPI_Type_dup(other.mpiColType_, &mpiColType_);
    }
}

DataField &DataField::operator=(const DataField &other) noexcept {
    if (this == &other) {
        return *this;
    }

    if (mpiColType_ != MPI_DATATYPE_NULL) {
        MPI_Type_free(&mpiColType_);
        mpiColType_ = MPI_DATATYPE_NULL;
    }

    Array2d::operator=(other);
    meshWidth_ = other.meshWidth_;
    offset_ = other.offset_;
    fieldID_ = other.fieldID_;

    if (other.mpiColType_ != MPI_DATATYPE_NULL) {
        MPI_Type_dup(other.mpiColType_, &mpiColType_);
    } else {
        mpiColType_ = MPI_DATATYPE_NULL;
    }

    return *this;
}


DataField::~DataField() {
    if (mpiColType_ != MPI_DATATYPE_NULL) {
        MPI_Type_free(&mpiColType_);
    }
}

int DataField::rows() const {
    return size_[1];
}

int DataField::cols() const {
    return size_[0];
}

int DataField::beginJ() const {
    return -1;
}

int DataField::endJ() const {
    return size_[1] - 1; // size = cells + 2
}

int DataField::beginI() const {
    return -1;
}

int DataField::endI() const {
    return size_[0] - 1;
}

void DataField::setToZero() {
    for (auto &entry : data_)
        entry = 0;
}

int DataField::getID() const {
    return fieldID_;
}

double DataField::interpolateAt(double x, double y) const {
    x /= meshWidth_[0];
    y /= meshWidth_[1];

    // Als Test, Druck:

    x += offset_[0]; // p: 0.5, u: 0
    y += offset_[1]; // p: 0.5, u: 0.5

    assert(0 <= x && x < size_[0]);
    assert(0 <= y && y < size_[1]);

    const int iLower = static_cast<int>(std::floor(x));
    const int iUpper = static_cast<int>(std::ceil(x));
    const double alpha = x - iLower;

    const int jLower = static_cast<int>(std::floor(y));
    const int jUpper = static_cast<int>(std::ceil(y));
    const double beta = y - jLower;

    assert(0 <= iLower && iLower < size_[0]);
    assert(0 <= iUpper && iUpper < size_[0]);
    assert(0 <= jLower && jLower < size_[1]);
    assert(0 <= jUpper && jUpper < size_[1]);

    const double d11 = data_.at(jLower * size_[0] + iLower);
    const double d21 = data_.at(jLower * size_[0] + iUpper);
    const double d12 = data_.at(jUpper * size_[0] + iLower);
    const double d22 = data_.at(jUpper * size_[0] + iUpper);

    const double dx1 = (1 - alpha) * d11 + alpha * d21;
    const double dx2 = (1 - alpha) * d12 + alpha * d22;

    return (1 - beta) * dx1 + beta * dx2;
}

MPI_Datatype DataField::mpiColType() const {
    return mpiColType_;
}
