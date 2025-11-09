#pragma once

#include "grid/staggeredGrid.h"

#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

#include <memory>

/** Inteface class for writing simulation data output.
 */
class OutputWriter {
public:
    virtual ~OutputWriter() = default;
    //! constructor
    //! @param grid shared pointer to the discretization object that will contain all the
    //! data to be written to the file
    explicit OutputWriter(const std::shared_ptr<StaggeredGrid> grid);

    //! write current velocities to file, filename is output_<count>.vti
    virtual void writeFile(double currentTime) = 0;

protected:
    const std::shared_ptr<StaggeredGrid> grid_; //< a shared pointer to the discretization which contains all data that
                                                // will be written to the file
    int fileNo_;                                //< a counter that increments for every file, this number is part of the file name
                                                // of output files
};
