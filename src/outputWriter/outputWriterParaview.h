#pragma once

#include "outputWriter.h"

/** Write *.vti files that can be viewed with ParaView.
 *  The mesh that can be visualized in ParaView corresponds to the mesh of the computational domain.
 *  All values are given for the nodes of the mesh, i.e., the corners of each cell.
 *  This means, values will be interpolated because the values are stored at positions given by the staggered grid.
 */
class OutputWriterParaview : public OutputWriter {
public:
    //! constructor
    //! @param grid shared pointer to the discretization object that
    //! will contain all the data to be written to the file
    OutputWriterParaview(std::shared_ptr<StaggeredGrid> grid, const Partitioning &partitioning);

    //! write current velocities to file, filename is output_<count>.vti
    void writeFile(double currentTime) override;

private:
    vtkSmartPointer<vtkXMLImageDataWriter> vtkWriter_; //< vtk writer to write ImageData
};
