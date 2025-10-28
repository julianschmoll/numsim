#pragma once

#include "outputWriter.h"

/** Write *.vti files that can be viewed with ParaView.
 *  The mesh that can be visualized in ParaView corresponds to the mesh of the
 * computational domain. All values are given for the nodes of the mesh, i.e.,
 * the corners of each cell. This means, values will be interpolated because the
 * values are stored at positions given by the staggered grid.
 */
class outputWriterParaview final : public outputWriter {
public:
    //! constructor
    //! @param discretization shared pointer to the discretization object that
    //! will contain all the data to be written to the file
    explicit outputWriterParaview(const std::shared_ptr<discretization>& discretization);

    //! write current velocities to file, filename is output_<count>.vti
    void writeFile(double currentTime) override;

private:
    vtkSmartPointer<vtkXMLImageDataWriter>
    vtkWriter_; //< vtk writer to write ImageData
};
