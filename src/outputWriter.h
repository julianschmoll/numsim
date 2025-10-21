#pragma once

#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

#include <memory>

class OutputWriterText : public OutputWriter {
public:
  //! use constructor of base class
  using OutputWriter::OutputWriter;

  //! write current velocities to file, filename is output_<count>.txt
  void writeFile(double currentTime);

  //! write only current values of pressure to file, filename is
  //! pressure_<count>.txt
  void writePressureFile();
};

/** Write *.vti files that can be viewed with ParaView.
 *  The mesh that can be visualized in ParaView corresponds to the mesh of the
 * computational domain. All values are given for the nodes of the mesh, i.e.,
 * the corners of each cell. This means, values will be interpolated because the
 * values are stored at positions given by the staggered grid.
 */
class OutputWriterParaview : public OutputWriter {
public:
  //! constructor
  //! @param discretization shared pointer to the discretization object that
  //! will contain all the data to be written to the file
  OutputWriterParaview(std::shared_ptr<Discretization> discretization);

  //! write current velocities to file, filename is output_<count>.vti
  void writeFile(double currentTime);

private:
  vtkSmartPointer<vtkXMLImageDataWriter>
      vtkWriter_; //< vtk writer to write ImageData
};

/** Inteface class for writing simulation data output.
 */
class OutputWriter
{
public:
  //! constructor
  //! @param discretization shared pointer to the discretization object that will contain all the data to be written to the file
  OutputWriter(std::shared_ptr<Discretization> discretization);

  //! write current velocities to file, filename is output_<count>.vti
  virtual void writeFile(double currentTime) = 0;

protected:

  std::shared_ptr<Discretization> discretization_;  //< a shared pointer to the discretization which contains all data that will be written to the file
  int fileNo_;   //< a counter that increments for every file, this number is part of the file name of output files
};

