#include "outputWriter.h"

#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPointData.h>

OutputWriterParaview::OutputWriterParaview(
    std::shared_ptr<Discretization> discretization)
    : OutputWriter(discretization) {
  // Create a vtkWriter_
  vtkWriter_ = vtkSmartPointer<vtkXMLImageDataWriter>::New();
}

void OutputWriterParaview::writeFile(double currentTime) {
  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output_" << std::setw(4) << setfill('0') << fileNo_ << "."
           << vtkWriter_->GetDefaultFileExtension();

  // increment file no.
  fileNo_++;

  // assign the new file name to the output vtkWriter_
  vtkWriter_->SetFileName(fileName.str().c_str());

  // initialize data set that will be output to the file
  vtkSmartPointer<vtkImageData> dataSet = vtkSmartPointer<vtkImageData>::New();
  dataSet->SetOrigin(0, 0, 0);

  // set spacing of mesh
  const double dx = discretization_->meshWidth()[0];
  const double dy = discretization_->meshWidth()[1];
  const double dz = 1;
  dataSet->SetSpacing(dx, dy, dz);

  // set number of points in each dimension, 1 cell in z direction
  std::array<int, 2> nCells = discretization_->nCells();
  dataSet->SetDimensions(
      nCells[0] + 1, nCells[1] + 1,
      1); // we want to have points at each corner of each cell

  // add pressure field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arrayPressure = vtkDoubleArray::New();

  // the pressure is a scalar which means the number of components is 1
  arrayPressure->SetNumberOfComponents(1);

  // Set the number of pressure values and allocate memory for it. We already
  // know the number, it has to be the same as there are nodes in the mesh.
  arrayPressure->SetNumberOfTuples(dataSet->GetNumberOfPoints());

  arrayPressure->SetName("pressure");

  // loop over the nodes of the mesh and assign the interpolated p values in the
  // vtk data structure we only consider the cells that are the actual
  // computational domain, not the helper values in the "halo"

  int index = 0; // index for the vtk data structure, will be incremented in the
                 // inner loop
  for (int j = 0; j < nCells[1] + 1; j++) {
    for (int i = 0; i < nCells[0] + 1; i++, index++) {
      const double x = i * dx;
      const double y = j * dy;

      arrayPressure->SetValue(index, discretization_->p().interpolateAt(x, y));
    }
  }

  // now, we should have added as many values as there are points in the vtk
  // data structure
  assert(index == dataSet->GetNumberOfPoints());

  // add the field variable to the data set
  dataSet->GetPointData()->AddArray(arrayPressure);

  // add velocity field variable
  // ---------------------------
  vtkSmartPointer<vtkDoubleArray> arrayVelocity = vtkDoubleArray::New();

  // here we have two components (u,v), but ParaView will only allow vector
  // glyphs if we have an ℝ^3 vector, therefore we use a 3-dimensional vector
  // and set the 3rd component to zero
  arrayVelocity->SetNumberOfComponents(3);

  // set the number of values
  arrayVelocity->SetNumberOfTuples(dataSet->GetNumberOfPoints());

  arrayVelocity->SetName("velocity");

  // loop over the mesh where p is defined and assign the values in the vtk data
  // structure
  index = 0; // index for the vtk data structure
  for (int j = 0; j < nCells[1] + 1; j++) {
    const double y = j * dy;

    for (int i = 0; i < nCells[0] + 1; i++, index++) {
      const double x = i * dx;

      std::array<double, 3> velocityVector;
      velocityVector[0] = discretization_->u().interpolateAt(x, y);
      velocityVector[1] = discretization_->v().interpolateAt(x, y);
      velocityVector[2] = 0.0; // z-direction is 0

      arrayVelocity->SetTuple(index, velocityVector.data());
    }
  }
  // now, we should have added as many values as there are points in the vtk
  // data structure
  assert(index == dataSet->GetNumberOfPoints());

  // add the field variable to the data set
  dataSet->GetPointData()->AddArray(arrayVelocity);

  // add current time
  vtkSmartPointer<vtkDoubleArray> arrayTime = vtkDoubleArray::New();
  arrayTime->SetName("TIME");
  arrayTime->SetNumberOfTuples(1);
  arrayTime->SetTuple1(0, currentTime);
  dataSet->GetFieldData()->AddArray(arrayTime);

  // Remove unused memory
  dataSet->Squeeze();

  // Write the data
  vtkWriter_->SetInputData(dataSet);

  // vtkWriter_->SetDataModeToAscii();     // comment this in to get ascii text
  // files: those can be checked in an editor
  vtkWriter_->SetDataModeToBinary(); // set file mode to binary files: smaller
                                     // file sizes

  // finally write out the data
  vtkWriter_->Write();
}

void OutputWriterText::writeFile(double currentTime) {
  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/output_" << std::setw(4) << std::setfill('0') << fileNo_
           << ".txt";

  // increment file no.
  fileNo_++;

  // open file
  std::ofstream file(fileName.str().c_str());

  if (!file.is_open()) {
    std::cout << "Could not write to file \"" << fileName.str() << "\".";
    return;
  }

  // write time
  file << "t: " << currentTime << std::endl;

  // write mesh width
  file << "nCells: " << discretization_->nCells()[0] << "x"
       << discretization_->nCells()[1] << ", dx: " << discretization_->dx()
       << ", dy: " << discretization_->dy() << std::endl
       << std::endl;

  const int fieldWidth = 9; // number of characters to use for a single value

  // write u
  // ---------
  // write header lines
  file << "u (" << discretization_->u().size()[0] << "x"
       << discretization_->u().size()[1] << "): " << std::endl
       << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl
       << std::string(fieldWidth * (discretization_->u().size()[0] + 2) + 1,
                      '-')
       << std::endl;

  // write u values
  for (int j = discretization_->uJEnd() - 1; j >= discretization_->uJBegin();
       j--) {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd();
         i++) {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6)
           << discretization_->u(i, j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write v
  // ---------
  // write header lines
  file << "v (" << discretization_->v().size()[0] << "x"
       << discretization_->v().size()[1] << "): " << std::endl
       << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl
       << std::string(fieldWidth * (discretization_->v().size()[0] + 2) + 1,
                      '-')
       << std::endl;

  // write v values
  for (int j = discretization_->vJEnd() - 1; j >= discretization_->vJBegin();
       j--) {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd();
         i++) {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6)
           << discretization_->v(i, j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write p
  // ---------
  // write header lines
  file << "p (" << discretization_->p().size()[0] << "x"
       << discretization_->p().size()[1] << "): " << std::endl
       << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl
       << std::string(fieldWidth * (discretization_->p().size()[0] + 2) + 1,
                      '-')
       << std::endl;

  // write p values
  for (int j = discretization_->pJEnd() - 1; j >= discretization_->pJBegin();
       j--) {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd();
         i++) {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6)
           << discretization_->p(i, j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write f
  // ---------
  // write header lines
  file << "F (" << discretization_->u().size()[0] << "x"
       << discretization_->u().size()[1] << "): " << std::endl
       << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->uIBegin(); i < discretization_->uIEnd(); i++) {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl
       << std::string(fieldWidth * (discretization_->u().size()[0] + 2) + 1,
                      '-')
       << std::endl;

  // write f values
  for (int j = discretization_->uJEnd() - 1; j >= discretization_->uJBegin();
       j--) {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->uIBegin(); i < discretization_->uIEnd();
         i++) {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6)
           << discretization_->f(i, j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write g
  // ---------
  // write header lines
  file << "G (" << discretization_->v().size()[0] << "x"
       << discretization_->v().size()[1] << "): " << std::endl
       << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->vIBegin(); i < discretization_->vIEnd(); i++) {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl
       << std::string(fieldWidth * (discretization_->v().size()[0] + 2) + 1,
                      '-')
       << std::endl;

  // write g values
  for (int j = discretization_->vJEnd() - 1; j >= discretization_->vJBegin();
       j--) {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->vIBegin(); i < discretization_->vIEnd();
         i++) {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6)
           << discretization_->g(i, j);
    }
    file << std::endl;
  }
  file << std::endl;

  // write rhs
  // ---------
  // write header lines
  file << "rhs (" << discretization_->p().size()[0] << "x"
       << discretization_->p().size()[1] << "): " << std::endl
       << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl
       << std::string(fieldWidth * (discretization_->p().size()[0] + 2) + 1,
                      '-')
       << std::endl;

  // write rhs values
  for (int j = discretization_->pJEnd() - 1; j >= discretization_->pJBegin();
       j--) {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd();
         i++) {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6)
           << discretization_->rhs(i, j);
    }
    file << std::endl;
  }
  file << std::endl;
}

void OutputWriterText::writePressureFile() {
  // counter for files, counter value is part of the file name
  static int pressurefileNo = 0;

  // Assemble the filename
  std::stringstream fileName;
  fileName << "out/pressure_" << std::setw(4) << std::setfill('0')
           << pressurefileNo++ << ".txt";

  // open file
  std::ofstream file(fileName.str().c_str());

  if (!file.is_open()) {
    std::cout << "Could not write to file \"" << fileName.str() << "\".";
    return;
  }

  // write mesh width
  file << "nCells: " << discretization_->nCells()[0] << "x"
       << discretization_->nCells()[1] << ", dx: " << discretization_->dx()
       << ", dy: " << discretization_->dy() << std::endl
       << std::endl;

  const int fieldWidth = 9; // number of characters to use for a single value

  // write p
  // ---------
  // write header lines
  file << "p (" << discretization_->p().size()[0] << "x"
       << discretization_->p().size()[1] << "): " << std::endl
       << std::string(fieldWidth, ' ') << "|";
  for (int i = discretization_->pIBegin(); i < discretization_->pIEnd(); i++) {
    file << std::setw(fieldWidth) << i;
  }
  file << std::endl
       << std::string(fieldWidth * (discretization_->p().size()[0] + 2) + 1,
                      '-')
       << std::endl;

  // write p values
  for (int j = discretization_->pJEnd() - 1; j >= discretization_->pJBegin();
       j--) {
    file << std::setw(fieldWidth) << j << "|";
    for (int i = discretization_->pIBegin(); i < discretization_->pIEnd();
         i++) {
      file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6)
           << discretization_->p(i, j);
    }
    file << std::endl;
  }
  file << std::endl;
}

OutputWriter::OutputWriter(std::shared_ptr<Discretization> discretization)
    : discretization_(discretization), fileNo_(0) {
  // create "out" subdirectory if it does not yet exist
  int returnValue = system("mkdir -p out");
  if (returnValue != 0)
    std::cout << "Could not create subdirectory \"out\"." << std::endl;
}
