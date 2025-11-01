#pragma once

#include "outputWriter.h"

class OutputWriterText final : public OutputWriter {
public:
    //! use constructor of base class
    using OutputWriter::OutputWriter;

    //! write current velocities to file, filename is output_<count>.txt
    void writeFile(double currentTime) override;

    //! write only current values of pressure to file, filename is
    //! pressure_<count>.txt
    void writePressureFile();
};
