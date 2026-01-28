#include "outputWriterTextParallel.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

void OutputWriterTextParallel::writeFile(double currentTime) {
    // Assemble the filename
    std::stringstream fileName;
    fileName << folderName_ << "/output_" << std::setw(4) << std::setfill('0') << fileNo_ << "." << partitioning_.ownRank() << ".txt";

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
    file << "nCells: " << grid_->nCells()[0] << "x" << grid_->nCells()[1] << ", dx: " << grid_->dx() << ", dy: " << grid_->dy() << std::endl
         << std::endl;

    const int fieldWidth = 9; // number of characters to use for a single value

    // write u
    // ---------
    // write header lines
    file << "u (" << grid_->u().size()[0] << "x" << grid_->u().size()[1] << "): " << std::endl << std::string(fieldWidth, ' ') << "|";
    for (int i = grid_->u().beginI(); i < grid_->u().endI(); i++) {
        file << std::setw(fieldWidth) << i;
    }
    file << std::endl << std::string(fieldWidth * (grid_->u().size()[0] + 2) + 1, '-') << std::endl;

    // write u values
    for (int j = grid_->u().endJ() - 1; j >= grid_->u().beginJ(); j--) {
        file << std::setw(fieldWidth) << j << "|";
        for (int i = grid_->u().beginI(); i < grid_->u().endI(); i++) {
            file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << grid_->u(i, j);
        }
        file << std::endl;
    }
    file << std::endl;

    // write v
    // ---------
    // write header lines
    file << "v (" << grid_->v().size()[0] << "x" << grid_->v().size()[1] << "): " << std::endl << std::string(fieldWidth, ' ') << "|";
    for (int i = grid_->v().beginI(); i < grid_->v().endI(); i++) {
        file << std::setw(fieldWidth) << i;
    }
    file << std::endl << std::string(fieldWidth * (grid_->v().size()[0] + 2) + 1, '-') << std::endl;

    // write v values
    for (int j = grid_->v().endJ() - 1; j >= grid_->v().beginJ(); j--) {
        file << std::setw(fieldWidth) << j << "|";
        for (int i = grid_->v().beginI(); i < grid_->v().endI(); i++) {
            file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << grid_->v(i, j);
        }
        file << std::endl;
    }
    file << std::endl;

    // write p
    // ---------
    // write header lines
    file << "p (" << grid_->p().size()[0] << "x" << grid_->p().size()[1] << "): " << std::endl << std::string(fieldWidth, ' ') << "|";
    for (int i = grid_->p().beginI(); i < grid_->p().endI(); i++) {
        file << std::setw(fieldWidth) << i;
    }
    file << std::endl << std::string(fieldWidth * (grid_->p().size()[0] + 2) + 1, '-') << std::endl;

    // write p values
    for (int j = grid_->p().endJ() - 1; j >= grid_->p().beginJ(); j--) {
        file << std::setw(fieldWidth) << j << "|";
        for (int i = grid_->p().beginI(); i < grid_->p().endI(); i++) {
            file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << grid_->p(i, j);
        }
        file << std::endl;
    }
    file << std::endl;

    // write f
    // ---------
    // write header lines
    file << "F (" << grid_->u().size()[0] << "x" << grid_->u().size()[1] << "): " << std::endl << std::string(fieldWidth, ' ') << "|";
    for (int i = grid_->u().beginI(); i < grid_->u().endI(); i++) {
        file << std::setw(fieldWidth) << i;
    }
    file << std::endl << std::string(fieldWidth * (grid_->u().size()[0] + 2) + 1, '-') << std::endl;

    // write f values
    for (int j = grid_->u().endJ() - 1; j >= grid_->u().beginJ(); j--) {
        file << std::setw(fieldWidth) << j << "|";
        for (int i = grid_->u().beginI(); i < grid_->u().endI(); i++) {
            file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << grid_->f(i, j);
        }
        file << std::endl;
    }
    file << std::endl;

    // write g
    // ---------
    // write header lines
    file << "G (" << grid_->v().size()[0] << "x" << grid_->v().size()[1] << "): " << std::endl << std::string(fieldWidth, ' ') << "|";
    for (int i = grid_->v().beginI(); i < grid_->v().endI(); i++) {
        file << std::setw(fieldWidth) << i;
    }
    file << std::endl << std::string(fieldWidth * (grid_->v().size()[0] + 2) + 1, '-') << std::endl;

    // write g values
    for (int j = grid_->v().endJ() - 1; j >= grid_->v().beginJ(); j--) {
        file << std::setw(fieldWidth) << j << "|";
        for (int i = grid_->v().beginI(); i < grid_->v().endI(); i++) {
            file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << grid_->g(i, j);
        }
        file << std::endl;
    }
    file << std::endl;

    // write rhs
    // ---------
    // write header lines
    file << "rhs (" << grid_->p().size()[0] << "x" << grid_->p().size()[1] << "): " << std::endl << std::string(fieldWidth, ' ') << "|";
    for (int i = grid_->p().beginI(); i < grid_->p().endI(); i++) {
        file << std::setw(fieldWidth) << i;
    }
    file << std::endl << std::string(fieldWidth * (grid_->p().size()[0] + 2) + 1, '-') << std::endl;

    // write rhs values
    for (int j = grid_->p().endJ() - 1; j >= grid_->p().beginJ(); j--) {
        file << std::setw(fieldWidth) << j << "|";
        for (int i = grid_->p().beginI(); i < grid_->p().endI(); i++) {
            file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << grid_->rhs(i, j);
        }
        file << std::endl;
    }
    file << std::endl;

    // write forces
    // ---------
    // top forces
    // write header lines
    file << "Forces (" << grid_->topF().size()[0] << "x" << grid_->topF().size()[1] << "): " << std::endl << std::string(fieldWidth, ' ') << "|";
    for (int i = grid_->topF().beginI(); i < grid_->topF().endI(); i++) {
        file << std::setw(fieldWidth) << i;
    }
    file << std::endl << std::string(fieldWidth * (grid_->topF().size()[0] + 2) + 1, '-') << std::endl;

    // write force values together in one field
    for (int j = grid_->p().endJ() - 1; j >= grid_->p().beginJ(); j--) {
        file << std::setw(fieldWidth) << j << "|";
        if (j == grid_->p().endJ() - 1) {
            for (int i = grid_->topF().beginI(); i < grid_->topF().endI(); i++) {
                file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << grid_->topF(i);
            }
            file << std::endl;
        } else if (j == grid_->p().beginJ()) {
            for (int i = grid_->bottomF().beginI(); i < grid_->bottomF().endI(); i++) {
                file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << grid_->bottomF(i);
            }
            file << std::endl;
        } else {
            for (int i = grid_->p().beginI(); i < grid_->p().endI(); i++) {
                file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << 0;
            }
            file << std::endl;
        }
    }
    file << std::endl;
}

void OutputWriterTextParallel::writePressureFile() const {
    // counter for files, counter value is part of the file name
    static int pressurefileNo = 0;

    // Assemble the filename
    std::stringstream fileName;
    fileName << "out/pressure_" << std::setw(4) << std::setfill('0') << pressurefileNo++ << "." << partitioning_.ownRank() << ".txt";

    // open file
    std::ofstream file(fileName.str().c_str());

    if (!file.is_open()) {
        std::cout << "Could not write to file \"" << fileName.str() << "\".";
        return;
    }

    // write mesh width
    file << "nCells: " << grid_->nCells()[0] << "x" << grid_->nCells()[1] << ", dx: " << grid_->dx() << ", dy: " << grid_->dy() << std::endl
         << std::endl;

    const int fieldWidth = 9; // number of characters to use for a single value

    // write p
    // ---------
    // write header lines
    file << "p (" << grid_->p().size()[0] << "x" << grid_->p().size()[1] << "): " << std::endl << std::string(fieldWidth, ' ') << "|";
    for (int i = grid_->p().beginI(); i < grid_->p().endI(); i++) {
        file << std::setw(fieldWidth) << i;
    }
    file << std::endl << std::string(fieldWidth * (grid_->p().size()[0] + 2) + 1, '-') << std::endl;

    // write p values
    for (int j = grid_->p().endJ() - 1; j >= grid_->p().beginJ(); j--) {
        file << std::setw(fieldWidth) << j << "|";
        for (int i = grid_->p().beginI(); i < grid_->p().endI(); i++) {
            file << std::setw(fieldWidth) << std::setprecision(fieldWidth - 6) << grid_->p(i, j);
        }
        file << std::endl;
    }
    file << std::endl;
}