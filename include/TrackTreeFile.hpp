#pragma once

#include <iostream>
#include <vector>
#include <string>

#include <TTree.h>
#include <TFile.h>

#include "TRIM2SQLite.hpp"

class TrackTreeFile
{
public:
    TrackTreeFile(const std::string &_fileName);

    ~TrackTreeFile();

    bool IsOpen() const;

    void Fill(std::vector<TRIM2SQLite::CollisionRecord> _track);

    bool IsWritten() const { return fWritten; };

    void Write();

    void Clear();

private:
    TFile fFile;
    TTree fTree;

    bool fWritten;

    std::string ion;
    int mass_number;
    double ene0;
    std::vector<int> track_id, collision_id;
    std::vector<double> vX, vY, vZ, vEne;
    // difference between next step (X == depth)
    std::vector<double> vdR;
    std::vector<double> vdX0, vdY0, vdZ0;
    std::vector<double> vdX1, vdY1, vdZ1;
    std::vector<double> vdEne;
    std::vector<double> vTh;        // THeta == scattering angle
    std::vector<std::string> vAtom; //Atom hit
    std::vector<double> vRecoil;    // Recoil energy
};
