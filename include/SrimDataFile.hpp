#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstring>

namespace Mm
{
class SrimDataFile
{
public:
    SrimDataFile();
    ~SrimDataFile();
    bool ReadFile(const std::string &_filename);

    void Show()const;
    const std::vector<double>& GetEkin()const{return fvEkin;}
    const std::vector<double>& GetEmLoss()const{return fvEmLoss;}
    const std::vector<double>& GetHdLoss()const{return fvHdLoss;}
    const std::vector<double>& GetRange()const{return fvRange;}
    const std::vector<double>& GetTransverseStraggling()const{return fvTransverseStraggling;}
    const std::vector<double>& GetLongitudinalStraggling()const{return fvLongitudinalStraggling;}

private:
    double fZ;
    double fMass;
    double fDensity;

    std::vector<double> fvEkin;
    std::vector<double> fvEmLoss;
    std::vector<double> fvHdLoss;
    std::vector<double> fvRange;
    std::vector<double> fvTransverseStraggling;
    std::vector<double> fvLongitudinalStraggling;

    bool fDebug;
};


} // namespace Mm