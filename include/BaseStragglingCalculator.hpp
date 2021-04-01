#pragma once
#include <iostream>

namespace Mm
{
class BaseStragglingCalculator
{
public:
    BaseStragglingCalculator();
    ~BaseStragglingCalculator();
    virtual double GetTransverseStraggling() const = 0;
    virtual double GetLongitudinalStraggling() const = 0;

    inline void SetZ(int _fZ) { fZ = _fZ; };
    inline void SetA(int _fA) { fA = _fA; };
    inline void SetMass(double _fMass) { fMass = _fMass; };

    inline int GetZ() const { return fZ; };
    inline int GetA() const { return fA; };
    inline double GetMass() const { return fMass; };

private:
    int fZ;
    int fA;
    double fMass;
};
} // namespace Mm