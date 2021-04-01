#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

#include <sqlite3.h>

class RANGE_3D
{
public:
    RANGE_3D(const std::string &_filename);

    bool Next();

    bool Good() const { return fGood; };
    std::string GetIonName() const { return fIonName; };
    int GetAtomicNumber() const { return fAtomicNumber; };
    double GetIonMassAMU() const { return fIonMassAMU; };
    int GetMassNumber() const { return std::round(fIonMassAMU); };
    double GetIncidentEnergy() const { return fIonEnergy; };
    double GetIncidentAngle() const { return fIonIncidentAngle; };
    int GetIonNumber() const { return fIonNumber; };
    double GetDepth() const { return fDepthX; };
    double GetLateralY() const { return fLateralY; };
    double GetLateralZ() const { return fLateralZ; };

private:
    std::ifstream fIfs;
    bool fGood;
    std::string fIonName;
    int fAtomicNumber;
    double fIonMassAMU;
    double fIonEnergy; // eV
    double fIonIncidentAngle;
    int fIonNumber;
    double fDepthX, fLateralY, fLateralZ; // cm
};

//COLLISON.txt
//     ~~~ <-- !!!
class COLLISON
{
public:
    COLLISON(const std::string &_filename);

    bool Next();

    bool Good() const { return fGood; };
    std::string GetIonName() const { return fIonName; };
    double GetIonMassAMU() const { return fIonMassAMU; };
    int GetMassNumber() const { return std::round(fIonMassAMU); };
    double GetIonEnergy() const { return fIonEnergy; };

    int GetIonNumber() const { return fIonNumber; };
    const std::vector<double> &GetEnergy() const { return fvEnergy; };
    const std::vector<double> &GetDepth() const { return fvDepth; };
    const std::vector<double> &GetLateralY() const { return fvLateralY; };
    const std::vector<double> &GetLateralZ() const { return fvLateralZ; };
    const std::vector<double> &GetSe() const { return fvSe; };
    const std::vector<std::string> &GetAtomHit() const { return fvAtomHit; };
    const std::vector<double> &GetRecoilEnergy() const { return fvRecoilEnergy; };

private:
    std::ifstream fIfs;
    bool fGood;
    std::string fIonName;
    double fIonMassAMU;
    double fIonEnergy; // keV
    int fIonNumber;
    std::vector<double> fvEnergy, fvDepth, fvLateralY, fvLateralZ, fvSe;
    std::vector<std::string> fvAtomHit;
    std::vector<double> fvRecoilEnergy;

    void Clear()
    {
        fvEnergy.clear();
        fvDepth.clear();
        fvLateralY.clear();
        fvLateralZ.clear();
        fvSe.clear();
        fvAtomHit.clear();
        fvRecoilEnergy.clear();
    }
};

class TRIM2SQLite
{
public:
    TRIM2SQLite(){};

    bool MakeSQLiteFile(const std::string &_path, const std::string &_outputname);

    class CollisionRecord
    {
    public:
        CollisionRecord(){};

        struct xyz
        {
            xyz() : fX(0), fY(0), fZ(0){};
            xyz(double _fX, double _fY, double _fZ)
                : fX(_fX), fY(_fY), fZ(_fZ){};

            double X() const { return fX; };
            double Y() const { return fY; };
            double Z() const { return fZ; };
            double fX, fY, fZ;
        };

        void SetTrackID(int _fTrackID) { fTrackID = _fTrackID; }
        void SetCollisionID(int _fCollisionID) { fCollisionID = _fCollisionID; }
        void SetIncidentEnergy(double _fIncidentEnergy) { fIncidentEnergy = _fIncidentEnergy; }
        void SetIncidentIon(const std::string &_fIncidentIon) { fIncidentIon = _fIncidentIon; }
        void SetMassNumber(int _fMassNumber) { fMassNumber = _fMassNumber; };
        void SetRecoilIon(const std::string &_fRecoilIon) { fRecoilIon = _fRecoilIon; }
        void SetRecoilEnergy(double _fRecoilEnergy) { fRecoilEnergy = _fRecoilEnergy; }
        void SetPosition(double _x, double _y, double _z)
        {
            fPosition = xyz(_x, _y, _z);
        }

        void SetIncidentDirection(double _dx, double _dy, double _dz)
        {
            double norm = std::sqrt(_dx * _dx + _dy * _dy + _dz * _dz);
            if (norm > 0)
                fIncidentDirection = xyz(_dx / norm, _dy / norm, _dz / norm);
            else
                fIncidentDirection = xyz(0, 0, 0);
        }
        void SetScatteringDirection(double _dx, double _dy, double _dz)
        {
            double norm = std::sqrt(_dx * _dx + _dy * _dy + _dz * _dz);
            if (norm > 0)
                fScatteringDirection = xyz(_dx / norm, _dy / norm, _dz / norm);
            else
                fScatteringDirection = xyz(0, 0, 0);
        }

        void SetPosition(xyz &_vec)
        {
            SetPosition(_vec.X(), _vec.Y(), _vec.Z());
        }

        void SetIncidentDirection(xyz &_vec)
        {
            SetIncidentDirection(_vec.X(), _vec.Y(), _vec.Z());
        }

        void SetScatteringDirection(xyz &_vec)
        {
            SetScatteringDirection(_vec.X(), _vec.Y(), _vec.Z());
        }

        void SetDistanceToNextCollision(double _fDistanceToNextCollision)
        {
            fDistanceToNextCollision = _fDistanceToNextCollision;
        }
        void SetEnergyLoss(double _fEnergyLoss)
        {
            fEnergyLoss = _fEnergyLoss;
        }

        int GetTrackID() const { return fTrackID; }
        int GetCollisionID() const { return fCollisionID; }

        double GetIncidentEnergy() const { return fIncidentEnergy; }
        const std::string &GetIncidentIon() const { return fIncidentIon; }
        int GetMassNumber() const { return fMassNumber; };
        const std::string &GetRecoilIon() const { return fRecoilIon; }
        double GetRecoilEnergy() const { return fRecoilEnergy; }
        const xyz &GetPosition() const
        {
            return fPosition;
        }

        const xyz &GetIncidentDirection() const
        {
            return fIncidentDirection;
        }
        const xyz &GetScatteringDirection() const
        {
            return fScatteringDirection;
        }
        double GetDistanceToNextCollision() const
        {
            return fDistanceToNextCollision;
        }
        double GetEnergyLoss() const
        {
            return fEnergyLoss;
        }

    private:
        int fTrackID;
        int fCollisionID;
        double fIncidentEnergy;
        std::string fIncidentIon;
        int fMassNumber;
        std::string fRecoilIon;
        double fRecoilEnergy;
        xyz fPosition;
        xyz fIncidentDirection;
        xyz fScatteringDirection;
        double fDistanceToNextCollision;
        double fEnergyLoss; // De due to electrons, except for one from collision w/ atom
    };

    static const std::string NameOfTable;

private:
};
