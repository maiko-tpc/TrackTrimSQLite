#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <functional>

#include "TRIM2SQLite.hpp"
#include "CollisionDBHandler.hpp"
#include "VectorAndMatrix.hpp"

class TrackGenerator
{
public:
    using Vector = TRIM2SQLite::CollisionRecord::xyz;
    using Transform = Mm::Matrix<Vector>;
    using Collision = TRIM2SQLite::CollisionRecord;
    using CollisionCollection = std::vector<TRIM2SQLite::CollisionRecord>;
    using CollisionIterator = CollisionCollection::iterator;

    TrackGenerator()
        : fFileName(), fRandomGenerator([]() { return 0; }),
          fTrackIDMin(-1), fTrackIDMax(-1), fAccessibilityGood(false){};

    TrackGenerator(const std::string &_fFileName)
        : fFileName(), fRandomGenerator([]() { return 0; }),
          fTrackIDMin(-1), fTrackIDMax(-1), fAccessibilityGood(false)
    {
        auto ok = SetFileName(_fFileName);
        if (!ok)
        {
            std::cerr << "TrackGenerator() :: File not found." << std::endl;
        }
    };

    virtual CollisionCollection
    Generate(double _ekin, double _x, double _y, double _z,
             double _dx, double _dy, double _dz) = 0;

    bool SetFileName(const std::string &_fFileName)
    {
        fFileName = _fFileName;
        return CheckAccessibility();
    };

    std::string GetFileName() const { return fFileName; };
    void SetRandomGenerator(const std::function<double()> &_fRandomGenerator)
    {
        fRandomGenerator = _fRandomGenerator;
    };

    void SetTrackIDMin(int _fTrackIDMin)
    {
        fTrackIDMin = _fTrackIDMin;
    };

    void SetTrackIDMax(int _fTrackIDMax)
    {
        fTrackIDMax = _fTrackIDMax;
    };

    double GetTrackIDMin() const
    {
        return fTrackIDMin;
    };

    double GetTrackIDMax() const
    {
        return fTrackIDMax;
    };

    CollisionCollection GetTrackRandom()
    {
        CollisionDBHandler db(GetFileName());
        int nTracks = db.GetNumberOfTracks();
        int iTrack;
        int iMin = fTrackIDMin >= 0 ? fTrackIDMin : 0;
        int iMax = (0 <= fTrackIDMax && fTrackIDMax < nTracks) ? fTrackIDMax : nTracks - 1;

        if (iMin <= iMax)
            iTrack = GetRandomInteger(iMin, iMax);
        else
        {
            std::cerr << "TrackGenerator::GetTracsRandom() :: Invalid index range -> use all tracks." << std::endl;
            iTrack = GetRandomInteger(0, nTracks - 1);
        }

        return db.GetTrack(iTrack);
    };

    CollisionCollection GetTrack(int _trackID)
    {
        CollisionDBHandler db(GetFileName());

        return db.GetTrack(_trackID);
    };

    bool IsAccesible() const
    {
        return fAccessibilityGood;
    };

    bool CheckAccessibility()
    {
        try
        {
            CollisionDBHandler db(GetFileName());
            fAccessibilityGood = true;
        }
        catch (...)
        {
            fAccessibilityGood = false;
        }

        return fAccessibilityGood;
    }

private:
    std::string fFileName;
    std::function<double()> fRandomGenerator;
    int fTrackIDMin, fTrackIDMax;
    bool fAccessibilityGood;

protected:
    double Random()
    {
        double val = fRandomGenerator();
        if (val < 0)
            return 0;
        else if (val > 1)
            return 1;
        else
            return val;
    };

    int GetRandomInteger(int _min, int _max)
    {
        int min = _min <= _max ? _min : _max;
        int max = _min <= _max ? _max : _min;

        int a = min + Random() * (max - min + 1);
        return a < max ? a : max;
    }
};

class TrackGeneratorA : public TrackGenerator
{
public:
    TrackGeneratorA()
        : TrackGenerator(){};

    TrackGeneratorA(const std::string &_fFileName)
        : TrackGenerator(_fFileName){};

    virtual CollisionCollection
    Generate(double _ekin, double _x, double _y, double _z,
             double _dx, double _dy, double _dz)
    {
        if (!IsAccesible())
        {
            throw std::runtime_error("Generate() :: DB " + GetFileName() + " is not accessible...");
        }

        if (_dx == 0 && _dy == 0 && _dz == 0)
        {
            throw std::runtime_error("Generate() :: Zero vector input");
        }

        auto track = GetTrackRandom();
        EraseHighEnergyCollisionsFromTrack(track, _ekin);

        Transform mat;
        for (auto it = track.begin(); it != track.end(); ++it)
        {
            if (it == track.begin())
            {
                mat = SetFirstCollision(it, _ekin, _x, _y, _z,
                                        _dx, _dy, _dz);
            }

            else
            {
                AdaptCollision(it, mat);
            }
        }

        return track;
    };

protected:
    void EraseHighEnergyCollisionsFromTrack(CollisionCollection &_track,
                                            double _ekin);

    Transform SetFirstCollision(CollisionIterator _it,
                                double _ekin, double _x, double _y, double _z,
                                double _dx, double _dy, double _dz) const;

    // Adapt collision at _it to the collision just before
    // _it must not be the first iterator
    void AdaptCollision(CollisionIterator _it, Transform _mat);
};

// phi-rotation at each collisions
class TrackGeneratorB : public TrackGeneratorA
{
public:
    TrackGeneratorB()
        : TrackGeneratorA(){};

    TrackGeneratorB(const std::string &_fFileName)
        : TrackGeneratorA(_fFileName){};

    virtual CollisionCollection
    Generate(double _ekin, double _x, double _y, double _z,
             double _dx, double _dy, double _dz)
    {

        if (!IsAccesible())
        {
            throw std::runtime_error("Generate() :: DB " + GetFileName() + " is not accessible...");
        }

        if (_dx == 0 && _dy == 0 && _dz == 0)
        {
            throw std::runtime_error("Generate()::Zero vector input");
        }

        auto track = GetTrackRandom();
        EraseHighEnergyCollisionsFromTrack(track, _ekin);

        Transform mat;
        for (auto it = track.begin(); it != track.end(); ++it)
        {
            if (it == track.begin())
            {
                mat = SetFirstCollision(it, _ekin, _x, _y, _z,
                                        _dx, _dy, _dz);

                PhiRotationRandom(it, mat);
            }

            else
            {
                AdaptCollision(it, mat);
                PhiRotationRandom(it, mat);
            }
        }

        return track;
    };

protected:
    void PhiRotationRandom(CollisionIterator _it, Transform &_mat);
};

class TrackGeneratorC : public TrackGeneratorB
{
public:
    TrackGeneratorC()
        : TrackGeneratorB(), fTransferProbability(0), fEnergyMarginRatio(0.5){};

    TrackGeneratorC(const std::string &_fFileName)
        : TrackGeneratorB(_fFileName), fTransferProbability(0), fEnergyMarginRatio(0.5){};

    virtual CollisionCollection
    Generate(double _ekin, double _x, double _y, double _z,
             double _dx, double _dy, double _dz)
    {

        if (!IsAccesible())
        {
            throw std::runtime_error("Generate() :: DB " + GetFileName() + " is not accessible...");
        }

        if (_dx == 0 && _dy == 0 && _dz == 0)
        {
            throw std::runtime_error("Generate()::Zero vector input");
        }

        auto track = GetTrackRandom();
        EraseHighEnergyCollisionsFromTrack(track, _ekin);

        Transform mat;
        for (auto it = track.begin(); it != track.end(); ++it)
        {
            if (it == track.begin())
            {
                mat = SetFirstCollision(it, _ekin, _x, _y, _z,
                                        _dx, _dy, _dz);

                PhiRotationRandom(it, mat);
            }

            else
            {
                AdaptCollision(it, mat);
                PhiRotationRandom(it, mat);

                // Kinetic energy after current collision
                double ene = it->GetIncidentEnergy() - it->GetRecoilEnergy();
                // Energy loss before next collision in current collision sequence
                double de = it->GetEnergyLoss();
                // Distance to next collision in current collision sequence
                double dr = it->GetDistanceToNextCollision();
                // Direction to next collision
                auto dirScattering = it->GetScatteringDirection();

                //Transfer?
                if (de > 0 &&
                    dr > 0 &&
                    Mm::Norm(dirScattering) > 0 &&
                    ene - de > 0 &&
                    Random() < GetTransferProbability())
                {
                    Transfer(track, it, mat);
                }
            }
        }

        return track;
    };

    double GetTransferProbability() const
    {
        return fTransferProbability;
    };

    void SetTransferProbability(double _fTransferProbability)
    {
        if (0 <= _fTransferProbability &&
            _fTransferProbability <= 1)
            fTransferProbability = _fTransferProbability;
    };

    double GetEnergyMarginRatio() const { return fEnergyMarginRatio; };
    void SetEnergyMarginRatio(double _fEnergyMarginRatio)
    {
        if (_fEnergyMarginRatio < 0)
            fEnergyMarginRatio = 0;
        else if (1 < _fEnergyMarginRatio)
            fEnergyMarginRatio = 1;
        else
            fEnergyMarginRatio = _fEnergyMarginRatio;
    };

protected:
    CollisionCollection GetTransferDestinationCandidates(double _ene,
                                                         double _ene_min, double _ene_max);

    bool Transfer(CollisionCollection &_track, CollisionIterator &_it, Transform &_mat);

private:
    double fTransferProbability;
    double fEnergyMarginRatio;
};
