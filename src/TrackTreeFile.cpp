#include "TrackTreeFile.hpp"

namespace
{

    template <typename T>
    std::vector<T> MakeScatteringAngleVector(const std::vector<T> &_vdx,
                                             const std::vector<T> &_vdy,
                                             const std::vector<T> &_vdz)
    {

        std::vector<T> ret;
        const int vec_size = _vdx.size();
        if (vec_size != _vdy.size() || vec_size != _vdz.size())
        {
            throw std::runtime_error("MakeScatteringAngleVector :: size mismatch");
        }
        ret.reserve(vec_size);

        for (int i = 0; i < vec_size; ++i)
        {
            if (i == 0)
            {
                ret.push_back(0); // Dont know before injection
            }
            else // last
            {

                double angle = -1000;

                double x_f = _vdx.at(i - 1);
                double x_l = _vdx.at(i);
                double y_f = _vdy.at(i - 1);
                double y_l = _vdy.at(i);
                double z_f = _vdz.at(i - 1);
                double z_l = _vdz.at(i);

                double inner_prd = x_f * x_l + y_f * y_l + z_f * z_l;
                double amp_f = sqrt(x_f * x_f + y_f * y_f + z_f * z_f);
                double amp_l = sqrt(x_l * x_l + y_l * y_l + z_l * z_l);

                if (amp_f > 0 && amp_l > 0)
                    angle = acos(inner_prd / amp_f / amp_l) / M_PI * 180;

                ret.push_back(angle);
            }
        }
        return ret;
    }
}

TrackTreeFile::TrackTreeFile(const std::string &_fileName)
    : fFile(_fileName.c_str(), "recreate"), fTree("tr", "tracks"), fWritten(false)
{

    fTree.Branch("ion", &ion);
    fTree.Branch("mass_number", &mass_number, "mass_number/I");
    fTree.Branch("ene0", &ene0, "ene0/D");
    fTree.Branch("track_id", &track_id);
    fTree.Branch("collision_id", &collision_id);
    fTree.Branch("x", &vX);
    fTree.Branch("y", &vY);
    fTree.Branch("z", &vZ);
    fTree.Branch("dr", &vdR);
    fTree.Branch("ene", &vEne);
    fTree.Branch("dx0", &vdX0);
    fTree.Branch("dy0", &vdY0);
    fTree.Branch("dz0", &vdZ0);
    fTree.Branch("dx1", &vdX1);
    fTree.Branch("dy1", &vdY1);
    fTree.Branch("dz1", &vdZ1);
    fTree.Branch("dene", &vdEne);
    fTree.Branch("th", &vTh);
    fTree.Branch("recoil", &vRecoil);
    fTree.Branch("atom", &vAtom);
};

TrackTreeFile::~TrackTreeFile()
{
    if (!IsWritten())
    {
        Write();
    }
};

bool TrackTreeFile::IsOpen() const { return fFile.IsOpen(); };

void TrackTreeFile::Fill(std::vector<TRIM2SQLite::CollisionRecord> _track)
{
    Clear();

    if (_track.size() != 0)
    {
        ion = _track.front().GetIncidentIon();
        mass_number = _track.front().GetMassNumber();
        ene0 = _track.front().GetIncidentEnergy();
    }

    for (auto &col : _track)
    {

        track_id.push_back(col.GetTrackID());
        collision_id.push_back(col.GetCollisionID());

        auto pos = col.GetPosition();
        vX.push_back(pos.X());
        vY.push_back(pos.Y());
        vZ.push_back(pos.Z());
        vEne.push_back(col.GetIncidentEnergy());

        vdR.push_back(col.GetDistanceToNextCollision());

        auto dx0 = col.GetIncidentDirection();
        vdX0.push_back(dx0.X());
        vdY0.push_back(dx0.Y());
        vdZ0.push_back(dx0.Z());

        auto dx1 = col.GetScatteringDirection();
        vdX1.push_back(dx1.X());
        vdY1.push_back(dx1.Y());
        vdZ1.push_back(dx1.Z());

        vdEne.push_back(col.GetEnergyLoss());
        vRecoil.push_back(col.GetRecoilEnergy());
        vAtom.push_back(col.GetRecoilIon());
    }

    vTh = MakeScatteringAngleVector(vdX0, vdY0, vdZ0);
    fTree.Fill();
};

void TrackTreeFile::Write()
{
    if (IsOpen())
    {
        fTree.Write();
        fWritten = true;
    }
};

void TrackTreeFile::Clear()
{

    ion = "";
    mass_number = 0;
    ene0 = 0;
    track_id.clear();
    collision_id.clear();

    vX.clear();
    vY.clear();
    vZ.clear();
    vEne.clear();
    vdR.clear();
    vdX0.clear();
    vdY0.clear();
    vdZ0.clear();
    vdX1.clear();
    vdY1.clear();
    vdZ1.clear();
    vTh.clear();
    vdEne.clear();
    vAtom.clear();
    vRecoil.clear();
}