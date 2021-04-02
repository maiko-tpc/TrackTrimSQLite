#include "TRIM2SQLite.hpp"

namespace
{

    bool SanitizeEndOfLine(std::string &_str)
    {
        if (_str.length() > 0 && _str.back() == '\r')
        {
            _str.resize(_str.size() - 1);
            return true;
        }
        else
        {
            return false;
        }
    }

    void EraseBlank(std::string &_str)
    {
        for (auto it = _str.begin();;)
        {
            if (*it == ' ')
            {
                it = _str.erase(it);
            }
            else if (it == _str.end())
            {
                break;
            }
            else
            {
                ++it;
            }
        }
    }

    template <typename T>
    std::vector<T> CalculateDifference(const std::vector<T> &_in)
    {
        std::vector<T> ret;
        ret.reserve(_in.size());
        for (auto it = _in.begin(); it != _in.end(); ++it)
        {
            if (it + 1 != _in.end())
            {
                ret.push_back(*(it + 1) - *(it));
            }
            else // last
            {
                ret.push_back(0);
            }
        }
        return ret;
    }

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

RANGE_3D::RANGE_3D(const std::string &_filename)
    : fIfs(_filename), fGood(false), fIonName(), fAtomicNumber(0), fIonMassAMU(0),
      fIonEnergy(0), fIonNumber(-1), fDepthX(0), fLateralY(0), fLateralZ(0)
{
    if (fIfs.good())
    {
        std::string sLine;
        //Header
        while (std::getline(fIfs, sLine))
        {
            SanitizeEndOfLine(sLine);

            if (sLine.substr(0, 5) == "Ion =" && sLine.substr(16, 9) == "Ion Mass=")
            //if(sLine.substr(0, 5) == "Ion =")
            {
                auto name = sLine.substr(5, 3);
                EraseBlank(name);
                fIonName = name;
                auto atomic_number = sLine.substr(10, 2);

                auto ion_mass = sLine.substr(26, 8);
                fAtomicNumber = std::stoi(atomic_number);
                fIonMassAMU = std::stoi(ion_mass);
            }
            else if (sLine.substr(0, 9) == "Energy  =")
            {
                auto tmp = sLine.substr(9, 13);
                EraseBlank(tmp);
                auto unit = sLine.substr(23, 3);

                if (unit == "keV")
                {
                    fIonEnergy = std::stod(tmp) * 1e+3;
                }
                else
                {
                    std::cerr << "RANGE_3D(): Unknown unit for ion energy. -> " << unit << std::endl;
                    fGood = false;
                    break;
                }
            }
            else if (sLine.substr(0, 22) == "Ion Angle to Surface =")
            {

                auto tmp = sLine.substr(22, 5);
                EraseBlank(tmp);
                auto unit = sLine.substr(28, 7);
                if (unit == "degrees")
                    fIonIncidentAngle = std::stod(tmp);
                else
                {
                    std::cerr << "RANGE_3D(): Unknown unit for incident angle. -> " << unit << std::endl;
                    fIonIncidentAngle = 0;
                    fGood = false;
                }
            }
            else if (sLine.substr(0, 7) == "-------")
            {
                fGood = true;
                break;
            }
        }
    }
};

bool RANGE_3D::Next()
{
    if (!fGood)
        return false;

    std::string sLine;
    if (std::getline(fIfs, sLine))
    {
        SanitizeEndOfLine(sLine);
        std::istringstream ssLine(sLine);
        if (ssLine >> fIonNumber >> fDepthX >> fLateralY >> fLateralZ)
        {
            fDepthX /= 1.e8;   //cm
            fLateralY /= 1.e8; //cm
            fLateralZ /= 1.e8; //cm
            fGood = true;
        }
        else
        {
            fGood = false;
        }
    }
    else
    {
        fGood = false;
    }
    return fGood;
}

COLLISON::COLLISON(const std::string &_filename)
    : fIfs(_filename), fGood(false), fIonName(), fIonMassAMU(0), fIonEnergy(0), fIonNumber(-1)
{
    if (fIfs.good())
    {
        std::string sLine;
        //Header
        while (std::getline(fIfs, sLine))
        {

            SanitizeEndOfLine(sLine);
            if (sLine.substr(0, 18) == "o     Ion Name   =")
            {
                auto tmp = sLine.substr(19, 11);
                EraseBlank(tmp);
                fIonName = tmp;
            }
            else if (sLine.substr(0, 18) == "o     Ion Mass   =")
            {
                auto tmp = sLine.substr(18, 11);

                auto amu = sLine.substr(30, 3);
                if (amu == "amu")
                    fIonMassAMU = std::stod(tmp);
                else
                {
                    std::cerr << "COLLISON(): Unknown unit for ion mass. -> " << amu << std::endl;
                    fIonMassAMU = 0;
                    fGood = false;
                    break;
                }
            }

            else if (sLine.substr(0, 18) == "o     Ion Energy =")
            {
                auto tmp = sLine.substr(18, 11);
                auto unit = sLine.substr(30, 3);
                if (unit == "keV")
                {
                    fIonEnergy = std::stod(tmp) * 1.0e+3; // eV
                }
                else
                {
                    std::cerr << "COLLISON(): Unknown unit for ion energy. -> " << unit << std::endl;
                    fGood = false;
                    break;
                }
            }
            else if (sLine == std::string(102, '-'))
            {
                fGood = true;
                break;
            }
        }
    }
};

bool COLLISON::Next()
{
    if (!fGood)
        return false;

    Clear();
    std::string sLine;
    fGood = false;
    fIonNumber = -1;
    while (std::getline(fIfs, sLine))
    {
        SanitizeEndOfLine(sLine);
        if (sLine == std::string(102, '='))
        {
            fGood = true;
            break;
        }
        int ionNumber = std::stoi(sLine.substr(1, 5));
        if (fIonNumber > 0 && fIonNumber != ionNumber)
        {
            std::cout << " ion number mismatch -> skipping this line" << std::endl;
            continue;
        }
        else
        {
            fIonNumber = ionNumber;
        }

        fvEnergy.push_back(std::stod(sLine.substr(7, 9)) * 1.0e+3); //eV

        fvDepth.push_back(std::stod(sLine.substr(17, 10)) / 1.e8); //cm

        fvLateralY.push_back(std::stod(sLine.substr(28, 10)) / 1.e8); //cm

        fvLateralZ.push_back(std::stod(sLine.substr(39, 10)) / 1.e8); //cm

        fvSe.push_back(std::stod(sLine.substr(50, 7)) * 1.0e+3); //eV

        {
            std::string atom_tmp = sLine.substr(58, 4);
            //std::cout << "     /" << sLine.substr(58, 4) << "/" << std::endl;
            EraseBlank(atom_tmp);
            fvAtomHit.push_back(atom_tmp);
        }
        fvRecoilEnergy.push_back(std::stod(sLine.substr(63, 10))); //eV
    }

    //ender
    if (fGood)
    {
        fGood = false;
        while (std::getline(fIfs, sLine))
        {
            SanitizeEndOfLine(sLine);
            if (sLine == std::string(102, '-'))
            {
                fGood = true;
                break;
            }
        }
    }

    return fGood;
}

const std::string TRIM2SQLite::NameOfTable = "collisions";

bool TRIM2SQLite::MakeSQLiteFile(const std::string &_path,
                                 const std::string &_outputname)
{

    COLLISON col(_path + "/COLLISON.txt");
    RANGE_3D rng(_path + "/RANGE_3D.txt");

    if (!col.Good() || !rng.Good())
    {
        return false;
    }

    sqlite3 *pDB;
    auto err = sqlite3_open_v2(_outputname.c_str(), &pDB, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, nullptr);
    if (err != SQLITE_OK)
    {
        sqlite3_close(pDB);
        pDB = nullptr;
        throw std::runtime_error("TRIM2SQLite::MakeSQLiteFile() DB file can not be opend.");
    }

    char *errMsg = nullptr;
    std::ostringstream qCreateTable;
    qCreateTable << "CREATE TABLE IF NOT EXISTS " << NameOfTable << " "
                 << "("
                 << "track_id INTEGER, "
                 << "collision_id INTEGER, "
                 << "e_inc REAL,  "
                 << "incident_ion TEXT,  "
                 << "incident_ion_mass INTEGER,  "
                 << "recoil_ion TEXT,  "
                 << "e_rec REAL,  "
                 << "x REAL, "
                 << "y REAL, "
                 << "z REAL, "
                 << "dx0 REAL, "
                 << "dy0 REAL, "
                 << "dz0 REAL, "
                 << "dx1 REAL, "
                 << "dy1 REAL, "
                 << "dz1 REAL, "
                 << "dr REAL, "
                 << "de REAL,"
                 << "PRIMARY KEY (track_id, collision_id)"
                 << ");";

    err = sqlite3_exec(pDB, qCreateTable.str().c_str(),
                       NULL, NULL, &errMsg);

    std::ostringstream qInsert;
    qInsert << "INSERT INTO " << NameOfTable << " "
            << "("
            << "track_id, "
            << "collision_id, "
            << "e_inc,  "
            << "incident_ion,  "
            << "incident_ion_mass,  "
            << "recoil_ion,  "
            << "e_rec,  "
            << "x, "
            << "y, "
            << "z, "
            << "dx0, "
            << "dy0, "
            << "dz0, "
            << "dx1, "
            << "dy1, "
            << "dz1, "
            << "dr, "
            << "de"
            << ") "
            << "VALUES "
            << "("
            << ":TID, "
            << ":CID, "
            << ":EIN,  "
            << ":IION,  "
            << ":AION,  "
            << ":RION,  "
            << ":ERC,  "
            << ":X, "
            << ":Y, "
            << ":Z, "
            << ":DX0, "
            << ":DY0, "
            << ":DZ0, "
            << ":DX1, "
            << ":DY1, "
            << ":DZ1, "
            << ":DR, "
            << ":DE); ";

    sqlite3_stmt *query;
    err = sqlite3_prepare_v2(pDB,
                             qInsert.str().c_str(),
                             -1, &query, nullptr);
    if (err != SQLITE_OK)
    {
        sqlite3_close(pDB);
        sqlite3_finalize(query);
        pDB = nullptr;
        throw std::runtime_error("logic error");
    }

    int trackID = 0;
    while (col.Next() && rng.Next())
    {

        // Ion number mismatch occurs if a track do NOT stop in the medium.
        // -> Skip unstopped tracks
        if (col.GetIonNumber() != rng.GetIonNumber())
        {

            int nCol = col.GetIonNumber();
            int nRng = rng.GetIonNumber();

            std::cerr << "TRIM2SQLite::MakeSQLiteFile() :  Track #" << nCol
                      << " is not stopped in the medium." << std::endl;

            // Skip collisions till nCol == nRng
            while (true)
            {
                if (!col.Next() || col.GetIonNumber() > nRng)
                {
                    std::cerr << "                             Unexpected input !" << std::endl;
                    std::cerr << "                             Aborted. " << std::endl;
                    sqlite3_finalize(query);
                    sqlite3_close(pDB);
                    return false;
                }
                if (col.GetIonNumber() == nRng)
                    break;
            }

            std::cerr << "                             track #" << nCol << " -- " << nRng
                      << " is ignored and track# in DB is reassigned. " << std::endl;
        }

        std::string ion;
        int mass_number;
        double ene0;

        std::vector<double> vX, vY, vZ, vEne;
        // difference between next step (X == depth)
        std::vector<double> vdX, vdY, vdZ, vdEne;
        std::vector<std::string> vAtom; //Atom hit
        std::vector<double> vRecoil;    // Recoil energy

        int nRecord;

        ion = col.GetIonName();
        mass_number = col.GetMassNumber();
        ene0 = col.GetIonEnergy();

        auto x = col.GetDepth();
        nRecord = 2 + x.size();
        vX.resize(nRecord);
        vX.front() = 0;
        std::copy(x.begin(), x.end(), vX.begin() + 1);
        vX.back() = rng.GetDepth();

        auto y = col.GetLateralY();
        vY.resize(nRecord);
        vY.front() = 0;
        std::copy(y.begin(), y.end(), vY.begin() + 1);
        vY.back() = rng.GetLateralY();

        auto z = col.GetLateralZ();
        vZ.resize(nRecord);
        vZ.front() = 0;
        std::copy(z.begin(), z.end(), vZ.begin() + 1);
        vZ.back() = rng.GetLateralZ();

        auto ene = col.GetEnergy();
        vEne.resize(nRecord);
        vEne.front() = rng.GetIncidentEnergy();
        std::copy(ene.begin(), ene.end(), vEne.begin() + 1);
        vEne.back() = 0;

        auto atom = col.GetAtomHit();
        vAtom.resize(nRecord);
        vAtom.front() = "";
        std::copy(atom.begin(), atom.end(), vAtom.begin() + 1);
        vAtom.back() = "";

        auto recoil = col.GetRecoilEnergy();
        vRecoil.resize(nRecord);
        vRecoil.front() = 0;
        std::copy(recoil.begin(), recoil.end(), vRecoil.begin() + 1);
        vRecoil.back() = 0;

        vdX = CalculateDifference(vX);
        vdY = CalculateDifference(vY);
        vdZ = CalculateDifference(vZ);
        vdEne = CalculateDifference(vEne);

        for (int collisionID = 0; collisionID < nRecord; ++collisionID)
        {

            CollisionRecord rec;

            rec.SetTrackID(trackID);
            rec.SetCollisionID(collisionID);
            rec.SetIncidentEnergy(vEne.at(collisionID));
            rec.SetIncidentIon(ion);
            rec.SetMassNumber(mass_number);
            rec.SetRecoilIon(vAtom.at(collisionID));
            rec.SetRecoilEnergy(vRecoil.at(collisionID));

            double X, Y, Z;
            X = vX.at(collisionID);
            Y = vY.at(collisionID);
            Z = vZ.at(collisionID);
            rec.SetPosition(X, Y, Z);

            double dX0, dY0, dZ0; //direction before collision
            if (collisionID == 0) // assume injection along X-axis
            {
                dX0 = 1;
                dY0 = 0;
                dZ0 = 0;
            }
            else
            {
                dX0 = vdX.at(collisionID - 1);
                dY0 = vdY.at(collisionID - 1);
                dZ0 = vdZ.at(collisionID - 1);
            }
            rec.SetIncidentDirection(dX0, dY0, dZ0);

            double dX = vdX.at(collisionID);
            double dY = vdY.at(collisionID);
            double dZ = vdZ.at(collisionID);
            rec.SetScatteringDirection(dX, dY, dZ);

            double dr = std::sqrt(dX * dX + dY * dY + dZ * dZ);
            rec.SetDistanceToNextCollision(dr);
            // Subtract energy loss due to atomic collision
            double de = std::abs(vdEne.at(collisionID)) - vRecoil.at(collisionID);
            if (de < 0) //Some times dE (= E-E'-E-r) < 0.
                de = 0;
            rec.SetEnergyLoss(de);

            err = sqlite3_bind_int(query, sqlite3_bind_parameter_index(query, ":TID"), rec.GetTrackID());
            err = sqlite3_bind_int(query, sqlite3_bind_parameter_index(query, ":CID"), rec.GetCollisionID());
            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":EIN"), rec.GetIncidentEnergy());

            err = sqlite3_bind_text(query, sqlite3_bind_parameter_index(query, ":IION"), rec.GetIncidentIon().c_str(), -1, nullptr);
            err = sqlite3_bind_int(query, sqlite3_bind_parameter_index(query, ":AION"), rec.GetMassNumber());
            err = sqlite3_bind_text(query, sqlite3_bind_parameter_index(query, ":RION"), rec.GetRecoilIon().c_str(), -1, nullptr);
            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":ERC"), rec.GetRecoilEnergy());

            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":X"), rec.GetPosition().fX);
            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":Y"), rec.GetPosition().fY);
            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":Z"), rec.GetPosition().fZ);
            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":DX0"), rec.GetIncidentDirection().fX);
            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":DY0"), rec.GetIncidentDirection().fY);
            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":DZ0"), rec.GetIncidentDirection().fZ);
            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":DX1"), rec.GetScatteringDirection().fX);
            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":DY1"), rec.GetScatteringDirection().fY);
            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":DZ1"), rec.GetScatteringDirection().fZ);
            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":DR"), rec.GetDistanceToNextCollision());
            err = sqlite3_bind_double(query, sqlite3_bind_parameter_index(query, ":DE"), rec.GetEnergyLoss());

            err = sqlite3_step(query);
            // if (err != SQLITE_OK)
            // {
            //     sqlite3_close(pDB);
            //     sqlite3_finalize(query);
            //     pDB = nullptr;
            //     throw std::runtime_error("logic error " + std::to_string(err));
            // }
            err = sqlite3_reset(query);
            //std::cout << rec.GetTrackID() << " - " << rec.GetCollisionID() << std::endl;
        }

        std::cout << trackID << std::endl;

        ++trackID;
    }

    sqlite3_finalize(query);
    sqlite3_close(pDB);

    return true;
}