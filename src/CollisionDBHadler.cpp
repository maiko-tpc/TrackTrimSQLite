#include "CollisionDBHandler.hpp"
#include "my_sqlite3util.hpp"

int CollisionDBHandler::GetNumberOfTracks(int _trackID)
{

    std::string sQuery = "SELECT COUNT(DISTINCT track_id) FROM collisions ";

    if (_trackID >= 0) //track ID is specified
    {
        sQuery += "WHERE track_id = " + std::to_string(_trackID);
    }

    sQuery += ";";

    int n;
    try
    {
        n = ExecuteCountQuery(sQuery);
    }
    catch (std::runtime_error e)
    {
        throw std::runtime_error("CollisionDBHandler::GetNumberOfTracks():" + std::string(e.what()));
    }

    return n;
}

int CollisionDBHandler::GetNumberOfCollisions(int _trackID)
{
    std::string sQuery = "SELECT COUNT(*) FROM collisions ";

    if (_trackID >= 0) //track ID is specified
    {
        sQuery += "WHERE track_id = " + std::to_string(_trackID);
    }

    sQuery += ";";

    int n;
    try
    {
        n = ExecuteCountQuery(sQuery);
    }
    catch (std::runtime_error e)
    {
        throw std::runtime_error("CollisionDBHandler::GetNumberOfCollisions():" + std::string(e.what()));
    }

    return n;
}

std::vector<TRIM2SQLite::CollisionRecord>
CollisionDBHandler::GetTrack(int _trackID)
{
    std::string sQuery = "SELECT * FROM collisions ";
    sQuery += "WHERE track_id = " + std::to_string(_trackID) + "　ORDER BY collision_id;";
    std::vector<TRIM2SQLite::CollisionRecord> ret;
    ExecuteSelectQuery(sQuery, ret);
    return ret;
}

std::vector<TRIM2SQLite::CollisionRecord>
CollisionDBHandler::GetCollisions(const std::string &_constraint,
                                  int _limit)
{
    std::string sQuery = "SELECT * FROM collisions ";

    if (_constraint != "")
        sQuery += "WHERE " + _constraint + " ";
    if (_limit > 0)
        sQuery += " LIMIT " + std::to_string(_limit);
    //std::cout << sQuery << std::endl;
    std::vector<TRIM2SQLite::CollisionRecord> ret;
    ExecuteSelectQuery(sQuery, ret);
    return ret;
}

int CollisionDBHandler::ExecuteCountQuery(const std::string &_query)
{
    sqlite3_stmt *query;
    auto err = sqlite3_prepare_v2(fpDB,
                                  _query.c_str(),
                                  -1, &query, nullptr);

    //    std::cout << "prepare -> " << err << std::endl;
    if (err != SQLITE_OK)
    {
        sqlite3_finalize(query);
        throw std::runtime_error("ExecuteCountQuery : Query preparation error.");
    }

    err = sqlite3_step(query);
    //  std::cout << "step    -> " << err << " " << _query << std::endl;

    int n = 0;
    if (err == SQLITE_ROW)
    {
        n = sqlite3_column_int(query, 0);
        sqlite3_finalize(query);
    }
    else
    {
        sqlite3_finalize(query);
        throw std::runtime_error("ExecuteCountQuery : Query execution error.");
    }

    return n;
}

// _query = SELECT * from collisions (+constraint)
int CollisionDBHandler::ExecuteSelectQuery(const std::string &_query,
                                           std::vector<TRIM2SQLite::CollisionRecord> &_rec)
{

    _rec.clear();

    sqlite3_stmt *query;
    auto err = sqlite3_prepare_v2(fpDB,
                                  _query.c_str(),
                                  -1, &query, nullptr);

    //std::cout << "prepare -> " << err << std::endl;
    if (err != SQLITE_OK)
    {
        sqlite3_finalize(query);
        throw std::runtime_error("ExecuteSelectQuery : Query preparation error.");
    }

    err = sqlite3_step(query);
    //std::cout << "step    -> " << err << " " << _query << std::endl;

    int n = 0;
    while (err == SQLITE_ROW)
    {
        _rec.push_back(TRIM2SQLite::CollisionRecord());

        int track_id = sqlite3_column_int(query, 0);
        int collision_id = sqlite3_column_int(query, 1);

        double e_inc = sqlite3_column_double(query, 2);
        std::string incident_ion = my_sqlite3_column_string(query, 3);
        int incident_ion_mass = sqlite3_column_int(query, 4);
        std::string recoil_ion = my_sqlite3_column_string(query, 5);
        double e_rec = sqlite3_column_double(query, 6);

        double x = sqlite3_column_double(query, 7);
        double y = sqlite3_column_double(query, 8);
        double z = sqlite3_column_double(query, 9);

        double dx0 = sqlite3_column_double(query, 10);
        double dy0 = sqlite3_column_double(query, 11);
        double dz0 = sqlite3_column_double(query, 12);

        double dx1 = sqlite3_column_double(query, 13);
        double dy1 = sqlite3_column_double(query, 14);
        double dz1 = sqlite3_column_double(query, 15);

        double dr = sqlite3_column_double(query, 16);
        double de = sqlite3_column_double(query, 17);

        _rec.back().SetTrackID(track_id);
        _rec.back().SetCollisionID(collision_id);
        _rec.back().SetIncidentEnergy(e_inc);
        _rec.back().SetIncidentIon(incident_ion);
        _rec.back().SetMassNumber(incident_ion_mass);
        _rec.back().SetRecoilIon(recoil_ion);
        _rec.back().SetRecoilEnergy(e_rec);
        _rec.back().SetPosition(x, y, z);
        _rec.back().SetIncidentDirection(dx0, dy0, dz0);
        _rec.back().SetScatteringDirection(dx1, dy1, dz1);
        _rec.back().SetDistanceToNextCollision(dr);
        _rec.back().SetEnergyLoss(de);

        err = sqlite3_step(query);
    }

    sqlite3_finalize(query);

    return _rec.size();
}