#pragma once
#include <iostream>
#include <vector>
#include <string>

#include <sqlite3.h>

#include "TRIM2SQLite.hpp"

class CollisionDBHandler
{
public:
    CollisionDBHandler(const std::string &_filename) : fpDB(nullptr)
    {
        auto err = sqlite3_open_v2(_filename.c_str(), &fpDB,
                                   SQLITE_OPEN_READONLY, nullptr);
        if (err != SQLITE_OK)
        {
            sqlite3_close(fpDB);
            fpDB = nullptr;
            throw std::runtime_error("CollisionDBHandler: DB file not found.");
        }
    };

    ~CollisionDBHandler()
    {
        sqlite3_close(fpDB);
    };

    int GetNumberOfTracks(int _trackID = -1);
    int GetNumberOfCollisions(int _trackID = -1);
    std::vector<TRIM2SQLite::CollisionRecord> GetTrack(int _trackID);

    std::vector<TRIM2SQLite::CollisionRecord>
    GetCollisions(const std::string &_constraint,
                  int _limit = -1);

private:
    int ExecuteCountQuery(const std::string &_query);
    int ExecuteSelectQuery(const std::string &_query,
                           std::vector<TRIM2SQLite::CollisionRecord> &_rec);
    sqlite3 *fpDB;
};
