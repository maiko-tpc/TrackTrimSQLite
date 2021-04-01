#include <iostream>

#include "TRIM2SQLite.hpp"

int main(int argc, char *argv[])
{

    if (argc < 3)
    {
        std::cerr << argv[0] << " [input_directory] [output_name]" << std::endl;
        return 1;
    }
    TRIM2SQLite t2s;
    // Example
    // t2s.MakeSQLiteFile("../input/TRIM/1000/3H/10", "hoge.sqlite");
    t2s.MakeSQLiteFile(argv[1], argv[2]);
    return 0;
}