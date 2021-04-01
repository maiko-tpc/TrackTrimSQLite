#pragma once 

#include <iostream>
#include <string>
#include <sstream>


#include <sqlite3.h>

std::string my_sqlite3_column_string(sqlite3_stmt * _query, int i){
    std::ostringstream tmp;
    tmp << sqlite3_column_text(_query, i);
    return tmp.str();
}
