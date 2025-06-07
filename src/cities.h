#pragma once

#include <vector>
#include <string>

struct City {
    int id;
    std::string name;
    double x;
    double y;
};

std::vector<City> createCities(); 