#include "cities.h"

std::vector<City> createCities() { 
    std::vector<City> cities;
    cities.push_back({0, "Paris", 48.8566, 2.3522});
    cities.push_back({1, "London", 51.5074, -0.1278});
    cities.push_back({2, "Berlin", 52.5200, 13.4050});
    cities.push_back({3, "Madrid", 40.4168, -3.7038});
    cities.push_back({4, "Rome", 41.9028, 12.4964});
    return cities;
} 