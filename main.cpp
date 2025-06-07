#include <fmt/core.h>
#include <vector>
#include <cmath>
#include <limits>
#include "src/cities.h"

struct DistanceMatrix {
    std::vector<double> distances;
    size_t size;
    
    DistanceMatrix(size_t n) : size(n), distances(n * n, 0.0) {}
    
    double getDistance(size_t i, size_t j) const {
        return distances[i * size + j];
    }
    
    void setDistance(size_t i, size_t j, double distance) {
        distances[i * size + j] = distance;
    }
    
    void computeFromCities(const std::vector<City>& cities) {
        for (size_t i = 0; i < cities.size(); ++i) {
            for (size_t j = 0; j < cities.size(); ++j) {
                if (i == j) {
                    setDistance(i, j, 0.0);
                } else {
                    double dx = cities[i].x - cities[j].x;
                    double dy = cities[i].y - cities[j].y;
                    double distance = std::sqrt(dx * dx + dy * dy);
                    setDistance(i, j, distance);
                }
            }
        }
    }
};

struct TSPInstance {
    std::vector<City> cities;
    DistanceMatrix distanceMatrix;

    TSPInstance(const std::vector<City>& cities) : cities(cities), distanceMatrix(cities.size()) {
        distanceMatrix.computeFromCities(cities);
    }
}; 

struct TSPSolution { 
    std::vector<int> tour;
}; 

void printCity(const City& city) {
    fmt::print("{} (x: {}, y: {})\n", city.id, city.name, city.x, city.y);
}

TSPSolution solveTSP(const TSPInstance& instance) {
    const size_t numCities = instance.cities.size();
    
    TSPSolution solution;
    solution.tour.reserve(numCities + 1); // +1 for return to start
    
    std::vector<bool> visited(numCities, false);
    
    // Start at city 0
    int currentCityId = 0;
    solution.tour.push_back(currentCityId);
    visited[currentCityId] = true;
    
    // Visit remaining cities using nearest neighbor
    for (size_t step = 1; step < numCities; ++step) {
        double minDistance = std::numeric_limits<double>::max();
        int nearestCityId = -1;
        
        // Find nearest unvisited city
        for (size_t candidate = 0; candidate < numCities; ++candidate) {
            if (visited[candidate]) {
                continue;
            }
            
            double distance = instance.distanceMatrix.getDistance(currentCityId, candidate);
            if (distance < minDistance) {
                minDistance = distance;
                nearestCityId = static_cast<int>(candidate);
            }
        }
        
        // Move to nearest city
        currentCityId = nearestCityId;
        solution.tour.push_back(currentCityId);
        visited[currentCityId] = true;
    }
    
    // Return to starting city to complete the tour
    solution.tour.push_back(0);
    
    return solution;
}

double calculateTourDistance(const TSPSolution& solution, const TSPInstance& instance) {
    if (solution.tour.size() < 2) {
        return 0.0;
    }
    
    double totalDistance = 0.0;
    for (size_t i = 0; i < solution.tour.size() - 1; ++i) {
        int from = solution.tour[i];
        int to = solution.tour[i + 1];
        totalDistance += instance.distanceMatrix.getDistance(from, to);
    }
    
    return totalDistance;
}

void printSolution(const TSPSolution& solution, const TSPInstance& instance) {
    fmt::print("\nTSP Solution (Nearest Neighbor):\n");
    fmt::print("Tour: ");
    for (size_t i = 0; i < solution.tour.size(); ++i) {
        fmt::print("{}", solution.tour[i]);
        if (i < solution.tour.size() - 1) {
            fmt::print(" -> ");
        }
    }
    fmt::print("\n");
    
    fmt::print("City names: ");
    for (size_t i = 0; i < solution.tour.size(); ++i) {
        int cityId = solution.tour[i];
        fmt::print("{}", instance.cities[cityId].name);
        if (i < solution.tour.size() - 1) {
            fmt::print(" -> ");
        }
    }
    fmt::print("\n");
    
    double distance = calculateTourDistance(solution, instance);
    fmt::print("Total distance: {:.2f}\n", distance);
}

int main()
{
    auto instance = TSPInstance(createCities());
    
    fmt::print("TSP Instance:\n");
    for (const auto& city : instance.cities) {
        printCity(city);
    }
    
    auto solution = solveTSP(instance);
    printSolution(solution, instance);
    
    return 0;
}