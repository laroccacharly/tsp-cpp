#include <fmt/core.h>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <Highs.h>
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

TSPSolution solveTSPWithHiGHS(const TSPInstance& instance) {
    const size_t n = instance.cities.size();
    
    // Create HiGHS model
    HighsModel model;
    
    // Number of variables: n*n for x_ij + (n-1) for u_i (MTZ variables)
    const size_t num_x_vars = n * n;
    const size_t num_u_vars = n - 1;  // u_0 is not needed (reference city)
    const size_t total_vars = num_x_vars + num_u_vars;
    
    // Number of constraints: 2*n (degree constraints) + n*(n-1) (MTZ constraints)
    const size_t num_degree_constraints = 2 * n;
    const size_t num_mtz_constraints = n * (n - 1);
    const size_t total_constraints = num_degree_constraints + num_mtz_constraints;
    
    // Set up model dimensions
    model.lp_.num_col_ = total_vars;
    model.lp_.num_row_ = total_constraints;
    model.lp_.sense_ = ObjSense::kMinimize;
    model.lp_.offset_ = 0;
    
    // Initialize objective coefficients (col_cost_)
    model.lp_.col_cost_.resize(total_vars, 0.0);
    
    // Set costs for x_ij variables (first n*n variables)
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            size_t var_index = i * n + j;
            if (i == j) {
                model.lp_.col_cost_[var_index] = 0.0;  // No cost for staying at same city
            } else {
                model.lp_.col_cost_[var_index] = instance.distanceMatrix.getDistance(i, j);
            }
        }
    }
    // u_i variables have cost 0 (already initialized)
    
    // Set variable bounds
    model.lp_.col_lower_.resize(total_vars);
    model.lp_.col_upper_.resize(total_vars);
    
    // x_ij variables: bounds [0,1]
    for (size_t i = 0; i < num_x_vars; ++i) {
        model.lp_.col_lower_[i] = 0.0;
        model.lp_.col_upper_[i] = 1.0;
    }
    
    // u_i variables: bounds [0, n-1]
    for (size_t i = num_x_vars; i < total_vars; ++i) {
        model.lp_.col_lower_[i] = 0.0;
        model.lp_.col_upper_[i] = static_cast<double>(n - 1);
    }
    
    // Set integrality constraints
    model.lp_.integrality_.resize(total_vars);
    
    // x_ij variables are binary
    for (size_t i = 0; i < num_x_vars; ++i) {
        model.lp_.integrality_[i] = HighsVarType::kInteger;
    }
    
    // u_i variables are continuous
    for (size_t i = num_x_vars; i < total_vars; ++i) {
        model.lp_.integrality_[i] = HighsVarType::kContinuous;
    }
    
    // Set up constraint bounds
    model.lp_.row_lower_.resize(total_constraints);
    model.lp_.row_upper_.resize(total_constraints);
    
    // Degree constraints: each should equal 1
    for (size_t i = 0; i < num_degree_constraints; ++i) {
        model.lp_.row_lower_[i] = 1.0;
        model.lp_.row_upper_[i] = 1.0;
    }
    
    // MTZ constraints: should be <= n-1
    for (size_t i = num_degree_constraints; i < total_constraints; ++i) {
        model.lp_.row_lower_[i] = -1.0e30;  // No lower bound
        model.lp_.row_upper_[i] = static_cast<double>(n - 1);
    }
    
    // Set up constraint matrix (column-wise format)
    model.lp_.a_matrix_.format_ = MatrixFormat::kColwise;
    
    std::vector<int> a_start;
    std::vector<int> a_index;
    std::vector<double> a_value;
    
    a_start.reserve(total_vars + 1);
    
    // Build constraint matrix column by column
    
    // Columns for x_ij variables
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            a_start.push_back(static_cast<int>(a_index.size()));
            
            size_t var_index = i * n + j;
            
            if (i != j) {  // Only add constraints for i != j
                // Outgoing constraint for city i (row i)
                a_index.push_back(static_cast<int>(i));
                a_value.push_back(1.0);
                
                // Incoming constraint for city j (row n + j)
                a_index.push_back(static_cast<int>(n + j));
                a_value.push_back(1.0);
                
                // MTZ constraints: for each pair (k,l) where k!=0, l!=0, k!=l
                // u_k - u_l + n*x_kl <= n-1
                // This becomes: -u_l + u_k + n*x_kl <= n-1
                if (i != 0 && j != 0) {
                    size_t mtz_constraint_idx = num_degree_constraints + i * (n - 1) + (j > 0 ? j - 1 : j);
                    a_index.push_back(static_cast<int>(mtz_constraint_idx));
                    a_value.push_back(static_cast<double>(n));
                }
            }
        }
    }
    
    // Columns for u_i variables (i = 1, 2, ..., n-1)
    for (size_t i = 1; i < n; ++i) {  // u_0 is not included
        a_start.push_back(static_cast<int>(a_index.size()));
        
        // u_i appears in MTZ constraints
        for (size_t j = 1; j < n; ++j) {
            if (i != j) {
                // Constraint: u_i - u_j + n*x_ij <= n-1
                // u_i coefficient: +1 in constraint for (i,j)
                size_t mtz_constraint_idx = num_degree_constraints + i * (n - 1) + (j > 0 ? j - 1 : j);
                a_index.push_back(static_cast<int>(mtz_constraint_idx));
                a_value.push_back(1.0);
                
                // u_i coefficient: -1 in constraint for (j,i)
                mtz_constraint_idx = num_degree_constraints + j * (n - 1) + (i > 0 ? i - 1 : i);
                a_index.push_back(static_cast<int>(mtz_constraint_idx));
                a_value.push_back(-1.0);
            }
        }
    }
    
    // Final entry in a_start
    a_start.push_back(static_cast<int>(a_index.size()));
    
    model.lp_.a_matrix_.start_ = a_start;
    model.lp_.a_matrix_.index_ = a_index;
    model.lp_.a_matrix_.value_ = a_value;
    
    // Create and run HiGHS solver
    Highs highs;
    highs.setOptionValue("output_flag", false);  // Suppress output
    
    HighsStatus return_status = highs.passModel(model);
    if (return_status != HighsStatus::kOk) {
        throw std::runtime_error("Failed to pass model to HiGHS");
    }
    
    return_status = highs.run();
    if (return_status != HighsStatus::kOk) {
        throw std::runtime_error("Failed to solve TSP with HiGHS");
    }
    
    const HighsModelStatus& model_status = highs.getModelStatus();
    if (model_status != HighsModelStatus::kOptimal) {
        throw std::runtime_error("HiGHS did not find optimal solution");
    }
    
    // Extract solution
    const HighsSolution& solution = highs.getSolution();
    
    // Build tour from x_ij values
    TSPSolution tsp_solution;
    tsp_solution.tour.clear();
    
    // Start from city 0
    int current_city = 0;
    tsp_solution.tour.push_back(current_city);
    
    std::vector<bool> visited(n, false);
    visited[current_city] = true;
    
    // Follow the tour
    for (size_t step = 1; step < n; ++step) {
        int next_city = -1;
        
        // Find the next city where x_current_city,next_city = 1
        for (size_t j = 0; j < n; ++j) {
            if (visited[j]) continue;
            
            size_t var_index = current_city * n + j;
            if (solution.col_value[var_index] > 0.5) {  // x_ij = 1
                next_city = static_cast<int>(j);
                break;
            }
        }
        
        if (next_city == -1) {
            throw std::runtime_error("Invalid tour found in HiGHS solution");
        }
        
        tsp_solution.tour.push_back(next_city);
        visited[next_city] = true;
        current_city = next_city;
    }
    
    // Return to starting city
    tsp_solution.tour.push_back(0);
    
    return tsp_solution;
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
    
    // Solve with nearest neighbor heuristic
    auto nn_solution = solveTSP(instance);
    printSolution(nn_solution, instance);
    
    // Solve with HiGHS (optimal solution)
    try {
        fmt::print("\n" + std::string(50, '=') + "\n");
        auto optimal_solution = solveTSPWithHiGHS(instance);
        
        fmt::print("\nTSP Solution (HiGHS Optimal):\n");
        fmt::print("Tour: ");
        for (size_t i = 0; i < optimal_solution.tour.size(); ++i) {
            fmt::print("{}", optimal_solution.tour[i]);
            if (i < optimal_solution.tour.size() - 1) {
                fmt::print(" -> ");
            }
        }
        fmt::print("\n");
        
        fmt::print("City names: ");
        for (size_t i = 0; i < optimal_solution.tour.size(); ++i) {
            int cityId = optimal_solution.tour[i];
            fmt::print("{}", instance.cities[cityId].name);
            if (i < optimal_solution.tour.size() - 1) {
                fmt::print(" -> ");
            }
        }
        fmt::print("\n");
        
        double optimal_distance = calculateTourDistance(optimal_solution, instance);
        double nn_distance = calculateTourDistance(nn_solution, instance);
        
        fmt::print("Total distance: {:.2f}\n", optimal_distance);
        fmt::print("\nComparison:\n");
        fmt::print("Nearest Neighbor distance: {:.2f}\n", nn_distance);
        fmt::print("Optimal distance: {:.2f}\n", optimal_distance);
        fmt::print("Improvement: {:.2f} ({:.1f}%)\n", 
                   nn_distance - optimal_distance, 
                   ((nn_distance - optimal_distance) / nn_distance) * 100.0);
                   
    } catch (const std::exception& e) {
        fmt::print("Error solving TSP with HiGHS: {}\n", e.what());
    }
    
    return 0;
}