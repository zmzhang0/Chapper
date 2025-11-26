#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <cmath>
#include <queue>
#include <vector>

#include "SFML/Graphics.hpp"
#include "VBM/field.hpp"
#include "VBM/flat_hash_map.hpp"

namespace vbm {

    struct Node {
        int x, y, z;
        double f;
        bool operator<(const Node& other) const { return f > other.f; }
    };

    struct point {
        int x, y, z;
    };

    struct vector3f {
        float x, y, z;
    };


    sf::Color getColor(double value);

    class Solver {
    public:
        using Map = ska::flat_hash_map<size_t, double>;

        explicit Solver(std::shared_ptr< Field<double>>& sharedVisibilityField, std::shared_ptr< Field<double>>& sharedSpeedField_);

        // Deconstructor
        ~Solver() = default;

        const Field<size_t>& visibilityBasedSolver();
        const Field<double>& computeDistanceFunction();
        const Field<double>& computeDistanceFunction(const Field<vector3f>& gradientField);


        // Get global iterator (number of iterations that had to be completed)
        inline int getNbOfIterations() const { return nb_of_iterations_; };

        // Get number of sources
        inline int getNbOfSources() const { return initial_frontline.size() / 3; };


        // Get the Distance Function
        Field<double> getDistanceFunction(){ return gScore_; };

        // get the CameFrom field
        Field<size_t> getCameFromField(){ return cameFrom_; };

        // get the light sources
        point* getLightSources() const {
            return lightSources_.get();
        }

    private:
        void reset();
        void reconstructPath(const Node& current, const std::string& methodName);
        inline int indexAt(const size_t x, const size_t y, const size_t z) const {
            return x * ny_ * nz_ + y * nz_ + z;
        };
        inline point coordinatesAt(const size_t index) const {
            const int x = index / (ny_ * nz_);
            const int y = index / nz_ % ny_;
            const int z = index % nz_;
            return { x, y, z };
        }

        /*!
         * @brief queues sources
         */
        inline void queuePotentialSources(std::vector<size_t>& potentialSources,
            const int neighbour_x,
            const int neighbour_y, 
            const int neighbour_z) const;

        /*!
         * @brief gets distances
         */
        inline void getPotentialDistances(
            const std::vector<size_t>& potentialSources,
            std::vector<std::pair<double, size_t>>& potentialDistances,
            const int neighbour_x, const int neighbour_y, const int neighbour_z);

        /*!
         * @brief gets distances
         */
        void getPotentialDistancesField(
            const std::vector<size_t>& potentialSources,
            std::vector<std::pair<double, size_t>>& potentialDistances,
            const int neighbour_x, const int neighbour_y, const int neighbour_z);

        inline void getPotentialDistancesSpeedField(
            const std::vector<size_t>& potentialSources,
            std::vector<std::pair<double, size_t>>& potentialDistances,
            const int neighbour_x, const int neighbour_y, const int neighbour_z);

        inline const double evaluateDistance(const int source_x,
            const int source_y, const int source_z,
            const int target_x, const int target_y, const int target_z) const {
            return sqrt((double)(source_x - target_x) * (source_x - target_x) +
                (source_y - target_y) * (source_y - target_y) +
                (source_z - target_z) * (source_z - target_z));
        };

        inline const double evaluateDistance(const int source_x,
            const int source_y, const int source_z,
            const int target_x, const int target_y, const int target_z, const vector3f& gradient) const {
            
            // 计算欧拉距离

            vector3f dir = vector3f(target_x - source_x,
                target_y - source_y,
                target_z - source_z);

            double euclidean_distance = std::sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);

            dir.x /= euclidean_distance;
            dir.y /= euclidean_distance;
            dir.z /= euclidean_distance;

            vector3f grad = vector3f(gradient.x, gradient.y, gradient.z);
            double grad_length = std::sqrt(grad.x * grad.x + grad.y * grad.y + grad.z * grad.z);

            if (grad_length < 1e-5) {
                return euclidean_distance; // 没有方向信息，退化为各向同性
            }

            grad.x /= grad_length;
            grad.y /= grad_length;
            grad.z /= grad_length;

            double dot_product = dir.x * grad.x + dir.y * grad.y + dir.z * grad.z;
            double alpha = 10.0;  // 越大越方向敏感 *********超参数*********
            double penalty = 1.0 + alpha * (1.0 - dot_product);

            return euclidean_distance * penalty;
        }

        inline const double evaluateManhattanDistance(const int source_x,
            const int source_y, const int source_z,
            const int target_x, const int target_y, const int target_z) const {
            return (double)(abs(source_x - target_x) +
                abs(source_y - target_y) +
                abs(source_z - target_z));
        };


        inline const double evaluateDistanceSpeedField(const int source_x,
            const int source_y, const int source_z,
            const int target_x, const int target_y, const int target_z) const {
            return sqrt((double)(source_x - target_x) * (source_x - target_x) +
                (source_y - target_y) * (source_y - target_y) +
                (source_z - target_z) * (source_z - target_z)) *
                sharedSpeedField_->get(target_x, target_y, target_z);
        };
        inline void createNewPivot(const int x, const int y, const int z,
            const int neighbour_x, const int neighbour_y, const int neighbour_z);
        
        void saveResults(const std::vector<point>& path,
            const std::string& methodName) const;
       
        void saveVisibilityBasedSolverImage(const Field<double>& gScore) const;
        void saveDistanceFunctionImage(const Field<double>& gScore) const;
        void saveCameFromImage(const Field<size_t>& cameFrom) const;

        /*!
         * @brief Updates accessibility/visibility to a point using PDE advection.
         * @param [in] lightSourceNumber number of the lightsource whose visibility of
         * the point we are checking.
         * @param [in] lightSource_x x position of the lightsource.
         * @param [in] lightSource_y y position of the lightsource.
         * @param [in] lightSource_z z position of the lightsource.
         * @param [in] x position of our queried pixel.
         * @param [in] y position of our queried pixel.
         * @param [in] z position of our queried pixel.
         */
        void updatePointVisibility(const size_t lightSourceNumber,
            const int lightSource_x, const int lightSource_y, const int lightSource_z,
            const int x, const int y, const int z);

        inline const int hashFunction(const int x, const int y, const int z,
            const int lightSourceNumber) const {
            const int prime1 = 73856093;
            const int prime2 = 19349663;
            const int prime3 = 83492791;
            const int prime4 = 49979687;

            const long long key = (x * prime1) ^ (y * prime2) ^ (z * prime3) ^ (lightSourceNumber * prime4);

            const int modulus = 2147483647; // 增大模数以减少冲突
            return key % modulus;
        }


        /*!
         * @brief Checks if requested cell is in grid
         * @param [in] x x position of the cell
         * @param [in] y y position of the cell
         * @param [in] z z position of the cell
         * @return true if cell is in grid, false otherwise
         */
        inline const bool isValid(const int x, const int y, const int z) const {
            return x >= 0 && x < nx_ && y >= 0 && y < ny_ && z >= 0 && z < nz_;
        }

        std::shared_ptr<Field<double>> sharedVisibilityField_;
        std::shared_ptr<Field<double>> sharedSpeedField_;
        std::vector<int> initial_frontline;

        Field<double> gScore_;
        Field<double> fScore_;
        Field<size_t> cameFrom_;
        Field<bool> inOpenSet_;
        Field<bool> updated_;

        std::unique_ptr<point[]> lightSources_;

        const double lightStrength_ = 1.0;
        int nb_of_iterations_ = 0;
        size_t nb_of_sources_ = 0;

        // Neighbours in 3D
        // [1 0; 0 1; -1 0; 0 -1; 1 1; -1 1; -1 -1; 1 -1] flattened out
        const int neighbours_[78] = {
            1,  0,  0,  -1,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,  1,  0,  0, -1,  // 6 direct neighbors
            1,  1,  0,  -1,  1,  0,  1, -1,  0, -1, -1,  0,  // 4 diagonal neighbors on the xy plane
            1,  0,  1,  -1,  0,  1,  1,  0, -1, -1,  0, -1,  // 4 diagonal neighbors on the xz plane
            0,  1,  1,   0, -1,  1,  0,  1, -1,  0, -1, -1,  // 4 diagonal neighbors on the yz plane
            1,  1,  1,  -1,  1,  1,  1, -1,  1, -1, -1,  1,  // 4 diagonal neighbors on the x, y, and z axes (positive direction)
            1,  1, -1,  -1,  1, -1,  1, -1, -1, -1, -1, -1};



        size_t nx_, ny_, nz_;  // 三维网格大小
        // Heap openSet_;
        std::unique_ptr<std::priority_queue<Node>> openSet_;

        double visibilityThreshold_ = 0.5;

        // Flat hash maps
        Map visibilityHashMap_;

        bool Outsilent = false;
        bool Outdetails= false;

        // Unique pointer to image holder
        std::unique_ptr<sf::Image> uniqueLoadedImage_;
    };

} // namespace vbm
#endif // SOLVER_HPP
