#include "VBM/solver.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace vbm {

    template <typename T> auto durationInMicroseconds(T start, T end) {
        return std::chrono::duration_cast<std::chrono::microseconds>(end - start)
            .count();
    }

    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    sf::Color getColor(double value) {
        // jet colormap for SFML visualization/plot
        const int color_index = 255 * value;
        double r, g, b;
        if (color_index < 32) {
            r = 0;
            g = 0;
            b = 0.5156 + 0.0156 * color_index;
        }
        else if (color_index < 96) {
            r = 0;
            g = 0.0156 + 0.9844 * (color_index - 32.0) / 64;
            b = 1;
        }
        else if (color_index < 158) {
            r = 0.0156 + (color_index - 96.0) / 64;
            g = 1;
            b = 0.9844 - (color_index - 96.0) / 64;
        }
        else if (color_index < 223) {
            r = 1;
            g = 1 - (color_index - 158.0) / 65;
            b = 0;
        }
        else {
            r = (2 - (color_index - 223.0) / 32) / 2.0;
            g = 0;
            b = 0;
        }
        return sf::Color(static_cast<sf::Uint8>(r * 255),
            static_cast<sf::Uint8>(g * 255),
            static_cast<sf::Uint8>(b * 255));
    }

    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    // 需要传递阐参数初始化
    Solver::Solver(std::shared_ptr< Field<double>>& sharedVisibilityField, std::shared_ptr< Field<double>>& sharedSpeedField) {
        // Init sharedVisibilityField_, sharedSpeedField_
        sharedVisibilityField_ = sharedVisibilityField;
        sharedSpeedField_ = sharedSpeedField;

        nx_ = sharedVisibilityField_->nx();
        ny_ = sharedVisibilityField_->ny();
        nz_ = sharedVisibilityField_->nz();

        visibilityThreshold_ = 0.5;

        // Init maps
        reset();
    }

    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    void Solver::reset() {
        gScore_.reset(nx_, ny_, nz_, std::numeric_limits<double>::infinity());
        fScore_.reset(nx_, ny_, nz_, std::numeric_limits<double>::infinity());
        cameFrom_.reset(nx_, ny_, nz_, 0);
        inOpenSet_.reset(nx_, ny_, nz_, false);
        updated_.reset(nx_, ny_, nz_, false);

        lightSources_.reset(new point[nx_ * ny_ * nz_]);

        visibilityHashMap_.clear();
        openSet_.reset();

        // Reserve openSet_
        std::vector<Node> container;
        container.reserve(nx_ * ny_ * nz_);
        std::priority_queue<Node, std::vector<Node>, std::less<Node>> heap(
            std::less<Node>(), std::move(container));
        openSet_ = std::make_unique<std::priority_queue<Node>>(heap);

        nb_of_sources_ = 0;
        nb_of_iterations_ = 0;

        // Reserve hash map
        visibilityHashMap_.reserve(nx_ * ny_ * nz_);
    }

    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    const Field<size_t>& Solver::visibilityBasedSolver() {
        reset();
        auto startTime = std::chrono::high_resolution_clock::now();

        // Init
        double d = 0;
        int x = 0, y = 0, z = 0, neighbour_x = 0, neighbour_y = 0, neighbour_z = 0;

        // Init frontline

        if (nz_ == 1)
        {
            for (int i = 0; i < nx_; i++)
                for (int j = 0; j < ny_; j++) {
                    if (i == 0 || i == nx_ - 1 || j == 0 || j == ny_ - 1) {
                        initial_frontline.push_back(i);
                        initial_frontline.push_back(j);
                        initial_frontline.push_back(0);
                    }
                }
        }
        else
        {
            for (int i = 0; i < nx_; i++)
                for (int j = 0; j < ny_; j++)
                    for (int k = 0; k < nz_; k++) {
                        if (i == 0 || i == nx_ - 1 || j == 0 || j == ny_ - 1 || k == 0 || k == nz_ - 1) {
                            initial_frontline.push_back(i);
                            initial_frontline.push_back(j);
                            initial_frontline.push_back(k);
                        }
                    }
        }


        if (initial_frontline.size() % 3 != 0) {
            std::cout << "###################### ERROR：Visibility-based solver output "
                "######################"
                << std::endl;
            std::cout << "Initial frontline must be of size that is a multiple of 3 "
                "for visibility-based solver"
                << std::endl;

            return cameFrom_;
        }


        for (size_t i = 0; i < initial_frontline.size(); i += 3) {
            x = initial_frontline[i];
            y = initial_frontline[i + 1];
            z = initial_frontline[i + 2];
            // check if starting positions are inside the map
            if (x >= nx_ || y >= ny_ || z >= nz_ || x < 0 || y < 0 || z < 0) {
                std::cout << "###################### ERROR：Visibility-based solver output "
                    "######################"
                    << std::endl;
                std::cout << "At least one of the starting positions is outside the map"
                    << std::endl;
                return cameFrom_;
            }

            if (sharedVisibilityField_->get(x, y, z) < 1) {
                std::cout << "###################### ERROR：Visibility-based solver output "
                    "######################"
                    << std::endl;
                std::cout << "At least one of the starting positions is invalid/occupied"
                    << std::endl;
                return cameFrom_;
            }

            d = 0;
            gScore_(x, y, z) = d;
            updated_(x, y, z) = true;
            cameFrom_(x, y, z) = nb_of_sources_;
            lightSources_[nb_of_sources_] = { x, y, z };

            openSet_->push(Node{ x, y, z, d });
            const auto key = hashFunction(x, y, z, nb_of_sources_);
            visibilityHashMap_[key] = lightStrength_;
            ++nb_of_sources_;
            ++nb_of_iterations_;
        }


        // For queing unique sources from neighbours of neighbour
        std::vector<size_t> potentialSources;
        potentialSources.reserve(10);
        std::vector<std::pair<double, size_t>> potentialDistances;
        potentialDistances.reserve(10);

        double distance = 0;

        while (openSet_->size() > 0) {
            auto& current = openSet_->top();
            x = current.x;
            y = current.y;
            z = current.z;
            openSet_->pop();

            // Expand frontline at current & update neighbours
            for (size_t j = 0; j < 78; j += 3) {
                // NOTE those are always positive
                neighbour_x = x + neighbours_[j];
                neighbour_y = y + neighbours_[j + 1];
                neighbour_z = z + neighbours_[j + 2];

                // Box check
                if (!isValid(neighbour_x, neighbour_y, neighbour_z)) {
                    continue;
                }
                if (updated_(neighbour_x, neighbour_y, neighbour_z)) {
                    continue;
                };
                if (sharedVisibilityField_->get(neighbour_x, neighbour_y, neighbour_z) < 1) {
                    if (0) {// 光线是否可以穿透障碍物
                        gScore_(neighbour_x, neighbour_y, neighbour_z) =
                            gScore_(x, y, z) +
                            evaluateDistanceSpeedField(x, y, z, neighbour_x, neighbour_y, neighbour_z);
                        openSet_->push(Node{ neighbour_x, neighbour_y,neighbour_z,
                                            gScore_(neighbour_x, neighbour_y,neighbour_z) });

                        const auto key = hashFunction(x, y, z, nb_of_sources_);
                        visibilityHashMap_[key] = sharedVisibilityField_->get(neighbour_x, neighbour_y, neighbour_z);
                    }
                    cameFrom_(neighbour_x, neighbour_y, neighbour_z) = cameFrom_(x, y, z);
                    updated_(neighbour_x, neighbour_y, neighbour_z) = true;
                    continue;
                }

                // in case only 1 source so far, no need to queue potential parents
                if (nb_of_sources_ == 1) {
                    potentialSources.clear();
                    potentialSources.push_back(0);
                }
                else {
                    queuePotentialSources(potentialSources, neighbour_x, neighbour_y, neighbour_z);
                }
                getPotentialDistancesField(potentialSources, potentialDistances,
                    neighbour_x, neighbour_y, neighbour_z);

                auto minimum_element =
                    std::min_element(potentialDistances.begin(), potentialDistances.end(),
                        [](const auto& lhs, const auto& rhs) {
                            return lhs.first < rhs.first;
                        });
                distance = minimum_element->first;

                // 找到了该点的光源后，但是所有光源都不可见，则需要创建新的中枢点
                if (distance == std::numeric_limits<double>::infinity()) {
                    createNewPivot(x, y, z, neighbour_x, neighbour_y, neighbour_z);
                }
                else {
                    // use source giving least distance
                    gScore_(neighbour_x, neighbour_y, neighbour_z) = minimum_element->first;
                    cameFrom_(neighbour_x, neighbour_y, neighbour_z) = minimum_element->second;
                }
                openSet_->push(
                    Node{ neighbour_x, neighbour_y,neighbour_z, gScore_(neighbour_x, neighbour_y,neighbour_z) });
                updated_(neighbour_x, neighbour_y, neighbour_z) = true;
                ++nb_of_iterations_;
            }
        };

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto executionDuration = durationInMicroseconds(startTime, stopTime);

        if (!Outsilent) {
            std::cout << std::endl << "###################### Visibility-based solver output "
                "######################"
                << std::endl;
            if (!Outdetails) {
                std::cout << "Execution time in s: " << (double)executionDuration / 1000.0 / CLOCKS_PER_SEC << "s"
                    << std::endl;
                std::cout << "Load factor: " << visibilityHashMap_.load_factor()
                    << std::endl;
                std::cout << "Iterations: " << nb_of_iterations_ << std::endl;
                std::cout << "Nb of sources: " << nb_of_sources_ << std::endl;
            }

            std::cout << "######################################################"
                "######################"
                << std::endl << std::endl;
        }

        saveVisibilityBasedSolverImage(gScore_);
        saveCameFromImage(cameFrom_);
        return cameFrom_;
    }


    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    const Field<double>& Solver::computeDistanceFunction() {
        reset();
        auto startTime = std::chrono::high_resolution_clock::now();

        for (size_t i = 0; i < nx_; ++i) {
            for (size_t j = 0; j < ny_; ++j)
                for (size_t k = 0; k < nz_; ++k) {
                    if (sharedVisibilityField_->get(i, j, k) <= visibilityThreshold_) {
                        initial_frontline.push_back(i);
                        initial_frontline.push_back(j);
                        initial_frontline.push_back(k);
                    }
                }
        }
        // Init
        int x = 0, y = 0, z = 0, neighbour_x = 0, neighbour_y = 0, neighbour_z = 0;
        double g = 0;
        // Fill in data from initial frontline
        for (size_t i = 0; i < initial_frontline.size(); i += 3) {
            x = initial_frontline[i];
            y = initial_frontline[i + 1];
            z = initial_frontline[i + 2];
            g = 0;
            openSet_->push(Node{ x, y,z, g });

            gScore_(x, y, z) = g;
            updated_(x, y, z) = true;
            cameFrom_(x, y, z) = nb_of_sources_;
            lightSources_[nb_of_sources_] = { x, y,z };
            const auto key = hashFunction(x, y, z, nb_of_sources_);
            visibilityHashMap_[key] = lightStrength_;
            ++nb_of_sources_;
        }

        // For queing unique sources from neighbours of neighbour
        std::vector<size_t> potentialSources;
        potentialSources.reserve(10);
        std::vector<std::pair<double, size_t>> potentialDistances;
        potentialDistances.reserve(10);

        double distance = 0;

        while (openSet_->size() > 0) {
            auto& current = openSet_->top();
            x = current.x;
            y = current.y;
            z = current.z;

            openSet_->pop();

            // Expand frontline at current & update neighbours
            for (size_t j = 0; j < 78; j += 3) {
                // NOTE those are always positive
                neighbour_x = x + neighbours_[j];
                neighbour_y = y + neighbours_[j + 1];
                neighbour_z = z + neighbours_[j + 2];

                // Box check
                if (!isValid(neighbour_x, neighbour_y, neighbour_z)) {
                    continue;
                }
                if (updated_(neighbour_x, neighbour_y, neighbour_z)) {
                    continue;
                };

                queuePotentialSources(potentialSources, neighbour_x, neighbour_y, neighbour_z);
                potentialDistances.clear();
                for (int k = 0; k < potentialSources.size(); ++k) {
                    int potentialSource = potentialSources[k];
                    int LS_x = lightSources_[potentialSource].x;
                    int LS_y = lightSources_[potentialSource].y;
                    int LS_z = lightSources_[potentialSource].z;
                    distance = gScore_(LS_x, LS_y, LS_z) +
                        evaluateDistance(LS_x, LS_y, LS_z, neighbour_x, neighbour_y, neighbour_z);
                    potentialDistances.push_back(
                        std::pair<double, int>{distance, potentialSource});
                }

                auto minimum_element =
                    std::min_element(potentialDistances.begin(), potentialDistances.end(),
                        [](const auto& lhs, const auto& rhs) {
                            return lhs.first < rhs.first;
                        });
                distance = minimum_element->first;
                // use source giving least distance
                gScore_(neighbour_x, neighbour_y, neighbour_z) = distance;
                cameFrom_(neighbour_x, neighbour_y, neighbour_z) = minimum_element->second;
                openSet_->push(Node{ neighbour_x, neighbour_y,neighbour_z, distance });
                updated_(neighbour_x, neighbour_y, neighbour_z) = true;
            }
        };

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto executionDuration = durationInMicroseconds(startTime, stopTime);

        if (!Outsilent) {
            std::cout << std::endl << "################# VBM:Distance function computation "
                "output #####################"
                << std::endl;

            if (!Outdetails) {
                std::cout << "Constructed distance function" << std::endl;
                std::cout << "Execution time in s: " << (double)executionDuration / 1000.0 / CLOCKS_PER_SEC << "s"
                    << std::endl;
            }
            std::cout << "####################################################"
                "############################"
                << std::endl << std::endl;
        }

        //saveDistanceFunctionImage(gScore_);

        return gScore_;
    }

    // 带方向感知的距离场计算
    // 额外输入一个梯度场，计算每个点到光源的距离，并根据梯度场的方向调整距离
    const Field<double>& Solver::computeDistanceFunction(const Field<vector3f>& gradientField) {
        reset();
        auto startTime = std::chrono::high_resolution_clock::now();

        for (size_t i = 0; i < nx_; ++i) {
            for (size_t j = 0; j < ny_; ++j)
                for (size_t k = 0; k < nz_; ++k) {
                    if (sharedVisibilityField_->get(i, j, k) <= visibilityThreshold_) {
                        initial_frontline.push_back(i);
                        initial_frontline.push_back(j);
                        initial_frontline.push_back(k);
                    }
                }
        }
        // Init
        int x = 0, y = 0, z = 0, neighbour_x = 0, neighbour_y = 0, neighbour_z = 0;
        double g = 0;
        // Fill in data from initial frontline
        for (size_t i = 0; i < initial_frontline.size(); i += 3) {
            x = initial_frontline[i];
            y = initial_frontline[i + 1];
            z = initial_frontline[i + 2];
            g = 0;
            openSet_->push(Node{ x, y,z, g });

            gScore_(x, y, z) = g;
            updated_(x, y, z) = true;
            cameFrom_(x, y, z) = nb_of_sources_;
            lightSources_[nb_of_sources_] = { x, y,z };
            const auto key = hashFunction(x, y, z, nb_of_sources_);
            visibilityHashMap_[key] = lightStrength_;
            ++nb_of_sources_;
        }

        // For queing unique sources from neighbours of neighbour
        std::vector<size_t> potentialSources;
        potentialSources.reserve(10);
        std::vector<std::pair<double, size_t>> potentialDistances;
        potentialDistances.reserve(10);

        double distance = 0;

        while (openSet_->size() > 0) {
            auto& current = openSet_->top();
            x = current.x;
            y = current.y;
            z = current.z;

            openSet_->pop();

            // Expand frontline at current & update neighbours
            for (size_t j = 0; j < 78; j += 3) {
                // NOTE those are always positive
                neighbour_x = x + neighbours_[j];
                neighbour_y = y + neighbours_[j + 1];
                neighbour_z = z + neighbours_[j + 2];

                // Box check
                if (!isValid(neighbour_x, neighbour_y, neighbour_z)) {
                    continue;
                }
                if (updated_(neighbour_x, neighbour_y, neighbour_z)) {
                    continue;
                };

                queuePotentialSources(potentialSources, neighbour_x, neighbour_y, neighbour_z);
                potentialDistances.clear();
                for (int k = 0; k < potentialSources.size(); ++k) {
                    int potentialSource = potentialSources[k];
                    int LS_x = lightSources_[potentialSource].x;
                    int LS_y = lightSources_[potentialSource].y;
                    int LS_z = lightSources_[potentialSource].z;

                    distance = gScore_(LS_x, LS_y, LS_z) +
                        evaluateDistance(LS_x, LS_y, LS_z, neighbour_x, neighbour_y, neighbour_z, gradientField(neighbour_x, neighbour_y, neighbour_z));
                    potentialDistances.push_back(
                        std::pair<double, int>{distance, potentialSource});
                }

                auto minimum_element =
                    std::min_element(potentialDistances.begin(), potentialDistances.end(),
                        [](const auto& lhs, const auto& rhs) {
                            return lhs.first < rhs.first;
                        });
                distance = minimum_element->first;
                // use source giving least distance
                gScore_(neighbour_x, neighbour_y, neighbour_z) = distance;
                cameFrom_(neighbour_x, neighbour_y, neighbour_z) = minimum_element->second;
                openSet_->push(Node{ neighbour_x, neighbour_y,neighbour_z, distance });
                updated_(neighbour_x, neighbour_y, neighbour_z) = true;
            }
        };

        auto stopTime = std::chrono::high_resolution_clock::now();
        auto executionDuration = durationInMicroseconds(startTime, stopTime);

        if (!Outsilent) {
            std::cout << std::endl << "################# VBM:Distance function computation "
                "output #####################"
                << std::endl;

            if (!Outdetails) {
                std::cout << "Constructed distance function" << std::endl;
                std::cout << "Execution time in s: " << (double)executionDuration / 1000.0 / CLOCKS_PER_SEC << "s"
                    << std::endl;
            }
            std::cout << "####################################################"
                "############################"
                << std::endl << std::endl;
        }

        //saveDistanceFunctionImage(gScore_);

        return gScore_;
    }

    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    inline void Solver::queuePotentialSources(std::vector<size_t>& potentialSources,
        const int neighbour_x,
        const int neighbour_y,
        const int neighbour_z) const {
        size_t potentialSource_x = 0, potentialSource_y = 0, potentialSource_z = 0, lightSource_num = 0;
        potentialSources.clear();
        // Queue sources from updated neighbours of neighbour
        for (size_t k = 0; k < 78; k += 3) {
            // NOTE those are always positive
            potentialSource_x = neighbour_x + neighbours_[k];
            potentialSource_y = neighbour_y + neighbours_[k + 1];
            potentialSource_z = neighbour_z + neighbours_[k + 2];

            // Box check
            if (!isValid(potentialSource_x, potentialSource_y, potentialSource_z)) {
                continue;
            };
            if (!updated_(potentialSource_x, potentialSource_y, potentialSource_z)) {
                continue;
            };

            lightSource_num = cameFrom_(potentialSource_x, potentialSource_y, potentialSource_z);
            // Pick only unique sources (no repitition in potentialSources)
            if (std::find(potentialSources.begin(), potentialSources.end(),
                lightSource_num) == potentialSources.end()) {
                potentialSources.push_back(lightSource_num);
            }
        }
    }

    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    void Solver::getPotentialDistances(
        const std::vector<size_t>& potentialSources,
        std::vector<std::pair<double, size_t>>& potentialDistances,
        const int neighbour_x, const int neighbour_y, const int neighbour_z) {
        size_t LS_x = 0, LS_y = 0, LS_z = 0, potentialSource = 0;
        double distance = 0;
        potentialDistances.clear();
        for (size_t k = 0; k < potentialSources.size(); ++k) {
            potentialSource = potentialSources[k];
            LS_x = lightSources_[potentialSource].x;
            LS_y = lightSources_[potentialSource].y;
            LS_z = lightSources_[potentialSource].z;

            // update visibility from source
            updatePointVisibility(potentialSource, LS_x, LS_y, LS_z, neighbour_x,
                neighbour_y, neighbour_z);
            distance = std::numeric_limits<double>::infinity();
            const auto key = hashFunction(neighbour_x, neighbour_y, neighbour_z, potentialSource);
            if (visibilityHashMap_.at(key) >= visibilityThreshold_) {
                distance = gScore_(LS_x, LS_y, LS_z) +
                    evaluateDistance(LS_x, LS_y, LS_z, neighbour_x, neighbour_y, neighbour_z);
            }
            potentialDistances.push_back(
                std::pair<double, int>{distance, potentialSource});
        }
    }


    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    void Solver::getPotentialDistancesField(
        const std::vector<size_t>& potentialSources,
        std::vector<std::pair<double, size_t>>& potentialDistances,
        const int neighbour_x, const int neighbour_y, const int neighbour_z) {
        size_t LS_x = 0, LS_y = 0, LS_z = 0, potentialSource = 0;
        double distance = 0;
        potentialDistances.clear();
        for (size_t k = 0; k < potentialSources.size(); ++k) {
            potentialSource = potentialSources[k];
            LS_x = lightSources_[potentialSource].x;
            LS_y = lightSources_[potentialSource].y;
            LS_z = lightSources_[potentialSource].z;

            // update visibility from source
            updatePointVisibility(potentialSource, LS_x, LS_y, LS_z,
                neighbour_x, neighbour_y, neighbour_z);
            distance = std::numeric_limits<double>::infinity();
            const auto key = hashFunction(neighbour_x, neighbour_y, neighbour_z, potentialSource);
            if (visibilityHashMap_.at(key) >= visibilityThreshold_) {
                distance =
                    gScore_(LS_x, LS_y, LS_z) +
                    evaluateDistance(LS_x, LS_y, LS_z, neighbour_x, neighbour_y, neighbour_z);
            }
            potentialDistances.push_back(
                std::pair<double, int>{distance, potentialSource});
        }
    }


    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    void Solver::getPotentialDistancesSpeedField(
        const std::vector<size_t>& potentialSources,
        std::vector<std::pair<double, size_t>>& potentialDistances,
        const int neighbour_x, const int neighbour_y, const int neighbour_z) {
        size_t LS_x = 0, LS_y = 0, LS_z = 0, potentialSource = 0;
        double distance = 0;
        potentialDistances.clear();
        for (size_t k = 0; k < potentialSources.size(); ++k) {
            potentialSource = potentialSources[k];
            LS_x = lightSources_[potentialSource].x;
            LS_y = lightSources_[potentialSource].y;
            LS_z = lightSources_[potentialSource].z;

            // update visibility from source
            updatePointVisibility(potentialSource, LS_x, LS_y, LS_z,
                neighbour_x, neighbour_y, neighbour_z);
            distance = std::numeric_limits<double>::infinity();
            const auto key = hashFunction(neighbour_x, neighbour_y, neighbour_z, potentialSource);
            if (visibilityHashMap_.at(key) >= visibilityThreshold_) {
                distance =
                    gScore_(LS_x, LS_y, LS_z) +
                    evaluateDistanceSpeedField(LS_x, LS_y, LS_z, neighbour_x, neighbour_y, neighbour_z);
            }
            potentialDistances.push_back(
                std::pair<double, int>{distance, potentialSource});
        }
    }

    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    void Solver::createNewPivot(const int x, const int y, const int z,
        const int neighbour_x, const int neighbour_y, const int neighbour_z) {
        int pivot_neighbour_x, pivot_neighbour_y, pivot_neighbour_z;

        // Pushback parent as a new lightSource
        lightSources_[nb_of_sources_] = { x, y,z }; // {x, y};
        // Pusback pivot & update light source visibility
        const auto key = hashFunction(x, y, z, nb_of_sources_);
        visibilityHashMap_[key] = lightStrength_;
        // Update maps of new pivot_
        // Update neighbours of initial frontline points - both distance & visibility
        for (size_t p = 0; p < 78; p += 3) {
            // NOTE those are always positive
            pivot_neighbour_x = x + neighbours_[p];
            pivot_neighbour_y = y + neighbours_[p + 1];
            pivot_neighbour_z = z + neighbours_[p + 2];

            // Box check
            if (!isValid(pivot_neighbour_x, pivot_neighbour_y, pivot_neighbour_z)) {
                continue;
            }
            // Update neighbour visibility
            updatePointVisibility(nb_of_sources_, x, y, z,
                pivot_neighbour_x, pivot_neighbour_y, pivot_neighbour_z);
        }
        cameFrom_(neighbour_x, neighbour_y, neighbour_z) = nb_of_sources_;
        gScore_(neighbour_x, neighbour_y, neighbour_z) =
            gScore_(x, y, z) + evaluateDistance(x, y, z, neighbour_x, neighbour_y, neighbour_z);
        ++nb_of_sources_;
    }

    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    void Solver::updatePointVisibility(const size_t lightSourceNumber,
        const int LS_x, const int LS_y, const int LS_z,
        const int x, const int y, const int z) {
        // Variable initialization
        double v = 0;
        double c = 0;
        double c1 = 0, c2 = 0;

        // check if out of range, return
        // if (!isValid(x, y)) {
        //   visibilityHashMap_[hashFunction(x, y, lightSourceNumber)] = 0;
        //   return;
        // }

        // Check if visibility value already exists
        auto key = hashFunction(x, y, z, lightSourceNumber);
        if (visibilityHashMap_.count(key)) {
            return;
        }
        if (sharedVisibilityField_->get(x, y, z) < visibilityThreshold_) {
            visibilityHashMap_[key] = 0;
            return;
        }

        if (x == LS_x && y == LS_y) {
            if (z - LS_z > 0) {
                key = hashFunction(x, y, z - 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                }
                v = visibilityHashMap_.at(key);
            }
            else {
                key = hashFunction(x, y, z + 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                }
                v = visibilityHashMap_.at(key);
            }
        }
        else if (y == LS_y && z == LS_z) {
            if (x - LS_x > 0) {
                key = hashFunction(x - 1, y, z, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                }
                v = visibilityHashMap_.at(key);
            }
            else {
                key = hashFunction(x + 1, y, z, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                }
                v = visibilityHashMap_.at(key);
            }
        }
        else if (x == LS_x && z == LS_z) {
            if (y - LS_y > 0) {
                key = hashFunction(x, y - 1, z, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                }
                v = visibilityHashMap_.at(key);
            }
            else {
                key = hashFunction(x, y + 1, z, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                }
                v = visibilityHashMap_.at(key);
            }
        }
        else if (z == LS_z) {
            // Quadrant 1
            if (x - LS_x > 0 && y - LS_y > 0) {
                key = hashFunction(x - 1, y - 1, z, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y - 1, z);
                }
                if (x - LS_x == y - LS_y) {
                    v = visibilityHashMap_.at(key);
                }
                else if (x - LS_x < y - LS_y) {
                    const auto key_1 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    c = static_cast<double>(x - LS_x) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x > y - LS_y) {
                    const auto key_2 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    c = static_cast<double>(y - LS_y) / (x - LS_x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Quadrant 2
            else if (x - LS_x < 0 && y - LS_y > 0) {
                key = hashFunction(x + 1, y - 1, z, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y - 1, z);
                }
                if (LS_x - x == y - LS_y) {
                    v = visibilityHashMap_.at(key);
                }
                else if (LS_x - x < y - LS_y) {
                    const auto key_1 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    c = static_cast<double>(LS_x - x) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x > y - LS_y) {
                    const auto key_2 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    c = static_cast<double>(y - LS_y) / (LS_x - x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Quadrant 3
            else if (x - LS_x < 0 && y - LS_y < 0) {
                key = hashFunction(x + 1, y + 1, z, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y + 1, z);
                }
                if (LS_x - x == LS_y - y) {
                    v = visibilityHashMap_.at(key);
                }
                else if (LS_x - x < LS_y - y) {
                    const auto key_1 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    c = static_cast<double>(LS_x - x) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x > LS_y - y) {
                    const auto key_2 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    c = static_cast<double>(LS_y - y) / (LS_x - x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Quadrant 4
            else if (x - LS_x > 0 && y - LS_y < 0) {
                key = hashFunction(x - 1, y + 1, z, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y + 1, z);
                }
                if (x - LS_x == LS_y - y) {
                    v = visibilityHashMap_.at(key);
                }
                else if (x - LS_x < LS_y - y) {
                    const auto key_1 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    c = static_cast<double>(x - LS_x) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x > LS_y - y) {
                    const auto key_2 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    c = static_cast<double>(LS_y - y) / (x - LS_x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
            }
        }
        else if (y == LS_y) {
            // Quadrant 1
            if (x - LS_x > 0 && z - LS_z > 0) {
                key = hashFunction(x - 1, y, z - 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z - 1);
                }
                if (x - LS_x == z - LS_z) {
                    v = visibilityHashMap_.at(key);
                }
                else if (x - LS_x < z - LS_z) {
                    const auto key_1 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    c = static_cast<double>(x - LS_x) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x > z - LS_z) {
                    const auto key_2 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    c = static_cast<double>(z - LS_z) / (x - LS_x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Quadrant 2
            else if (x - LS_x < 0 && z - LS_z > 0) {
                key = hashFunction(x + 1, y, z - 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z - 1);
                }
                if (LS_x - x == z - LS_z) {
                    v = visibilityHashMap_.at(key);
                }
                else if (LS_x - x < z - LS_z) {
                    const auto key_1 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    c = static_cast<double>(LS_x - x) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x > z - LS_z) {
                    const auto key_2 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    c = static_cast<double>(z - LS_z) / (LS_x - x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Quadrant 3
            else if (x - LS_x < 0 && z - LS_z < 0) {
                key = hashFunction(x + 1, y, z + 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z + 1);
                }
                if (LS_x - x == LS_z - z) {
                    v = visibilityHashMap_.at(key);
                }
                else if (LS_x - x < LS_z - z) {
                    const auto key_1 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    c = static_cast<double>(LS_x - x) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x > LS_z - z) {
                    const auto key_2 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    c = static_cast<double>(LS_z - z) / (LS_x - x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Quadrant 4
            else if (x - LS_x > 0 && z - LS_z < 0) {
                key = hashFunction(x - 1, y, z + 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z + 1);
                }
                if (x - LS_x == LS_z - z) {
                    v = visibilityHashMap_.at(key);
                }
                else if (x - LS_x < LS_z - z) {
                    const auto key_1 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    c = static_cast<double>(x - LS_x) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x > LS_z - z) {
                    const auto key_2 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    c = static_cast<double>(LS_z - z) / (x - LS_x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
            }
        }
        else if (x == LS_x) {
            // Quadrant 1
            if (y - LS_y > 0 && z - LS_z > 0) {
                key = hashFunction(x, y - 1, z - 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z - 1);
                }
                if (y - LS_y == z - LS_z) {
                    v = visibilityHashMap_.at(key);
                }
                else if (y - LS_y < z - LS_z) {
                    const auto key_1 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    c = static_cast<double>(y - LS_y) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y > z - LS_z) {
                    const auto key_2 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    c = static_cast<double>(z - LS_z) / (y - LS_y);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Quadrant 2
            else if (y - LS_y < 0 && z - LS_z > 0) {
                key = hashFunction(x, y + 1, z - 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z - 1);
                }
                if (LS_y - y == z - LS_z) {
                    v = visibilityHashMap_.at(key);
                }
                else if (LS_y - y < z - LS_z) {
                    const auto key_1 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    c = static_cast<double>(LS_y - y) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y > z - LS_z) {
                    const auto key_2 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    c = static_cast<double>(z - LS_z) / (LS_y - y);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Quadrant 3
            else if (y - LS_y < 0 && z - LS_z < 0) {
                key = hashFunction(x, y + 1, z + 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z + 1);
                }
                if (LS_y - y == LS_z - z) {
                    v = visibilityHashMap_.at(key);
                }
                else if (LS_y - y < LS_z - z) {
                    const auto key_1 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    c = static_cast<double>(LS_y - y) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y > LS_z - z) {
                    const auto key_2 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    c = static_cast<double>(LS_z - z) / (LS_y - y);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Quadrant 4
            else if (y - LS_y > 0 && z - LS_z < 0) {
                key = hashFunction(x, y - 1, z + 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z + 1);
                }
                if (y - LS_y == LS_z - z) {
                    v = visibilityHashMap_.at(key);
                }
                else if (y - LS_y < LS_z - z) {
                    const auto key_1 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    c = static_cast<double>(y - LS_y) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y > LS_z - z) {
                    const auto key_2 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    c = static_cast<double>(LS_z - z) / (y - LS_y);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
            }
        }
        else {
            // Q1
            if ((x - LS_x > 0) && (y - LS_y > 0) && (z - LS_z > 0)) {
                key = hashFunction(x - 1, y - 1, z - 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y - 1, z - 1);
                }
                if (x - LS_x == y - LS_y && z - LS_z == x - LS_x) {
                    v = visibilityHashMap_.at(key);
                }
                else if (x - LS_x == y - LS_y && z - LS_z < x - LS_x)
                {
                    const auto key_1 = hashFunction(x - 1, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y - 1, z);
                    }
                    c = static_cast<double>(z - LS_z) / (x - LS_x);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x == y - LS_y && z - LS_z > x - LS_x)
                {
                    const auto key_2 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    c = static_cast<double>(x - LS_x) / (z - LS_z);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x == z - LS_z && y - LS_y < z - LS_z)
                {
                    const auto key_1 = hashFunction(x - 1, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z - 1);
                    }
                    c = static_cast<double>(y - LS_y) / (x - LS_x);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x == z - LS_z && y - LS_y > z - LS_z)
                {
                    const auto key_2 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    c = static_cast<double>(x - LS_x) / (y - LS_y);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y == z - LS_z && x - LS_x < y - LS_y)
                {
                    const auto key_1 = hashFunction(x, y - 1, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z - 1);
                    }
                    c = static_cast<double>(x - LS_x) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y == z - LS_z && x - LS_x > y - LS_y)
                {
                    const auto key_2 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    c = static_cast<double>(y - LS_y) / (x - LS_x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x < y - LS_y && y - LS_y < z - LS_z)
                {
                    const auto key_1 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    const auto key_2 = hashFunction(x, y - 1, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z - 1);
                    }
                    c1 = static_cast<double> (x - LS_x + y - LS_y) / 2 * (z - LS_z);
                    c2 = static_cast<double> (x - LS_x) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x < z - LS_z && z - LS_z < y - LS_y) {
                    const auto key_1 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    const auto key_2 = hashFunction(x, y - 1, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z - 1);
                    }
                    c1 = static_cast<double> (x - LS_x + z - LS_z) / 2 * (y - LS_y);
                    c2 = static_cast<double> (x - LS_x) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y < z - LS_z && z - LS_z < x - LS_x)
                {
                    const auto key_1 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    const auto key_2 = hashFunction(x - 1, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z - 1);
                    }
                    c1 = static_cast<double> (y - LS_y + z - LS_z) / 2 * (x - LS_x);
                    c2 = static_cast<double> (y - LS_y) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y < x - LS_x && x - LS_x < z - LS_z) {
                    const auto key_1 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    const auto key_2 = hashFunction(x - 1, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z - 1);
                    }
                    c1 = static_cast<double> (y - LS_y + x - LS_x) / 2 * (z - LS_z);
                    c2 = static_cast<double> (y - LS_y) / (x - LS_x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (z - LS_z < x - LS_x && x - LS_x < y - LS_y) {
                    const auto key_1 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    const auto key_2 = hashFunction(x - 1, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y - 1, z);
                    }
                    c1 = static_cast<double> (z - LS_z + x - LS_x) / 2 * (y - LS_y);
                    c2 = static_cast<double> (z - LS_z) / (x - LS_x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (z - LS_z < y - LS_y && y - LS_y < x - LS_x) {
                    const auto key_1 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    const auto key_2 = hashFunction(x - 1, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y - 1, z);
                    }
                    c1 = static_cast<double> (z - LS_z + y - LS_y) / 2 * (x - LS_x);
                    c2 = static_cast<double> (z - LS_z) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Q2
            else if ((x - LS_x < 0) && (y - LS_y > 0) && (z - LS_z > 0)) {
                key = hashFunction(x + 1, y - 1, z - 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y - 1, z - 1);
                }
                if (LS_x - x == y - LS_y && z - LS_z == LS_x - x) {
                    v = visibilityHashMap_.at(key);
                }
                else if (LS_x - x == y - LS_y && z - LS_z < LS_x - x)
                {
                    const auto key_1 = hashFunction(x + 1, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y - 1, z);
                    }
                    c = static_cast<double>(z - LS_z) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x == y - LS_y && z - LS_z > LS_x - x)
                {
                    const auto key_2 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    c = static_cast<double>(LS_x - x) / (z - LS_z);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x == z - LS_z && y - LS_y < z - LS_z)
                {
                    const auto key_1 = hashFunction(x + 1, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z - 1);
                    }
                    c = static_cast<double>(y - LS_y) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x == z - LS_z && y - LS_y > z - LS_z)
                {
                    const auto key_2 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    c = static_cast<double>(z - LS_z) / (y - LS_y);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y == z - LS_z && LS_x - x < y - LS_y)
                {
                    const auto key_1 = hashFunction(x, y - 1, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z - 1);
                    }
                    c = static_cast<double>(LS_x - x) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y == z - LS_z && LS_x - x > y - LS_y)
                {
                    const auto key_2 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    c = static_cast<double>(y - LS_y) / (LS_x - x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x < y - LS_y && y - LS_y < z - LS_z) {
                    const auto key_1 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    const auto key_2 = hashFunction(x, y - 1, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z - 1);
                    }
                    c1 = static_cast<double> (LS_x - x + y - LS_y) / 2 * (z - LS_z);
                    c2 = static_cast<double> (LS_x - x) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x < z - LS_z && z - LS_z < y - LS_y)
                {
                    const auto key_1 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    const auto key_2 = hashFunction(x, y - 1, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z - 1);
                    }
                    c1 = static_cast<double> (LS_x - x + z - LS_z) / 2 * (y - LS_y);
                    c2 = static_cast<double> (LS_x - x) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y < z - LS_z && z - LS_z < LS_x - x)
                {
                    const auto key_1 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    const auto key_2 = hashFunction(x + 1, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z - 1);
                    }
                    c1 = static_cast<double> (y - LS_y + z - LS_z) / 2 * (LS_x - x);
                    c2 = static_cast<double> (y - LS_y) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y < LS_x - x && LS_x - x < z - LS_z)
                {
                    const auto key_1 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    const auto key_2 = hashFunction(x + 1, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z - 1);
                    }
                    c1 = static_cast<double> (LS_x - x + y - LS_y) / 2 * (z - LS_z);
                    c2 = static_cast<double> (LS_x - x) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (z - LS_z < y - LS_y && y - LS_y < LS_x - x) {
                    const auto key_1 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    const auto key_2 = hashFunction(x + 1, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y - 1, z);
                    }
                    c1 = static_cast<double> (z - LS_z + y - LS_y) / 2 * (LS_x - x);
                    c2 = static_cast<double> (z - LS_z) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (z - LS_z < LS_x - x && LS_x - x < y - LS_y)
                {
                    const auto key_1 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    const auto key_2 = hashFunction(x + 1, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y - 1, z);
                    }
                    c1 = static_cast<double> (z - LS_z + LS_x - x) / 2 * (y - LS_y);
                    c2 = static_cast<double> (z - LS_z) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Q3
            else if ((x - LS_x < 0) && (y - LS_y < 0) && (z - LS_z > 0)) {
                key = hashFunction(x + 1, y + 1, z - 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y + 1, z - 1);
                }
                if (LS_x - x == LS_y - y && z - LS_z == LS_x - x) {
                    v = visibilityHashMap_.at(key);
                }
                else if (LS_x - x == LS_y - y && z - LS_z < LS_x - x)
                {
                    const auto key_1 = hashFunction(x + 1, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y + 1, z);
                    }
                    c = static_cast<double>(z - LS_z) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x == LS_y - y && z - LS_z > LS_x - x)
                {
                    const auto key_2 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    c = static_cast<double>(LS_x - x) / (z - LS_z);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x == z - LS_z && LS_y - y < z - LS_z)
                {
                    const auto key_1 = hashFunction(x + 1, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z - 1);
                    }
                    c = static_cast<double>(LS_y - y) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x == z - LS_z && LS_y - y > z - LS_z)
                {
                    const auto key_2 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    c = static_cast<double>(z - LS_z) / (LS_y - y);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y == z - LS_z && LS_x - x < LS_y - y)
                {
                    const auto key_1 = hashFunction(x, y + 1, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z - 1);
                    }
                    c = static_cast<double>(LS_x - x) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y == z - LS_z && LS_x - x > LS_y - y)
                {
                    const auto key_2 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    c = static_cast<double>(LS_y - y) / (LS_x - x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x < LS_y - y && y - LS_y < z - LS_z)
                {
                    const auto key_1 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    const auto key_2 = hashFunction(x, y + 1, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z - 1);
                    }
                    c1 = static_cast<double>(LS_x - x + LS_y - y) / 2 * (z - LS_z);
                    c2 = static_cast<double>(LS_x - x) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x < z - LS_z && z - LS_z < LS_y - y)
                {
                    const auto key_1 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    const auto key_2 = hashFunction(x, y + 1, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z - 1);
                    }
                    c1 = static_cast<double>(z - LS_z + LS_x - x) / 2 * (LS_y - y);
                    c2 = static_cast<double>(z - LS_z) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y < z - LS_z && z - LS_z < LS_x - x)
                {
                    const auto key_1 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    const auto key_2 = hashFunction(x + 1, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z - 1);
                    }
                    c1 = static_cast<double>(LS_y - y + z - LS_z) / 2 * (LS_x - x);
                    c2 = static_cast<double>(LS_y - y) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y < LS_x - x && LS_x - x < z - LS_z)
                {
                    const auto key_1 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    const auto key_2 = hashFunction(x + 1, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z - 1);
                    }
                    c1 = static_cast<double>(LS_y - y + LS_x - x) / 2 * (z - LS_z);
                    c2 = static_cast<double>(LS_y - y) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (z - LS_z < LS_y - y && LS_y - y < LS_x - x)
                {
                    const auto key_1 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    const auto key_2 = hashFunction(x + 1, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y + 1, z);
                    }
                    c1 = static_cast<double>(z - LS_z + LS_y - y) / 2 * (LS_x - x);
                    c2 = static_cast<double>(z - LS_z) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (z - LS_z < LS_x - x && LS_x - x < LS_y - y)
                {
                    const auto key_1 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    const auto key_2 = hashFunction(x + 1, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y + 1, z);
                    }
                    c1 = static_cast<double>(z - LS_z + LS_x - x) / 2 * (LS_y - y);
                    c2 = static_cast<double>(z - LS_z) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Q4
            else if ((x - LS_x > 0) && (y - LS_y < 0) && (z - LS_z > 0)) {
                const auto key = hashFunction(x - 1, y + 1, z - 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y + 1, z - 1);
                }
                if (x - LS_x == LS_y - y && LS_y - y == z - LS_z) {
                    v = visibilityHashMap_.at(key);
                }
                else if (x - LS_x == LS_y - y && z - LS_z < LS_y - y)
                {
                    const auto key_1 = hashFunction(x - 1, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y + 1, z);
                    }
                    c = static_cast<double>(z - LS_z) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x == LS_y - y && z - LS_z > LS_y - y)
                {
                    const auto key_2 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    c = static_cast<double>(LS_y - y) / (z - LS_z);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x == z - LS_z && LS_y - y < z - LS_z)
                {
                    const auto key_1 = hashFunction(x - 1, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z - 1);
                    }
                    c = static_cast<double>(LS_y - y) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x == z - LS_z && LS_y - y > z - LS_z)
                {
                    const auto key_2 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    c = static_cast<double>(z - LS_z) / (LS_y - y);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y == z - LS_z && x - LS_x < LS_y - y)
                {
                    const auto key_1 = hashFunction(x, y + 1, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z - 1);
                    }
                    c = static_cast<double>(x - LS_x) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y == z - LS_z && x - LS_x > LS_y - y)
                {
                    const auto key_2 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    c = static_cast<double>(LS_y - y) / (x - LS_x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x < LS_y - y && LS_y - y < z - LS_z)
                {
                    const auto key_1 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    const auto key_2 = hashFunction(x, y + 1, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z - 1);
                    }
                    c1 = static_cast<double>(x - LS_x + LS_y - y) / 2 * (z - LS_z);
                    c2 = static_cast<double>(x - LS_x) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x < z - LS_z && z - LS_z < LS_y - y)
                {
                    const auto key_1 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    const auto key_2 = hashFunction(x, y + 1, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z - 1);
                    }
                    c1 = static_cast<double>(x - LS_x + z - LS_z) / 2 * (LS_y - y);
                    c2 = static_cast<double>(x - LS_x) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y < z - LS_z && z - LS_z < x - LS_x)
                {
                    const auto key_1 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    const auto key_2 = hashFunction(x - 1, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z - 1);
                    }
                    c1 = static_cast<double>(LS_y - y + z - LS_z) / 2 * (x - LS_x);
                    c2 = static_cast<double>(LS_y - y) / (z - LS_z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y < x - LS_x && x - LS_x < z - LS_z)
                {
                    const auto key_1 = hashFunction(x, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z - 1);
                    }
                    const auto key_2 = hashFunction(x - 1, y, z - 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z - 1);
                    }
                    c1 = static_cast<double>(LS_y - y + x - LS_x) / 2 * (z - LS_z);
                    c2 = static_cast<double>(LS_y - y) / (x - LS_x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (z - LS_z < x - LS_x && x - LS_x < LS_y - y)
                {
                    const auto key_1 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    const auto key_2 = hashFunction(x - 1, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y + 1, z);
                    }
                    c1 = static_cast<double>(z - LS_z + x - LS_x) / 2 * (LS_y - y);
                    c2 = static_cast<double>(z - LS_z) / (x - LS_x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (z - LS_z < LS_y - y && LS_y - y < x - LS_x)
                {
                    const auto key_1 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    const auto key_2 = hashFunction(x - 1, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y + 1, z);
                    }
                    c1 = static_cast<double>(z - LS_z + LS_y - y) / 2 * (x - LS_x);
                    c2 = static_cast<double>(z - LS_z) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Q5
            else if ((x - LS_x > 0) && (y - LS_y > 0) && (z - LS_z < 0)) {
                key = hashFunction(x - 1, y - 1, z + 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y - 1, z + 1);
                }
                if (x - LS_x == y - LS_y && LS_z - z == x - LS_x) {
                    v = visibilityHashMap_.at(key);
                }
                else if (x - LS_x == y - LS_y && LS_z - z < x - LS_x)
                {
                    const auto key_1 = hashFunction(x - 1, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y - 1, z);
                    }
                    c = static_cast<double>(LS_z - z) / (x - LS_x);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x == y - LS_y && LS_z - z > x - LS_x)
                {
                    const auto key_2 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    c = static_cast<double>(x - LS_x) / (LS_z - z);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x == LS_z - z && y - LS_y < LS_z - z)
                {
                    const auto key_1 = hashFunction(x - 1, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z + 1);
                    }
                    c = static_cast<double>(y - LS_y) / (x - LS_x);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x == LS_z - z && y - LS_y > LS_z - z)
                {
                    const auto key_2 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    c = static_cast<double>(x - LS_x) / (y - LS_y);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y == LS_z - z && x - LS_x < y - LS_y)
                {
                    const auto key_1 = hashFunction(x, y - 1, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z + 1);
                    }
                    c = static_cast<double>(x - LS_x) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y == LS_z - z && x - LS_x > y - LS_y)
                {
                    const auto key_2 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    c = static_cast<double>(y - LS_y) / (x - LS_x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x < y - LS_y && y - LS_y < LS_z - z)
                {
                    const auto key_1 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    const auto key_2 = hashFunction(x, y - 1, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z + 1);
                    }
                    c1 = static_cast<double> (x - LS_x + y - LS_y) / 2 * (LS_z - z);
                    c2 = static_cast<double> (x - LS_x) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x < LS_z - z && LS_z - z < y - LS_y) {
                    const auto key_1 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    const auto key_2 = hashFunction(x, y - 1, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z + 1);
                    }
                    c1 = static_cast<double> (x - LS_x + LS_z - z) / 2 * (y - LS_y);
                    c2 = static_cast<double> (x - LS_x) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y < LS_z - z && LS_z - z < x - LS_x)
                {
                    const auto key_1 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    const auto key_2 = hashFunction(x - 1, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z + 1);
                    }
                    c1 = static_cast<double> (y - LS_y + LS_z - z) / 2 * (x - LS_x);
                    c2 = static_cast<double> (y - LS_y) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y < x - LS_x && x - LS_x < LS_z - z) {
                    const auto key_1 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    const auto key_2 = hashFunction(x - 1, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z + 1);
                    }
                    c1 = static_cast<double> (y - LS_y + x - LS_x) / 2 * (LS_z - z);
                    c2 = static_cast<double> (y - LS_y) / (x - LS_x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_z - z < x - LS_x && x - LS_x < y - LS_y) {
                    const auto key_1 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    const auto key_2 = hashFunction(x - 1, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y - 1, z);
                    }
                    c1 = static_cast<double> (LS_z - z + x - LS_x) / 2 * (y - LS_y);
                    c2 = static_cast<double> (LS_z - z) / (x - LS_x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_z - z < y - LS_y && y - LS_y < x - LS_x) {
                    const auto key_1 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    const auto key_2 = hashFunction(x - 1, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y - 1, z);
                    }
                    c1 = static_cast<double> (LS_z - z + y - LS_y) / 2 * (x - LS_x);
                    c2 = static_cast<double> (LS_z - z) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Q6
            else if ((x - LS_x < 0) && (y - LS_y > 0) && (z - LS_z < 0)) {
                key = hashFunction(x + 1, y - 1, z + 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y - 1, z + 1);
                }
                if (LS_x - x == y - LS_y && LS_z - z == LS_x - x) {
                    v = visibilityHashMap_.at(key);
                }
                else if (LS_x - x == y - LS_y && LS_z - z < LS_x - x)
                {
                    const auto key_1 = hashFunction(x + 1, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y - 1, z);
                    }
                    c = static_cast<double>(LS_z - z) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x == y - LS_y && LS_z - z > LS_x - x)
                {
                    const auto key_2 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    c = static_cast<double>(LS_x - x) / (LS_z - z);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x == LS_z - z && y - LS_y < LS_z - z)
                {
                    const auto key_1 = hashFunction(x + 1, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z + 1);
                    }
                    c = static_cast<double>(y - LS_y) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x == LS_z - z && y - LS_y > LS_z - z)
                {
                    const auto key_2 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    c = static_cast<double>(LS_z - z) / (y - LS_y);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y == LS_z - z && LS_x - x < y - LS_y)
                {
                    const auto key_1 = hashFunction(x, y - 1, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z + 1);
                    }
                    c = static_cast<double>(LS_x - x) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y == LS_z - z && LS_x - x > y - LS_y)
                {
                    const auto key_2 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    c = static_cast<double>(y - LS_y) / (LS_x - x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x < y - LS_y && y - LS_y < LS_z - z) {
                    const auto key_1 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    const auto key_2 = hashFunction(x, y - 1, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z + 1);
                    }
                    c1 = static_cast<double> (LS_x - x + y - LS_y) / 2 * (LS_z - z);
                    c2 = static_cast<double> (LS_x - x) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x < LS_z - z && LS_z - z < y - LS_y)
                {
                    const auto key_1 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    const auto key_2 = hashFunction(x, y - 1, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z + 1);
                    }
                    c1 = static_cast<double> (LS_x - x + LS_z - z) / 2 * (y - LS_y);
                    c2 = static_cast<double> (LS_x - x) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y < LS_z - z && LS_z - z < LS_x - x)
                {
                    const auto key_1 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    const auto key_2 = hashFunction(x + 1, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z + 1);
                    }
                    c1 = static_cast<double> (y - LS_y + LS_z - z) / 2 * (LS_x - x);
                    c2 = static_cast<double> (y - LS_y) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (y - LS_y < LS_x - x && LS_x - x < LS_z - z)
                {
                    const auto key_1 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    const auto key_2 = hashFunction(x + 1, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z + 1);
                    }
                    c1 = static_cast<double> (LS_x - x + y - LS_y) / 2 * (LS_z - z);
                    c2 = static_cast<double> (LS_x - x) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_z - z < y - LS_y && y - LS_y < LS_x - x) {
                    const auto key_1 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    const auto key_2 = hashFunction(x + 1, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y - 1, z);
                    }
                    c1 = static_cast<double> (LS_z - z + y - LS_y) / 2 * (LS_x - x);
                    c2 = static_cast<double> (LS_z - z) / (y - LS_y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_z - z < LS_x - x && LS_x - x < y - LS_y)
                {
                    const auto key_1 = hashFunction(x, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y - 1, z);
                    }
                    const auto key_2 = hashFunction(x + 1, y - 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y - 1, z);
                    }
                    c1 = static_cast<double> (LS_z - z + LS_x - x) / 2 * (y - LS_y);
                    c2 = static_cast<double> (LS_z - z) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Q7
            else if ((x - LS_x < 0) && (y - LS_y < 0) && (z - LS_z < 0)) {
                key = hashFunction(x + 1, y + 1, z + 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y + 1, z + 1);
                }
                if (LS_x - x == LS_y - y && LS_z - z == LS_x - x) {
                    v = visibilityHashMap_.at(key);
                }
                else if (LS_x - x == LS_y - y && LS_z - z < LS_x - x)
                {
                    const auto key_1 = hashFunction(x + 1, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y + 1, z);
                    }
                    c = static_cast<double>(LS_z - z) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x == LS_y - y && LS_z - z > LS_x - x)
                {
                    const auto key_2 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    c = static_cast<double>(LS_x - x) / (LS_z - z);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x == LS_z - z && LS_y - y < LS_z - z)
                {
                    const auto key_1 = hashFunction(x + 1, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z + 1);
                    }
                    c = static_cast<double>(LS_y - y) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x == LS_z - z && LS_y - y > LS_z - z)
                {
                    const auto key_2 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    c = static_cast<double>(LS_z - z) / (LS_y - y);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y == LS_z - z && LS_x - x < LS_y - y)
                {
                    const auto key_1 = hashFunction(x, y + 1, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z + 1);
                    }
                    c = static_cast<double>(LS_x - x) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y == LS_z - z && LS_x - x > LS_y - y)
                {
                    const auto key_2 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    c = static_cast<double>(LS_y - y) / (LS_x - x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x < LS_y - y && y - LS_y < LS_z - z)
                {
                    const auto key_1 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    const auto key_2 = hashFunction(x, y + 1, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z + 1);
                    }
                    c1 = static_cast<double>(LS_x - x + LS_y - y) / 2 * (LS_z - z);
                    c2 = static_cast<double>(LS_x - x) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_x - x < LS_z - z && LS_z - z < LS_y - y)
                {
                    const auto key_1 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    const auto key_2 = hashFunction(x, y + 1, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z + 1);
                    }
                    c1 = static_cast<double>(LS_z - z + LS_x - x) / 2 * (LS_y - y);
                    c2 = static_cast<double>(LS_z - z) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y < LS_z - z && LS_z - z < LS_x - x)
                {
                    const auto key_1 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    const auto key_2 = hashFunction(x + 1, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z + 1);
                    }
                    c1 = static_cast<double>(LS_y - y + LS_z - z) / 2 * (LS_x - x);
                    c2 = static_cast<double>(LS_y - y) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y < LS_x - x && LS_x - x < LS_z - z)
                {
                    const auto key_1 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    const auto key_2 = hashFunction(x + 1, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z + 1);
                    }
                    c1 = static_cast<double>(LS_y - y + LS_x - x) / 2 * (LS_z - z);
                    c2 = static_cast<double>(LS_y - y) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_z - z < LS_y - y && LS_y - y < LS_x - x)
                {
                    const auto key_1 = hashFunction(x + 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y, z);
                    }
                    const auto key_2 = hashFunction(x + 1, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y + 1, z);
                    }
                    c1 = static_cast<double>(LS_z - z + LS_y - y) / 2 * (LS_x - x);
                    c2 = static_cast<double>(LS_z - z) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_z - z < LS_x - x && LS_x - x < LS_y - y)
                {
                    const auto key_1 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    const auto key_2 = hashFunction(x + 1, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x + 1, y + 1, z);
                    }
                    c1 = static_cast<double>(LS_z - z + LS_x - x) / 2 * (LS_y - y);
                    c2 = static_cast<double>(LS_z - z) / (LS_x - x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
            }
            // Q8
            else if ((x - LS_x > 0) && (y - LS_y < 0) && (z - LS_z < 0)) {
                const auto key = hashFunction(x - 1, y + 1, z + 1, lightSourceNumber);
                if (!visibilityHashMap_.count(key)) {
                    updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y + 1, z + 1);
                }
                if (x - LS_x == LS_y - y && LS_y - y == LS_z - z) {
                    v = visibilityHashMap_.at(key);
                }
                else if (x - LS_x == LS_y - y && LS_z - z < LS_y - y)
                {
                    const auto key_1 = hashFunction(x - 1, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y + 1, z);
                    }
                    c = static_cast<double>(LS_z - z) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x == LS_y - y && LS_z - z > LS_y - y)
                {
                    const auto key_2 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    c = static_cast<double>(LS_y - y) / (LS_z - z);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x == LS_z - z && LS_y - y < LS_z - z)
                {
                    const auto key_1 = hashFunction(x - 1, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z + 1);
                    }
                    c = static_cast<double>(LS_y - y) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x == LS_z - z && LS_y - y > LS_z - z)
                {
                    const auto key_2 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    c = static_cast<double>(LS_z - z) / (LS_y - y);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y == LS_z - z && x - LS_x < LS_y - y)
                {
                    const auto key_1 = hashFunction(x, y + 1, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z + 1);
                    }
                    c = static_cast<double>(x - LS_x) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    v = v1 - c * (v1 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y == LS_z - z && x - LS_x > LS_y - y)
                {
                    const auto key_2 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    c = static_cast<double>(LS_y - y) / (x - LS_x);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v2 - c * (v2 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x < LS_y - y && LS_y - y < LS_z - z)
                {
                    const auto key_1 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    const auto key_2 = hashFunction(x, y + 1, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z + 1);
                    }
                    c1 = static_cast<double>(x - LS_x + LS_y - y) / 2 * (LS_z - z);
                    c2 = static_cast<double>(x - LS_x) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (x - LS_x < LS_z - z && LS_z - z < LS_y - y)
                {
                    const auto key_1 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    const auto key_2 = hashFunction(x, y + 1, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z + 1);
                    }
                    c1 = static_cast<double>(x - LS_x + LS_z - z) / 2 * (LS_y - y);
                    c2 = static_cast<double>(x - LS_x) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y < LS_z - z && LS_z - z < x - LS_x)
                {
                    const auto key_1 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    const auto key_2 = hashFunction(x - 1, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z + 1);
                    }
                    c1 = static_cast<double>(LS_y - y + LS_z - z) / 2 * (x - LS_x);
                    c2 = static_cast<double>(LS_y - y) / (LS_z - z);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_y - y < x - LS_x && x - LS_x < LS_z - z)
                {
                    const auto key_1 = hashFunction(x, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y, z + 1);
                    }
                    const auto key_2 = hashFunction(x - 1, y, z + 1, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z + 1);
                    }
                    c1 = static_cast<double>(LS_y - y + x - LS_x) / 2 * (LS_z - z);
                    c2 = static_cast<double>(LS_y - y) / (x - LS_x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_z - z < x - LS_x && x - LS_x < LS_y - y)
                {
                    const auto key_1 = hashFunction(x, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x, y + 1, z);
                    }
                    const auto key_2 = hashFunction(x - 1, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y + 1, z);
                    }
                    c1 = static_cast<double>(LS_z - z + x - LS_x) / 2 * (LS_y - y);
                    c2 = static_cast<double>(LS_z - z) / (x - LS_x);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
                else if (LS_z - z < LS_y - y && LS_y - y < x - LS_x)
                {
                    const auto key_1 = hashFunction(x - 1, y, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_1)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y, z);
                    }
                    const auto key_2 = hashFunction(x - 1, y + 1, z, lightSourceNumber);
                    if (!visibilityHashMap_.count(key_2)) {
                        updatePointVisibility(lightSourceNumber, LS_x, LS_y, LS_z, x - 1, y + 1, z);
                    }
                    c1 = static_cast<double>(LS_z - z + LS_y - y) / 2 * (x - LS_x);
                    c2 = static_cast<double>(LS_z - z) / (LS_y - y);
                    double v1 = visibilityHashMap_.at(key_1);
                    double v2 = visibilityHashMap_.at(key_2);
                    v = v1 - c1 * (v1 - v2) - c1 * c2 * (v2 - visibilityHashMap_.at(key));
                }
            }
        }
        v = v * sharedVisibilityField_->get(x, y, z);
        key = hashFunction(x, y, z, lightSourceNumber);
        visibilityHashMap_[key] = v;

    }

    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    void Solver::reconstructPath(const Node& current,
        const std::string& methodName) {
        std::vector<point> resultingPath;
        int x = current.x, y = current.y, z = current.z;
        double t = cameFrom_(x, y, z);
        double t_old = std::numeric_limits<double>::infinity();
        while (t != t_old) {
            resultingPath.push_back({ x, y });
            t_old = t;
            if (methodName == "vstar") {
                x = lightSources_[t].x;
                y = lightSources_[t].y;
            }
            else if (methodName == "astar") {
                auto p = coordinatesAt(t);
                x = p.x;
                y = p.y;
            }
            t = cameFrom_(x, y, z);
        }
        resultingPath.push_back({ x, y,z });
        std::reverse(resultingPath.begin(), resultingPath.end());

        saveResults(resultingPath, methodName);

        // compute total distance
        double totalDistance = 0;
        for (size_t i = 0; i < resultingPath.size() - 1; ++i) {
            totalDistance +=
                evaluateDistance(resultingPath[i].x, resultingPath[i].y, resultingPath[i].z,
                    resultingPath[i + 1].x, resultingPath[i + 1].y, resultingPath[i + 1].z);
        }
        if (!Outsilent) {
            std::cout << methodName << " path length: " << totalDistance << std::endl;
        }

        return;
    }



    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    void Solver::saveVisibilityBasedSolverImage(const Field<double>& gScore) const {
        const int width = nx_;
        const int height = ny_;

        sf::Image image;
        image.create(width, height);

        double minVal = std::numeric_limits<double>::max();
        double maxVal = std::numeric_limits<double>::min();

        // Find min and max values in gScore for normalization
        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                double val = gScore(i, j, 0);
                if (val == std::numeric_limits<double>::infinity())
                    continue;
                if (val < minVal)
                    minVal = val;
                if (val > maxVal)
                    maxVal = val;
            }
        }

        // Generate the image
        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                if (sharedVisibilityField_->get(i, j, 0) < 1) {
                    image.setPixel(i, ny_ - 1 - j, sf::Color::Black);
                }
                else {
                    const double normalized_value = gScore(i, j, 0) / (maxVal - minVal);
                    sf::Color color = getColor(normalized_value);
                    image.setPixel(i, ny_ - 1 - j, color);
                }
            }
        }

        // color all initial frontline points as green circles with radius 10
        sf::Color color;
        color.a = 1;
        int x0, y0;
        int radius = 1;
        for (size_t k = 0; k < initial_frontline.size(); k += 3) {
            x0 = initial_frontline[k];
            y0 = initial_frontline[k + 1];

            for (int i = -radius; i <= radius; ++i) {
                for (int j = -radius; j <= radius; ++j) {
                    if (i * i + j * j <= radius * radius) {
                        if (x0 + i >= 0 && x0 + i < nx_ && y0 + j >= 0 && y0 + j < ny_) {
                            image.setPixel(x0 + i, ny_ - 1 - (y0 + j), color.Green);
                        }
                    }
                }
            }
        }

        if (1) {
            // compute the step size based on the max and min values
            int number_of_contour_lines = 10;
            double stepSize = (maxVal - minVal) / number_of_contour_lines;

            std::vector<double> contourLevels;
            for (double level = minVal; level <= maxVal; level += stepSize) {
                contourLevels.push_back(level);
            }
            // Draw contour lines on the image
            for (double level : contourLevels) {
                for (int i = 0; i < width; ++i) {
                    for (int j = 0; j < height; ++j) {
                        double value = gScore(i, j, 0);
                        if (std::abs(value - level) <= stepSize / 15) {
                            image.setPixel(i, ny_ - 1 - j, sf::Color::Black);
                        }
                    }
                }
            }
        }

        std::string outputPath = "./output/visibilityBasedSolver.png";
        image.saveToFile(outputPath);
    }

    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    void Solver::saveDistanceFunctionImage(const Field<double>& gScore) const {
        const int width = nx_;
        const int height = ny_;

        sf::Image image;
        image.create(width, height);

        double minVal = std::numeric_limits<double>::max();
        double maxVal = std::numeric_limits<double>::min();

        // Find min and max values in gScore for normalization
        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                double val = gScore(i, j, 0);
                if (val == std::numeric_limits<double>::infinity())
                    continue;
                if (val < minVal)
                    minVal = val;
                if (val > maxVal)
                    maxVal = val;
            }
        }

        // Generate the image
        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                if (sharedVisibilityField_->get(i, j, 0) < 1) {
                    image.setPixel(i, height - 1 - j, sf::Color::Black);
                }
                else {
                    const double normalized_value = gScore(i, j, 0) / (maxVal - minVal);
                    sf::Color color = getColor(normalized_value);
                    image.setPixel(i, height - 1 - j, color);
                }
            }
        }

        // compute the step size based on the max and min values
        if (0) {
            int number_of_contour_lines = 30 / 3;
            double stepSize = (maxVal - minVal) / number_of_contour_lines;

            std::vector<double> contourLevels;
            for (double level = minVal; level <= maxVal; level += stepSize) {
                contourLevels.push_back(level);
            }

            // Draw contour lines on the image
            for (double level : contourLevels) {
                for (int i = 0; i < width; ++i) {
                    for (int j = 0; j < height; ++j) {
                        double value = gScore(i, j, 0);
                        if (std::abs(value - level) <= stepSize / 15) {
                            image.setPixel(i, height - 1 - j, sf::Color::Black);
                        }
                    }
                }
            }
        }

        std::string outputPath = "output/distanceFunction.png";
        image.saveToFile(outputPath);
    }

    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    void Solver::saveCameFromImage(const Field<size_t>& cameFrom) const
    {
        const int width = nx_;
        const int height = ny_;

        sf::Image image;
        image.create(width, height);

        double minVal = std::numeric_limits<double>::max();
        double maxVal = std::numeric_limits<double>::min();

        // Find min and max values in gScore for normalization
        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                double val = cameFrom(i, j, 0);
                if (val == std::numeric_limits<double>::infinity())
                    continue;
                if (val < minVal)
                    minVal = val;
                if (val > maxVal)
                    maxVal = val;
            }
        }

        // Generate the image
        for (int i = 0; i < width; ++i) {
            for (int j = 0; j < height; ++j) {
                if (sharedVisibilityField_->get(i, j, 0) < 1) {
                    image.setPixel(i, height - 1 - j, sf::Color::Black);
                }
                else {
                    const double normalized_value = cameFrom(i, j, 0) / (maxVal - minVal);
                    sf::Color color = getColor(normalized_value);
                    image.setPixel(i, height - 1 - j, color);
                }
            }
        }

        // compute the step size based on the max and min values
        if (0) {
            int number_of_contour_lines = 30 / 3;
            double stepSize = (maxVal - minVal) / number_of_contour_lines;

            std::vector<double> contourLevels;
            for (double level = minVal; level <= maxVal; level += stepSize) {
                contourLevels.push_back(level);
            }

            // Draw contour lines on the image
            for (double level : contourLevels) {
                for (int i = 0; i < width; ++i) {
                    for (int j = 0; j < height; ++j) {
                        double value = cameFrom(i, j, 0);
                        if (std::abs(value - level) <= stepSize / 15) {
                            image.setPixel(i, height - 1 - j, sf::Color::Black);
                        }
                    }
                }
            }
        }

        std::string outputPath = "output/CamFrom.png";
        image.saveToFile(outputPath);
    }

    /*****************************************************************************/
    /*****************************************************************************/
    /*****************************************************************************/
    void Solver::saveResults(const std::vector<point>& resultingPath,
        const std::string& methodName) const {
        namespace fs = std::filesystem;

        // Define the path to the output file
        std::string outputFilePath = "./output/" + methodName + ".txt";

        // Check if the directory exists, and create it if it doesn't
        fs::path directory = fs::path(outputFilePath).parent_path();
        if (!fs::exists(directory)) {
            if (!fs::create_directories(directory)) {
                std::cerr << "Failed to create directory " << directory.string()
                    << std::endl;
                return;
            }
        }

        if (methodName == "distanceFunction") {
            // save distance function
            if (1) {
                std::fstream of(outputFilePath, std::ios::out | std::ios::trunc);
                if (!of.is_open()) {
                    std::cerr << "Failed to open output file " << outputFilePath
                        << std::endl;
                    return;
                }
                std::ostream& os = of;
                for (int j = ny_ - 1; j >= 0; --j) {
                    for (size_t i = 0; i < nx_; ++i) {
                        os << gScore_(i, j, 0) << " ";
                    }
                    os << "\n";
                }
                of.close();
                if (!Outsilent) {
                    std::cout << "Saved " + methodName << std::endl;
                }
            }
            return;
        }

        if (methodName == "visibilityBased") {
            std::fstream of(outputFilePath, std::ios::out | std::ios::trunc);
            if (!of.is_open()) {
                std::cerr << "Failed to open output file " << outputFilePath << std::endl;
                return;
            }
            std::ostream& os = of;
            for (int j = ny_ - 1; j >= 0; --j) {
                for (size_t i = 0; i < nx_; ++i) {
                    os << gScore_(i, j, 0) << " ";
                }
                os << "\n";
            }
            of.close();
            if (!Outsilent) {
                std::cout << "Saved " + methodName + " solution" << std::endl;
            }

            if (1) {
                outputFilePath = "./output/" + methodName + "_cameFrom.txt";
                std::fstream of1(outputFilePath, std::ios::out | std::ios::trunc);
                if (!of1.is_open()) {
                    std::cerr << "Failed to open output file " << outputFilePath
                        << std::endl;
                    return;
                }
                std::ostream& os1 = of1;
                for (int j = ny_ - 1; j >= 0; --j) {
                    for (size_t i = 0; i < nx_; ++i) {
                        os1 << cameFrom_(i, j, 0) << " ";
                    }
                    os1 << "\n";
                }
                of1.close();
                if (!Outsilent) {
                    std::cout << "Saved " + methodName + " cameFrom_" << std::endl;
                }
            }

            if (1) {
                outputFilePath = "./output/" + methodName + "_lightSources.txt";
                std::fstream of3(outputFilePath, std::ios::out | std::ios::trunc);
                if (!of3.is_open()) {
                    std::cerr << "Failed to open output file " << outputFilePath
                        << std::endl;
                    return;
                }
                std::ostream& os = of3;
                for (size_t i = 0; i < nb_of_sources_; ++i) {
                    os << lightSources_[i].x << " " << ny_ - 1 - lightSources_[i].y;
                    os << "\n";
                }
                of3.close();
                if (!Outsilent) {
                    std::cout << "Saved " + methodName + " lightSources" << std::endl;
                }
            }
            return;
        }

        // Save resulting path
        if (1) {
            outputFilePath = "./output/" + methodName + "_path.txt";
            std::fstream of(outputFilePath, std::ios::out | std::ios::trunc);
            if (!of.is_open()) {
                std::cerr << "Failed to open output file " << outputFilePath << std::endl;
                return;
            }
            std::ostream& os = of;
            for (size_t i = 0; i < resultingPath.size(); ++i) {
                os << resultingPath[i].x << " " << ny_ - 1 - resultingPath[i].y;
                os << "\n";
            }
            of.close();
            if (!Outsilent) {
                std::cout << "Saved " + methodName + " path" << std::endl;
            }
        }

        // Save gScore_
        if (1) {
            outputFilePath = "./output/" + methodName + "_gScore.txt";
            std::fstream of1(outputFilePath, std::ios::out | std::ios::trunc);
            if (!of1.is_open()) {
                std::cerr << "Failed to open output file " << outputFilePath << std::endl;
                return;
            }
            std::ostream& os = of1;
            for (int j = ny_ - 1; j >= 0; --j) {
                for (size_t i = 0; i < nx_; ++i) {
                    os << gScore_(i, j, 0) << " ";
                }
                os << "\n";
            }
            of1.close();
            if (!Outsilent) {
                std::cout << "Saved " + methodName + " gScore" << std::endl;
            }
        }

        // Save fScore_
        if (1) {
            outputFilePath = "./output/" + methodName + "_fScore.txt";
            std::fstream of2(outputFilePath, std::ios::out | std::ios::trunc);
            if (!of2.is_open()) {
                std::cerr << "Failed to open output file " << outputFilePath << std::endl;
                return;
            }
            std::ostream& os = of2;
            for (int j = ny_ - 1; j >= 0; --j) {
                for (size_t i = 0; i < nx_; ++i) {
                    os << fScore_(i, j, 0) << " ";
                }
                os << "\n";
            }
            of2.close();
            if (!Outsilent) {
                std::cout << "Saved " + methodName + " fScore" << std::endl;
            }
        }

        // Save lightsources
        if (1) {
            outputFilePath = "./output/" + methodName + "_lightSources.txt";
            std::fstream of3(outputFilePath, std::ios::out | std::ios::trunc);
            if (!of3.is_open()) {
                std::cerr << "Failed to open output file " << outputFilePath << std::endl;
                return;
            }
            std::ostream& os = of3;
            for (size_t i = 0; i < nb_of_sources_; ++i) {
                os << lightSources_[i].x << " " << ny_ - 1 - lightSources_[i].y;
                os << "\n";
            }
            of3.close();
            if (!Outsilent) {
                std::cout << "Saved " + methodName + " lightSources" << std::endl;
            }
        }
    }

} // namespace vbm
