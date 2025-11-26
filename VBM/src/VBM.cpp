#include "VBM/VBM.hpp"

namespace vbm {

    void DistanceFunction(const size_t nx, const size_t ny, const size_t nz, const std::vector<bool>& visibilityField,
        const std::vector<float>& gradientX, const std::vector<float>& gradientY, const std::vector<float>& gradientZ, std::vector<float>& unsignedDis) {

        vbm::Field<double> field(nx, ny, nz, 1.0);
        for (size_t i = 0; i < nx * ny * nz; i++) {
            double value;
            if (visibilityField[i])
                value = 0.0f;
            else
                value = 1.0f;

            field.set(i, value);
        }

        std::shared_ptr<vbm::Field<double>> sharedVisbilityField = std::make_shared<vbm::Field<double>>(field);
        std::shared_ptr<vbm::Field<double>> sharedSpeedField = std::make_shared<vbm::Field<double>>(field);


        vbm::Solver sol = vbm::Solver(sharedVisbilityField, sharedSpeedField);

        // Set the gradient fields
        vbm::Field<vector3f> gradientField(nx, ny, nz, vector3f(0, 0, 0));
        for (size_t i = 0; i < nx * ny * nz; i++) {
            gradientField.set(i, vector3f(gradientX[i], gradientY[i], gradientZ[i]));
        }

        field = sol.computeDistanceFunction(gradientField);

        for (size_t i = 0; i < nx * ny * nz; i++) {
            unsignedDis[i] = field.get(i);
        }

        return;
    }


    void DistanceFunction(const size_t nx, const size_t ny, const size_t nz, const std::vector<bool>& visibilityField, std::vector<float>& unsignedDis) {

        vbm::Field<double> field(nx, ny, nz, 1.0);
        for (size_t i = 0; i < nx * ny * nz; i++) {
            double value;
            if (visibilityField[i])
                value = 0.0f;
            else
                value = 1.0f;

            field.set(i, value);
        }

        std::shared_ptr<vbm::Field<double>> sharedVisbilityField = std::make_shared<vbm::Field<double>>(field);
        std::shared_ptr<vbm::Field<double>> sharedSpeedField = std::make_shared<vbm::Field<double>>(field);


        vbm::Solver sol = vbm::Solver(sharedVisbilityField, sharedSpeedField);

        field = sol.computeDistanceFunction();

        for (size_t i = 0; i < nx * ny * nz; i++) {
            unsignedDis[i] = field.get(i);
        }

        return;
    }

    void DistanceFunction(const size_t nx, const size_t ny, const size_t nz, const std::vector<bool>& visibilityField, std::vector<int>& unsignedDis) {

        vbm::Field<double> field(nx, ny, nz, 1.0);
        for (size_t i = 0; i < nx * ny * nz; i++) {
            double value;
            if (visibilityField[i])
                value = 0.0f;
            else
                value = 1.0f;

            field.set(i, value);
        }

        std::shared_ptr<vbm::Field<double>> sharedVisbilityField = std::make_shared<vbm::Field<double>>(field);
        std::shared_ptr<vbm::Field<double>> sharedSpeedField = std::make_shared<vbm::Field<double>>(field);


        vbm::Solver sol = vbm::Solver(sharedVisbilityField, sharedSpeedField);

        field = sol.computeDistanceFunction();

        for (size_t i = 0; i < nx * ny * nz; i++) {
            unsignedDis[i] = field.get(i);
        }

        return;
    }

    void DistanceFunctionWithSources(const size_t nx, const size_t ny, const size_t nz, 
        const std::vector<bool>& visibilityField, std::vector<int>& unsignedDis,std::vector<int>& sources, std::vector<std::vector<int>>& sourcePoints) {

        vbm::Field<double> field(nx, ny, nz, 1.0);
        for (size_t i = 0; i < nx * ny * nz; i++) {
            double value;
            if (visibilityField[i])
                value = 0.0f;
            else
                value = 1.0f;

            field.set(i, value);
        }

        std::shared_ptr<vbm::Field<double>> sharedVisbilityField = std::make_shared<vbm::Field<double>>(field);
        std::shared_ptr<vbm::Field<double>> sharedSpeedField = std::make_shared<vbm::Field<double>>(field);


        vbm::Solver sol = vbm::Solver(sharedVisbilityField, sharedSpeedField);

        field = sol.computeDistanceFunction();

        // 计算source点的距离场
        point* LightSources = sol.getLightSources();

        int numofSources = sol.getNbOfSources();

        for (size_t i = 0; i < numofSources; i++)
        {
            int x = LightSources[i].x;
            int y = LightSources[i].y;
            int z = LightSources[i].z;
            sourcePoints.push_back({ x,y,z });
        }

        vbm::Field<size_t> sourceField = sol.getCameFromField();
        for (size_t i = 0; i < nx * ny * nz; i++) {
            unsignedDis[i] = field.get(i);
            sources[i] = sourceField.get(i);
        }

        return;
    }


    // 计算internal visaul hull
    void computeIVH(const size_t nx, const size_t ny, const size_t nz, const std::vector<bool>& visibilityField, std::vector<int>& internalCH)
    {
        // 初始化一下internalCH 以及可见性场
        vbm::Field<double> field(nx, ny, nz, 1.0);
        for (size_t i = 0; i < nx * ny * nz; i++) {
            double value;
            if (visibilityField[i])
                value = 0.0f;
            else
                value = 1.0f;

            field.set(i, value);
            internalCH[i] = 1;
        }

        

        std::shared_ptr<vbm::Field<double>> sharedVisbilityField = std::make_shared<vbm::Field<double>>(field);
        std::shared_ptr<vbm::Field<double>> sharedSpeedField = std::make_shared<vbm::Field<double>>(field);


        vbm::Solver sol = vbm::Solver(sharedVisbilityField, sharedSpeedField);

        vbm::Field<size_t> cameFrom = sol.visibilityBasedSolver();

        const int numofSources = sol.getNbOfSources();


        for (size_t i = 0; i < nx * ny * nz; i++) {
            size_t label = cameFrom.get(i);
            if (label < numofSources && visibilityField[i] != 1) //说明该点可见
                internalCH[i] = 0;
        }

        return;
    }


    void Vector2Image(const size_t nx, const size_t ny, const size_t nz, const size_t ind, const std::vector<int>& gScore, const std::string& outputPath) {

        sf::Image image;
        image.create(nx, ny);

        double minVal = std::numeric_limits<double>::max();
        double maxVal = std::numeric_limits<double>::min();

        // Find min and max values in gScore for normalization
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                int index = i * ny * nz + j * nz + ind;
                double val = gScore[index];
                if (val == std::numeric_limits<double>::infinity())
                    continue;
                if (val < minVal)
                    minVal = val;
                if (val > maxVal)
                    maxVal = val;
            }
        }

        // Generate the image
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                int index = i * ny * nz + j * nz + ind;
                if (gScore[index] == 0) {
                    image.setPixel(i, ny - 1 - j, sf::Color::Black);
                }
                else {
                    const double normalized_value = gScore[index] / (maxVal - minVal);
                    sf::Color color = getColor(normalized_value);
                    image.setPixel(i, ny - 1 - j, color);
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
                for (int i = 0; i < nx; ++i) {
                    for (int j = 0; j < ny; ++j) {
                        int index = i * ny * nz + j * nz + ind;
                        double value = gScore[index];
                        if (std::abs(value - level) <= stepSize / 15) {
                            image.setPixel(i, ny - 1 - j, sf::Color::Black);
                        }
                    }
                }
            }
        }

        image.saveToFile(outputPath);
    }
}

