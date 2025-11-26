#ifndef FIELD_H
#define FIELD_H

#include <cstddef>
#include <memory>

namespace vbm {

    template <typename T> class Field {
    public:
        using size_t = std::size_t;

        // 默认构造函数
        Field() : nx_(0), ny_(0), nz_(0), size_(0), data_(nullptr) {}

        // 构造函数，初始化三维字段
        explicit Field(const size_t nx, const size_t ny, const size_t nz, const T default_value)
            : nx_(nx), ny_(ny), nz_(nz), size_(nx* ny* nz) {
            data_ = std::make_unique<T[]>(size_);
            for (size_t i = 0; i < size_; ++i) {
                data_[i] = default_value;
            }
        }

        // 访问元素
        inline T& operator()(size_t x, size_t y, size_t z) {
            return data_[x * ny_ * nz_ + y * nz_ + z];
        }

        inline const T& operator()(size_t x, size_t y, size_t z) const {
            return data_[x * ny_ * nz_ + y * nz_ + z];
        }

        // 设置元素
        void set(const size_t x, const size_t y, const size_t z, const T value) {
            data_[x * ny_ * nz_ + y * nz_ + z] = value;
        }

        // 设置元素
        void set(const size_t index, const T value) {
            data_[index] = value;
        } 

        // 获取元素
        const T get(const size_t x, const size_t y, const size_t z) const {
            return data_[x * ny_ * nz_ + y * nz_ + z];
        }
        const T get(const size_t i) const {
            return data_[i];
        }

        // 获取字段的大小
        const size_t nx() const { return nx_; }
        const size_t ny() const { return ny_; }
        const size_t nz() const { return nz_; }

        // 复制构造函数
        Field(const Field& other)
            : nx_(other.nx_), ny_(other.ny_), nz_(other.nz_), size_(other.size_) {
            data_ = std::make_unique<T[]>(size_);
            for (size_t i = 0; i < size_; ++i) {
                data_[i] = other.data_[i];
            }
        }

        // 复制赋值操作符
        Field& operator=(const Field& other) {
            if (this != &other) {
                nx_ = other.nx_;
                ny_ = other.ny_;
                nz_ = other.nz_;
                size_ = other.size_;
                data_ = std::make_unique<T[]>(size_);
                for (size_t i = 0; i < size_; ++i) {
                    data_[i] = other.data_[i];
                }
            }
            return *this;
        }

        // 移动构造函数
        Field(Field&& other) noexcept
            : nx_(other.nx_), ny_(other.ny_), nz_(other.nz_), size_(other.size_) {
            data_ = std::move(other.data_);
        }

        // 移动赋值操作符
        Field& operator=(Field&& other) noexcept {
            if (this != &other) {
                nx_ = other.nx_;
                ny_ = other.ny_;
                nz_ = other.nz_;
                size_ = other.size_;
                data_ = std::move(other.data_);
            }
            return *this;
        }

        // 重置字段大小并初始化
        void reset(const size_t nx, const size_t ny, const size_t nz, const T default_value) {
            nx_ = nx;
            ny_ = ny;
            nz_ = nz;
            size_ = nx * ny * nz;
            data_ = std::make_unique<T[]>(size_);
            for (size_t i = 0; i < size_; ++i) {
                data_[i] = default_value;
            }
        }

        // 析构函数
        ~Field() = default;

    private:
        size_t nx_;  // x维大小
        size_t ny_;  // y维大小
        size_t nz_;  // z维大小
        size_t size_; // 总大小 nx * ny * nz
        std::unique_ptr<T[]> data_; // 数据指针
    };

} // namespace vbm

#endif // FIELD_H
