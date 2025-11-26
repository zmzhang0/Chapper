#pragma once

#include <map>
#include <list>
#include <mutex>
#include <utility>

template <typename Key, typename Value>
class LimitedSizeMap {
private:
    size_t maxSize; // 容量上限
    std::map<Key, Value> dataMap; // 存储数据的 map
    std::list<Key> accessOrder;  // 记录访问顺序
    std::mutex mapMutex;         // 保护数据结构的互斥锁

    // 删除最老的元素
    void evictIfNecessary() {
        if (dataMap.size() > maxSize) {
            Key oldestKey = accessOrder.front();
            accessOrder.pop_front();
            dataMap.erase(oldestKey);
        }
    }

public:
    // 构造函数，初始化容量上限
    explicit LimitedSizeMap(size_t capacity) : maxSize(capacity) {}

    // 设置maxSize
    void setCapacity(size_t capacity) {
        maxSize = capacity;
    }

    // 插入一个键值对
    void insert(const Key& key, const Value& value) {
        std::lock_guard<std::mutex> lock(mapMutex);

        // 如果键已存在，先删除旧值
        if (dataMap.find(key) != dataMap.end()) {
            accessOrder.remove(key);
        }

        // 插入新值并更新访问顺序
        dataMap[key] = value;
        accessOrder.push_back(key);

        // 检查是否需要驱逐旧元素
        evictIfNecessary();
    }

    // 查找一个键
    bool find(const Key& key, Value& result) {
        std::lock_guard<std::mutex> lock(mapMutex);

        auto it = dataMap.find(key);
        if (it != dataMap.end()) {
            // 更新访问顺序
            accessOrder.remove(key);
            accessOrder.push_back(key);
            result = it->second;
            return true;
        }
        return false;
    }

    // 获取当前元素数量
    size_t size() const {
        std::lock_guard<std::mutex> lock(mapMutex);
        return dataMap.size();
    }

    // 清空所有数据
    void clear() {
        std::lock_guard<std::mutex> lock(mapMutex);
        dataMap.clear();
        accessOrder.clear();
    }
};
