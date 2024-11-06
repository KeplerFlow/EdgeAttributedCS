#ifndef MAPPING_H
#define MAPPING_H

#include <string>
#include <unordered_map>
#include <fstream>
#include <iostream>

// 将字符串名称映射到整数索引和将整数索引映射回字符串名称
class Mapping {
private:
    std::unordered_map<std::string, int> nameIndex;
    std::unordered_map<int, std::string> indexName;

public:
    Mapping() {
        //从文件中读取每一行，每一行包含一个整数索引和对应的字符串名称
        std::cout << "Reading mapping..." << std::endl;
        std::ifstream fileReader("map.txt");
        if (!fileReader.is_open()) {
            std::cout << "File not found." << std::endl;
            exit(0);
        }

        std::string line;
        while (std::getline(fileReader, line)) {
            size_t pos = line.find('\t');
            std::string indexStr = line.substr(0, pos);
            std::string name = line.substr(pos + 1);
            int index = std::stoi(indexStr);

            nameIndex[name] = index;
            indexName[index] = name;
        }

        fileReader.close();
        std::cout << "Mapping read." << std::endl;
    }

    int getIndex(const std::string& name) {
        return nameIndex[name];
    }

    std::string getName(int index) {
        return indexName[index];
    }

    int getDimension() {
        return nameIndex.size();
    }
};

#endif /* MAPPING_H */