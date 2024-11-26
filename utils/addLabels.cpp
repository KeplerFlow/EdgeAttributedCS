#include <fstream>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <vector>
int main() {
    std::string inputFilename = "/home/asc23/lcy/graph/Final/Dataset/Email-Enron.txt";  // 输入文件名
    std::string outputFilename = "/home/asc23/lcy/graph/Final/Dataset/Email-Enron-labels.txt"; // 输出文件名

     std::ifstream inputFile(inputFilename);   // 打开输入文件
    std::ofstream outputFile(outputFilename); // 打开输出文件

    if (!inputFile.is_open() || !outputFile.is_open()) {
        std::cerr << "Error opening files!" << std::endl;
        return 1;
    }

    std::srand(std::time(nullptr)); // 设置随机数种子

    std::string line;
    while (std::getline(inputFile, line)) {
        if (!line.empty()) {
            int randomNum = std::rand() % 8; // 生成0到7之间的随机数

            // 使用stringstream来处理行
            std::stringstream ss(line);
            std::string word;
            std::vector<std::string> words;

            // 将每个单词（以空格或制表符分隔）存入vector
            while (ss >> word) {
                words.push_back(word);
            }

            // 输出处理过的行和随机数
            for (size_t i = 0; i < words.size(); ++i) {
                if (i > 0) {
                    outputFile << " ";
                }
                outputFile << words[i];
            }
            outputFile << " " << randomNum << std::endl; // 追加随机数
        }
    }

    inputFile.close();
    outputFile.close();

    std::cout << "Process completed. Output written to " << outputFilename << std::endl;

    return 0;
}
