#include <iostream>
#include <string>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <fstream>
#include <omp.h>

using namespace std;

namespace std
{
    template <>
    struct hash<pair<string, string>>
    {
        size_t operator()(const pair<string, string> &p) const noexcept
        {
            size_t h1 = hash<string>{}(p.first);
            size_t h2 = hash<string>{}(p.second);
            return h1 ^ (h2 << 1);
        }
    };
}

int match(string seq1, string seq2, int index1, int index2, string direction)
{
    int matchCount = 0;
    if (direction == "lijevo")
    {
        for (int i = index1; i <= index2; i++)
        {
            if (seq1[i] == seq2[i - 1])
                matchCount++;
            if (seq1[i] == seq2[i])
                matchCount--;
        }
    }
    else
    {
        for (int i = index1; i <= index2; i++)
        {
            if (seq1[i] == seq2[i + 1])
                matchCount++;
            if (seq1[i] == seq2[i])
                matchCount--;
        }
    }
    return matchCount;
}

tuple<string, string, int> provjeri_susjede(string seq1, string seq2, int index1, int index2)
{
    if (index1 > 0 && index2 < min(seq1.length(), seq2.length()) - 1)
    {
        if (seq1[index1 - 1] == '-' && seq1[index2 + 1] != '-' && seq2[index1 - 1] != '-' && seq2[index2 + 1] == '-')
        {
            if (match(seq1, seq2, index1, index2, "lijevo") >= 0)
            {
                seq1.erase(seq1.begin() + index1 - 1);
                seq2.erase(seq2.begin() + index2 + 1);
                provjeri_susjede(seq1, seq2, index1 - 1, index2);
            }
        }
        else if (seq1[index1 - 1] != '-' && seq1[index2 + 1] == '-' && seq2[index1 - 1] == '-' && seq2[index2 + 1] != '-')
        {
            if (match(seq1, seq2, index1, index2, "desno") >= 0)
            {
                seq1.erase(seq1.begin() + index2 + 1);
                seq2.erase(seq2.begin() + index1 - 1);
                provjeri_susjede(seq1, seq2, index1 - 1, index2);
            }
        }
    }
    return {seq1, seq2, index2};
}

pair<string, string> poravnaj(string seq1, string seq2)
{
    for (int i = 0; i < min(seq1.length(), seq2.length()); i++)
    {
        if (seq1[i] == '-' && seq2[i] == '-')
        {
            seq1.erase(seq1.begin() + i);
            seq2.erase(seq2.begin() + i);
            i--;
        }
    }

    for (int i = 1; i < min(seq1.length(), seq2.length()); i++)
    {
        if (seq1[i] == '-' && seq2[i] != '-' && seq1[i - 1] != '-' && seq2[i - 1] == '-')
        {
            seq1.erase(seq1.begin() + i);
            seq2.erase(seq2.begin() + i - 1);
            i--;
            auto [new_seq1, new_seq2, new_i] = provjeri_susjede(seq1, seq2, i, i);
            seq1 = new_seq1;
            seq2 = new_seq2;
            i = new_i;
        }
        else if (seq2[i] == '-' && seq1[i] != '-' && seq2[i - 1] != '-' && seq1[i - 1] == '-')
        {
            seq2.erase(seq2.begin() + i);
            seq1.erase(seq1.begin() + i - 1);
            i--;
            auto [new_seq1, new_seq2, new_i] = provjeri_susjede(seq1, seq2, i, i);
            seq1 = new_seq1;
            seq2 = new_seq2;
            i = new_i;
        }
    }

    return {seq1, seq2};
}

int main(int argc, char *argv[])
{
    string inputFileName = argv[1];
    string outputFileName = argv[2];
    ifstream file(inputFileName);
    if (!file.is_open())
    {
        cerr << "Error opening file: " << inputFileName << endl;
        return 1;
    }
    unordered_map<string, string> sequences;
    vector<string> names;
    string line;
    string id;
    while (getline(file, line))
    {
        if (line[0] == '>')
        {
            id = line.substr(1);
            names.push_back(id);
            sequences[id] = "";
            continue;
        }
        sequences[id] += line;
    }
    ofstream outputFile(outputFileName, ios::binary);
    if (!outputFile.is_open())
    {
        cerr << "Error opening file: " << outputFileName << endl;
        return 1;
    }
    // printf("names size: %d\n", names.size());
    int totalWrites = 0;
#pragma omp parallel
    {
        unordered_map<pair<string, string>, pair<string, string>> localMap;
        int numPairs = 0;

        for (int i = 0; i < names.size(); i++)
        {
#pragma omp for nowait
            for (int j = i + 1; j < names.size(); j++)
            {
                localMap[{names[i], names[j]}] = poravnaj(sequences[names[i]], sequences[names[j]]);
                numPairs++;
                if (numPairs % 100000 == 0)
                {
#pragma omp critical
                    {
                        printf("Merging results from thread %d... %d\n", omp_get_thread_num(), numPairs);
                        int cleanedCount = localMap.size();
                        outputFile.write(reinterpret_cast<char *>(&cleanedCount), sizeof(int));
                        for (auto &seq : localMap)
                        {
                            outputFile.write(reinterpret_cast<const char *>(seq.first.first.c_str()), seq.first.first.size() + 1);
                            outputFile.write(reinterpret_cast<const char *>(seq.first.second.c_str()), seq.first.second.size() + 1);
                            outputFile.write(reinterpret_cast<const char *>(seq.second.first.c_str()), seq.second.first.size() + 1);
                            outputFile.write(reinterpret_cast<const char *>(seq.second.second.c_str()), seq.second.second.size() + 1);
                        }
                        localMap.clear();
                        numPairs = 0;
                        totalWrites++;
                    }
                }
            }
        }
#pragma omp critical
        {
            printf("Merging results from thread %d...\n", omp_get_thread_num());
            int cleanedCount = localMap.size();
            outputFile.write(reinterpret_cast<char *>(&cleanedCount), sizeof(int));
            for (auto &seq : localMap)
            {
                outputFile.write(reinterpret_cast<const char *>(seq.first.first.c_str()), seq.first.first.size() + 1);
                outputFile.write(reinterpret_cast<const char *>(seq.first.second.c_str()), seq.first.second.size() + 1);
                outputFile.write(reinterpret_cast<const char *>(seq.second.first.c_str()), seq.second.first.size() + 1);
                outputFile.write(reinterpret_cast<const char *>(seq.second.second.c_str()), seq.second.second.size() + 1);
            }
            localMap.clear();
            totalWrites++;
        }
    }
    printf("Total writes: %d\n", totalWrites);
    file.close();
    outputFile.close();

    return 0;
}
