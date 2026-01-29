#include <iostream>
#include <string>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <fstream>
#include <omp.h>
#include <cstdio>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

// hash for pair<string,string>
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

bool allowed(char c)
{
    return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == '-';
}

int match(const string &seq1, const string &seq2, int index1, int index2, const string &direction)
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
    if (index1 > 0 && index2 < (int)min(seq1.length(), seq2.length()) - 1)
    {
        if (seq1[index1 - 1] == '-' && seq1[index2 + 1] != '-' && seq2[index1 - 1] != '-' && seq2[index2 + 1] == '-')
        {
            if (match(seq1, seq2, index1, index2, "lijevo") >= 0)
            {
                seq1.erase(seq1.begin() + index1 - 1);
                seq2.erase(seq2.begin() + index2 + 1);
                return provjeri_susjede(seq1, seq2, index1 - 1, index2);
            }
        }
        else if (seq1[index1 - 1] != '-' && seq1[index2 + 1] == '-' && seq2[index1 - 1] == '-' && seq2[index2 + 1] != '-')
        {
            if (match(seq1, seq2, index1, index2, "desno") >= 0)
            {
                seq1.erase(seq1.begin() + index2 + 1);
                seq2.erase(seq2.begin() + index1 - 1);
                return provjeri_susjede(seq1, seq2, index1 - 1, index2);
            }
        }
    }
    return {seq1, seq2, index2};
}

pair<string, string> poravnaj(string seq1, string seq2)
{
    for (size_t i = 0; i < min(seq1.length(), seq2.length()); ++i)
    {
        if (seq1[i] == '-' && seq2[i] == '-')
        {
            seq1.erase(seq1.begin() + i);
            seq2.erase(seq2.begin() + i);
            --i;
        }
    }

    for (int i = 1; i < (int)min(seq1.length(), seq2.length()); ++i)
    {
        if (seq1[i] == '-' && seq2[i] != '-' && seq1[i - 1] != '-' && seq2[i - 1] == '-')
        {
            seq1.erase(seq1.begin() + i);
            seq2.erase(seq2.begin() + i - 1);
            --i;
            auto [new_seq1, new_seq2, new_i] = provjeri_susjede(seq1, seq2, i, i);
            seq1 = move(new_seq1);
            seq2 = move(new_seq2);
            i = new_i;
        }
        else if (seq2[i] == '-' && seq1[i] != '-' && seq2[i - 1] != '-' && seq1[i - 1] == '-')
        {
            seq2.erase(seq2.begin() + i);
            seq1.erase(seq1.begin() + i - 1);
            --i;
            auto [new_seq1, new_seq2, new_i] = provjeri_susjede(seq1, seq2, i, i);
            seq1 = move(new_seq1);
            seq2 = move(new_seq2);
            i = new_i;
        }
    }
    return {seq1, seq2};
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        cerr << "Usage: writer <input.fasta> <output.bin>\n";
        return 1;
    }
    string inputFileName = argv[1];
    string outputFileName = argv[2];

    ifstream file(inputFileName);
    if (!file.is_open())
    {
        cerr << "Error opening file: " << inputFileName << "\n";
        return 1;
    }

    unordered_map<string, string> sequences;
    vector<string> names;
    string line;
    string id;
    while (getline(file, line))
    {
        if (line.empty())
            continue;
        if (line[0] == '>')
        {
            id = line.substr(1);
            names.push_back(id);
            sequences[id] = "";
        }
        else
        {
            string clean;
            for (char c : line)
            {
                if (allowed(c))
                    clean.push_back(c);
            }
            sequences[id] += clean;
        }
    }
    file.close();
    size_t N = 1000; // names.size();
    if (N < 2)
    {
        cerr << "Need at least 2 sequences\n";
        return 1;
    }

    // privremeni direktorij
    string tmpdir = outputFileName + ".parts";
    fs::create_directory(tmpdir);

    // Paralelno poravnanje
    int max_threads = omp_get_max_threads();
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        vector<pair<pair<string, string>, pair<string, string>>> local_entries;
#pragma omp for schedule(dynamic)
        for (int i = 0; i < (int)N; ++i)
        {
            for (int j = i + 1; j < (int)N; ++j)
            {
                auto aligned = poravnaj(sequences[names[i]], sequences[names[j]]);
                local_entries.emplace_back(make_pair(names[i], names[j]), aligned);
            }
        }

        // zapiši svoj dio u privremenu datoteku
        string partName = tmpdir + "/part_" + to_string(tid) + ".bin";
        ofstream partOut(partName, ios::binary);
        if (!partOut.is_open())
        {
#pragma omp critical
            cerr << "Thread " << tid << " could not open part file: " << partName << "\n";
        }
        else
        {
            int cleanedCount = (int)local_entries.size();
            partOut.write(reinterpret_cast<const char *>(&cleanedCount), sizeof(int));
            for (auto &entry : local_entries)
            {
                const string &n1 = entry.first.first;
                const string &n2 = entry.first.second;
                const string &s1 = entry.second.first;
                const string &s2 = entry.second.second;
                partOut.write(n1.c_str(), n1.size() + 1);
                partOut.write(n2.c_str(), n2.size() + 1);
                partOut.write(s1.c_str(), s1.size() + 1);
                partOut.write(s2.c_str(), s2.size() + 1);
            }
            partOut.close();
        }
    }

    // Spoji dijelove u konačni izlaz
    ofstream finalOut(outputFileName, ios::binary);
    if (!finalOut.is_open())
    {
        cerr << "Error opening final output: " << outputFileName << "\n";
        return 1;
    }

    for (int tid = 0; tid < max_threads; ++tid)
    {
        string partName = tmpdir + "/part_" + to_string(tid) + ".bin";
        if (!fs::exists(partName))
            continue;
        ifstream partIn(partName, ios::binary);
        finalOut << partIn.rdbuf();
        partIn.close();
        fs::remove(partName);
    }
    finalOut.close();
    fs::remove(tmpdir);

    cout << "Writer finished. Output: " << outputFileName << "\n";
    return 0;
}
