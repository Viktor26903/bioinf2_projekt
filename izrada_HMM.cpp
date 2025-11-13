#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>

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

string readCString(ifstream &input)
{
    string result;
    char ch;
    while (input.get(ch))
    {
        if (ch == '\0')
            break;
        result += ch;
    }
    return result;
}

void read_and_predict_hmm(string &input, unordered_map<pair<string, string>, pair<string, string>> &alignments, vector<vector<double>> &a, vector<unordered_map<string, double>> &e)
{
    ifstream file(input, ios::binary);
    if (!file.is_open())
    {
        cerr << "Error opening file: " << input << endl;
        exit(1);
    }
    int beg_match = 0, beg_ins = 0, beg_del = 0, match_match = 0, match_ins = 0, match_del = 0, match_end = 0, ins_match = 0, ins_ins = 0, ins_del = 0, ins_end = 0, del_match = 0, del_del = 0, del_ins = 0, del_end = 0;
    int beg = 0, ins = 0, del = 0, match = 0;
    while (true)
    {
        int cleanedCount;
        if (!file.read(reinterpret_cast<char *>(&cleanedCount), sizeof(int)))
            break;

        cout << "Block with " << cleanedCount << " entries:\n";

        for (int i = 0; i < cleanedCount; i++)
        {
            string name1 = readCString(file);
            string name2 = readCString(file);
            string seq1 = readCString(file);
            string seq2 = readCString(file);
            // cout << name1 << "\t" << name2 << "\n";
            alignments[{name1, name2}] = {seq1, seq2};
            // procijeni_hmm(seq1, seq2, beg, ins, del, match, match_match, match_ins, match_del, match_end, ins_match, ins_ins, ins_del, ins_end, del_match, del_del, del_ins, del_end);
            for (int i = 0; i < min(seq1.length(), seq2.length()); i++)
            {
                if (i == 0)
                {
                    if (seq1[i] != '-' && seq2[i] != '-')
                    {
                        beg_match++;
                        match++;
                    }
                    else if (seq1[i] == '-' && seq2[i] != '-')
                    {
                        beg_ins++;
                    }
                    else if (seq1[i] != '-' && seq2[i] == '-')
                    {
                        beg_del++;
                    }
                }
                else
                {
                    if (seq1[i] != '-' && seq2[i] != '-')
                    {
                        match++;
                        if (seq1[i - 1] != '-' && seq2[i - 1] != '-')
                            match_match++;
                        else if (seq1[i - 1] == '-' && seq2[i - 1] != '-')
                            ins_match++;
                        else if (seq1[i - 1] != '-' && seq2[i - 1] == '-')
                            del_match++;
                        if (i == min(seq1.length(), seq2.length()) - 1)
                            match_end++;
                    }
                    else if (seq1[i] == '-' && seq2[i] != '-')
                    {
                        ins++;
                        if (seq1[i - 1] != '-' && seq2[i - 1] != '-')
                            match_ins++;
                        else if (seq1[i - 1] == '-' && seq2[i - 1] != '-')
                            ins_ins++;
                        else if (seq1[i - 1] != '-' && seq2[i - 1] == '-')
                            del_ins++;
                        if (i == min(seq1.length(), seq2.length()) - 1)
                            ins_end++;
                    }
                    else if (seq1[i] != '-' && seq2[i] == '-')
                    {
                        del++;
                        if (seq1[i - 1] != '-' && seq2[i - 1] != '-')
                            match_del++;
                        else if (seq1[i - 1] == '-' && seq2[i - 1] != '-')
                            ins_del++;
                        else if (seq1[i - 1] != '-' && seq2[i - 1] == '-')
                            del_del++;
                        if (i == min(seq1.length(), seq2.length()) - 1)
                            del_end++;
                    }
                }
            }

            if (!file)
            {
                cerr << "Unexpected EOF while reading entry " << i << endl;
                exit(1);
            }
        }
        beg += cleanedCount;
    }

    a[0][0] = 0.0;
    a[0][1] = (double)beg_match / beg;
    a[0][2] = (double)beg_ins / beg;
    a[0][3] = (double)beg_del / beg;
    a[0][4] = 0.0;

    a[1][0] = 0.0;
    a[1][1] = (double)match_match / (match - match_end);
    a[1][2] = (double)match_ins / (match - match_end);
    a[1][3] = (double)match_del / (match - match_end);
    a[1][4] = (double)match_end / (match - match_end);

    a[2][0] = 0.0;
    a[2][1] = (double)ins_match / (ins - ins_end);
    a[2][2] = (double)ins_ins / (ins - ins_end);
    a[2][3] = (double)ins_del / (ins - ins_end);
    a[2][4] = (double)ins_end / (ins - ins_end);

    a[3][0] = 0.0;
    a[3][1] = (double)del_match / (del - del_end);
    a[3][2] = (double)del_ins / (del - del_end);
    a[3][3] = (double)del_del / (del - del_end);
    a[3][4] = (double)del_end / (del - del_end);

    file.close();
}

int main(int argc, char *argv[])
{
    string input = argv[1];
    unordered_map<pair<string, string>, pair<string, string>> alignemnts;
    vector<string> states = {"Begin", "Match", "Insertion", "Deletion", "End"};
    vector<double> pi = {1.0, 0.0, 0.0, 0.0, 0.0};
    vector<vector<double>> a(states.size(), vector<double>(states.size(), 0.0));
    vector<unordered_map<string, double>> e(states.size(), unordered_map<string, double>());

    read_and_predict_hmm(input, alignemnts, a, e);

    for (auto &i : a)
    {
        for (auto &j : i)
        {
            cout << j << "\t";
        }
        cout << endl;
    }

    cout << "Total alignments read: " << alignemnts.size() << endl;

    return 0;
}