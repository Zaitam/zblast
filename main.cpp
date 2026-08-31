#include <assert.h>
#include <bits/stdc++.h>
#include <filesystem>
using namespace std;

#define forsn(i, s, n) for (int i = (s); i < (int)(n); i++)
#define forn(i, n) forsn(i, 0, n)
#define sz(x) (int)(x).size()
#define fi first
#define se second
#define pb push_back

using BLOSUM_WEIGHT = int8_t;
using BLOSUM = array<array<BLOSUM_WEIGHT, 25>, 25>;

BLOSUM load_blosum(const string &filename) {
    BLOSUM res;
    fstream f(filename.c_str(), ios::in);
    char c;
    BLOSUM_WEIGHT sign = 1; // 1 / -1
    forn(i, 25) forn(j, 25) {
        f >> c;
        if (c == '-') j--, sign = -1;
        else res[i][j] = (BLOSUM_WEIGHT)(c - '0') * sign, sign = 1;
    }
    return res;
}

struct LetterCode {
    int8_t value;
    LetterCode& operator=(const int8_t &new_value) {
        value = new_value;
        return *this;
    }
};

struct LetterCodeMap {
    array<LetterCode, 25> map;

    LetterCode& operator[](const char &letter_index) {
        if (letter_index == '*') return map[24]; // ultima pos
        int8_t index = (int8_t)(letter_index - 'A');
        if (index > (int8_t)('O'-'A')) index--;
        if (index > (int8_t)('U'-'A')) index--;
        return map[index];
    }
};

LetterCodeMap load_letter_code(const string &filename) {
    LetterCodeMap res{};
    fstream f(filename.c_str(), ios::in);
    char c;
    unsigned int t;
    while(f >> c) {
        f >> t;
        res[c] = (int8_t)t;
    }
    return res;
}

struct WordIndex {
    uint16_t value;
    WordIndex(uint16_t index) : value(index) {};
    void updateIndex(const uint8_t &new_last) { value = (value % (25 * 25)) * 25 + new_last; }
};
struct Word {
    array<uint8_t, 3> word; // Wordsize 3

    WordIndex toIndex() {
        return WordIndex(word[0] * 25 * 25 + word[1] * 25 + word[2]); // inyectiva al menos
    };
};
struct WordPermutation {
    array<vector<uint16_t>, 25 * 25 * 25> permutations;
    BLOSUM &blosum;

    WordPermutation(BLOSUM &blosum) : blosum(blosum) {};

    int transformCost(const Word& from, const Word& to) {
        int cost = 0;
        forn(i, sz(to.word)) cost += blosum[from.word[i]][to.word[i]];
        return cost;
    };

    void generate_permutations(int threshold) {
        for(uint8_t ii = 0; ii < 25; ++ii)
            for(uint8_t ji = 0; ji < 25; ++ji)
                for(uint8_t ki = 0; ki < 25; ++ki) {
                    Word initial = {ii, ji, ki};
                    for(uint8_t it = 0; it < 25; ++it)
                        for(uint8_t jt = 0; jt < 25; ++jt)
                            for(uint8_t kt = 0; kt < 25; ++kt) {
                                Word target = {it, jt, kt};
                                if (transformCost(initial, target) >= threshold)
                                    permutations[initial.toIndex().value].pb(target.toIndex().value);
                            }
                }
    }

    vector<uint16_t>& operator[](const uint16_t &value) {
        return permutations[value];
    }

    void save_permutations(const string& filename, int threshold) {
        fstream f(filename.c_str(), ios::out);
        f << threshold << '\n';
        forn(i, 25*25*25) {
            f << i << " ";
            f << sz(permutations[i]);
            forn(j, sz(permutations[i]))
                f << " " << permutations[i][j];
            f << '\n';
        }
    }
    void load_permutations(const string& filename, int threshold) {
        fstream f(filename.c_str(), ios::in);
        int t;
        f >> t;
        //assert(t == threshold); // same k
        if (t != threshold) return generate_and_save(filename, threshold);
        forn(i, 25*25*25) {
            f >> t;
            //assert(i == t); // broken cache
            if (i != t) return generate_and_save(filename, threshold);
            f >> t;
            permutations[i].resize(t);
            forn(j, t) f >> permutations[i][j];
        }
    }

    void generate_and_save(const string& filename, const int threshold) {
        cout << "START: BUILD WORD PERMUTATION CACHE" << '\n';
        generate_permutations(threshold);
        save_permutations(filename, threshold);
        cout << "DONE: BUILD WORD PERMUTATION CACHE" << '\n';
    }

    void auto_get_permutations(const string& filename, const int threshold) {
        if (filesystem::exists(filename)) load_permutations(filename, threshold);
        else generate_and_save(filename, threshold);
    }

};
//const static string DATABASE_REG = "./DB.reg.db";
using AccessNumber = uint16_t;
using Position = uint16_t;
using WordAlign = pair<AccessNumber, Position>;

struct DB {
    array<vector<WordAlign>, 25 * 25 * 25> wordReg;
    LetterCodeMap &letter_code_map;

    DB(LetterCodeMap & letter_map) : letter_code_map(letter_map) {};

    void add(const string &sequence, const AccessNumber access_number) {
        Word initial = {0,
            (uint8_t)letter_code_map[sequence[0]].value,
            (uint8_t)letter_code_map[sequence[1]].value};
        WordIndex index = initial.toIndex();
        forsn(i, 2, sz(sequence))
            index.updateIndex(letter_code_map[sequence[i]].value), wordReg[index.value].pb({access_number, i});
        // TODO: Offload to DISK
    };


    void load_file_to_database(const string& filename) {
        fstream f(filename.c_str(), ios::in);
        int access;
        string sequence;
        while(f >> access) {
            f >> sequence;
            add(sequence, (AccessNumber)access);
        }
    }

    void generate_appearence_map(Word word, map<AccessNumber, long long>& appearence_map) { //rename
        generate_appearence_map(word.toIndex().value, appearence_map);
    }
    void generate_appearence_map(uint16_t wordIndex, map<AccessNumber, long long> &appearence_map) {
        for (auto match: wordReg[wordIndex])
            appearence_map[match.fi]++;
    }
};

int main() {
    BLOSUM blosum = load_blosum("assets/BLOSUM62");
    LetterCodeMap index_map = load_letter_code("assets/LETTER_INDEX");

    DB database(index_map);
    database.load_file_to_database("assets/SEQUENCE_large");

    WordPermutation word_permutation(blosum);
    word_permutation.auto_get_permutations("permutation.cache", 11);
    string query_string;
    cin >> query_string;

    Word initial = {0,
            (uint8_t)index_map[query_string[0]].value,
            (uint8_t)index_map[query_string[1]].value};
    WordIndex index = initial.toIndex();
    map<AccessNumber, long long> appearence_map;
    forsn(i, 2, sz(query_string)) {
        index.updateIndex(index_map[query_string[i]].value);
        for ( auto wordIndex: word_permutation[index.value])
            database.generate_appearence_map(wordIndex, appearence_map);
    }
    // TODO: Al

    cout << 2;
    //for (auto i : blosum) for (auto j : i) cout << (int)j << ", ";
}