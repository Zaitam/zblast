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
    int c;
    forn(i, 25) forn(j, 25) {
        f >> c;
        res[i][j] = (BLOSUM_WEIGHT)(c);
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
        if (index > (int8_t)('O'-'A') - 1) index--;
        if (index > (int8_t)('U'-'A') - 1) index--;
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
    }

    uint16_t toValue() {
        return word[0] * 25 * 25 + word[1] * 25 + word[2]; // inyectiva al menos
    }
};

struct WordPermutation {
    array<vector<uint16_t>, 25 * 25 * 25> permutations;
    bitset<25 * 25 * 25> computed{};
    BLOSUM &blosum;
    int threshold;

    array<int8_t, 25> best; // best[x] = best core contra x
    WordPermutation(BLOSUM &blosum, int threshold) : blosum(blosum), threshold(threshold) {
        forn(x, 25) best[x] = *max_element(blosum[x].begin(), blosum[x].end());
    };

    int transformCost(const Word& from, const Word& to) {
        int cost = 0;
        forn(i, sz(to.word)) cost += blosum[from.word[i]][to.word[i]];
        return cost;
    };

    void dfs(const Words &from, int k, int partial, uint16_t index, vector<uint16_t> &out) {
        if (k == 3){ out.pb(index); return;}
        int max_rest = 0;
        forsn(j, k+1, 3) max_rest += best[from.word[j]];
        forn(c, 25) {
            int next = partial + blosum[from.word[k]][c];
            if (next + max_rest < threshold) continue;
            dfs(from, k + 1, next, (uint16_t)(index * 25 + c), out);
        }
    }

    void generate_all_permutations() {
        for(uint8_t ii = 0; ii < 25; ++ii)
            for(uint8_t ji = 0; ji < 25; ++ji)
                for(uint8_t ki = 0; ki < 25; ++ki) {
                    Word initial = {ii, ji, ki};
                    const uint16_t w = initial.toValue();
                    for(uint8_t it = 0; it < 25; ++it) {
                        for(uint8_t jt = 0; jt < 25; ++jt)
                            for(uint8_t kt = 0; kt < 25; ++kt) {
                                Word target = {it, jt, kt};
                                if (transformCost(initial, target) >= threshold)
                                    permutations[w].pb(target.toValue());
                            }
                    }
                    computed[w] = true;
                }
    }

    vector<uint16_t>& get_permutations(const uint16_t& from_index) {
        Word from = {(uint8_t)(from_index / 625), (uint8_t)(from_index / 25 % 25), (uint8_t)(from_index % 25)};
        return get_permutations(from);
    }
    vector<uint16_t>& get_permutations(Word& from) {
        const uint16_t w = from.toValue();
        if (computed[w]) return permutations[from.toValue()];
        for(uint8_t it = 0; it < 25; ++it)
            for(uint8_t jt = 0; jt < 25; ++jt)
                for(uint8_t kt = 0; kt < 25; ++kt) {
                    Word target = {it, jt, kt};
                    if (transformCost(from, target) >= threshold)
                        permutations[w].pb(target.toValue());
                }
        computed[w] = true;
        return permutations[w];
    }

    vector<uint16_t>& operator[](const uint16_t &value) {
        if (computed[value]) return permutations[value];
        Word from = {(uint8_t)(value / (25 * 25)), (uint8_t)(value / 25 % 25), (uint8_t)(value % 25)};
        dfs(from, 0, 0, 0, permutations[value]);
        computed[value] = true;
        return permutations[value];
        //return get_permutations(value);
    }

    void save_permutations(const string& filename) {
        fstream f(filename.c_str(), ios::out);
        f << threshold << '\n';
        forn(i, 25*25*25) {
            f << i << " ";
            if (computed[i]) {
                f << 1 << " " << sz(permutations[i]);
                forn(j, sz(permutations[i]))
                    f << " " << permutations[i][j];
            } else f << 0;
            f << '\n';
        }
    }
    void load_permutations(const string& filename) {
        fstream f(filename.c_str(), ios::in);
        int t;
        f >> t;
        assert(t == threshold); // same k
        for (auto& p : permutations) p.clear();
        forn(i, 25*25*25) {
            f >> t;
            assert(i == t); // broken cache
            f >> t;
            if (t) {
                f >> t;
                computed[i] = true;
                permutations[i].resize(t);
                forn(j, t) f >> permutations[i][j];
            } else {
                computed[i] = false;
                permutations[i].resize(0);
            }
        }
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
        generate_appearence_map(word.toValue(), appearence_map);
    }
    void generate_appearence_map(uint16_t wordIndex, map<AccessNumber, long long> &appearence_map) {
        for (auto match: wordReg[wordIndex])
            appearence_map[match.fi]++;
    }

};

using Diagonal = int32_t;
using DiagonalKey = uint32_t; // access << 16 | dbPos - queryPos + queryLen

const size_t INITIAL_CAPACITY = 1 << 16; // para los buffers
struct Hit {
    DiagonalKey key;
    Position query_pos;
    Position db_pos;

    Hit() = default;
    Hit(AccessNumber access, Position query_pos, Position db_pos, Position query_len) :  query_pos(query_pos), db_pos(db_pos) {
        const Diagonal diagonal = (Diagonal)db_pos - (Diagonal)query_pos;
        key = ((DiagonalKey)access << 16)
            | (DiagonalKey)(uint16_t)(diagonal + (Diagonal)query_len);
    }

    Diagonal diagonal() const { return (Diagonal)db_pos - (Diagonal)query_pos; }
    AccessNumber access() const { return (key >> 16); }
};

struct HitBuffer {
    Hit *hits = nullptr;
    size_t count = 0, capacity = 0;

    // Singletonish
    HitBuffer() : hits(new Hit[INITIAL_CAPACITY]), capacity(INITIAL_CAPACITY) {};
    HitBuffer(const HitBuffer&) = delete; // no borrar doble
    HitBuffer& operator=(const HitBuffer&) = delete; // no borrar doble
    ~HitBuffer() {delete[] hits;};


    void reserve(size_t n) {
        if (n > capacity) { delete[] hits; hits = new Hit[n]; capacity = n; }
        count = 0;
    }

    void grow() {
        Hit *bigger = new Hit[capacity * 2];
        forn(i, count) bigger[i] = hits[i];
        delete[] hits;
        hits = bigger, capacity *= 2;
    }

    void add(const Hit &hit) {
        if (count == capacity) grow();
        hits[count++] = hit;
    }

    Hit& operator[](const size_t &i) {return hits[i];};

    struct Radix {
        Hit *scratch = nullptr;
        size_t capacity = 0;

        Radix() = default;
        explicit Radix(size_t n) : scratch(new Hit[n]), capacity(n) {}
        ~Radix() { delete[] scratch; }

        Radix(const Radix&) = delete; // no borrar doble
        Radix& operator=(const Radix&) = delete; // no borrar doble

        void reserve(size_t n) {
            if (n > capacity) { delete[] scratch; scratch = new Hit[n]; capacity = n; }
        }

        void sort(HitBuffer &buffer) { // TODO: NO GET KEY DIRECTO KEY
            // cuantas pasadas (ahora es 4 por el seize pero es algo que en compile time podemos saber  y calcular dinamicamente...)
            using Key = decltype(Hit::key);
            constexpr int PASSES = (int)sizeof(Key);  // uint32 -> 4, uint64 -> 8
            if (buffer.count < 2) return; // abs ordenar con 2 o 1 elementos
            reserve(buffer.count);
            size_t histogram[PASSES][256] = {}; // 8kb -> 2kb => entra en L1

            for (size_t i = 0; i < buffer.count; ++i)
                forn(p, PASSES) ++histogram[p][(buffer.hits[i].key >> (p*8)) & 0xFF];

            Hit *from = buffer.hits, *to = scratch;
            forn(p, PASSES) {
                if (histogram[p][(from[0].key >> (p * 8)) & 0xFF] == buffer.count) continue; // skip si bucket = total
                size_t *hist = histogram[p], offset = 0;
                size_t n;
                forn(bucket, 256) n = hist[bucket], hist[bucket] = offset, offset += n;
                for (size_t i = 0; i < buffer.count; i++)
                    to[hist[(from[i].key >> (p * 8)) & 0xFF]++] = from[i];
                swap(from, to);
            }
            if (from == buffer.hits) return; // paridad -> cambiar el buffer a que quede en el otro array sino
            swap(buffer.hits, scratch);
            swap(buffer.capacity, capacity);
        }
    };
};

// basicamente un pair<uint32_t, uint32_t> pero se puede manejar por referencia o nullptr
struct Seed {uint32_t first, second;};


struct SeedBuffer {

    Seed *seeds = nullptr;
    size_t count = 0, capacity = 0;

    SeedBuffer() : seeds(new Seed[INITIAL_CAPACITY]), capacity(INITIAL_CAPACITY) {}; // como HitBuffer -> hacer generico?
    SeedBuffer(const SeedBuffer&) = delete; // no borrar doble
    SeedBuffer& operator=(const SeedBuffer&) = delete; // no borrar doble
    ~SeedBuffer() { delete[] seeds; }

    void grow() {
        Seed *bigger = new Seed[capacity * 2];
        forn(i, count) bigger[i] = seeds[i];
        delete[] seeds;
        seeds = bigger, capacity *= 2;
    }

    void add(const Seed &seed) {
        if (count == capacity) grow();
        seeds[count++] = seed;
    }

    Seed& operator[](const size_t &i) { return seeds[i]; }
};

void two_hit(HitBuffer &buffer, SeedBuffer &seeds, uint32_t window, uint32_t min_gap) {
    seeds.count = 0;
    size_t i = 0;
    while (i < buffer.count) {
        const DiagonalKey key = buffer[i].key;
        size_t last = i++;
        while (i < buffer.count && buffer[i].key == key) {
            const uint32_t gap = (uint32_t)(buffer[i].db_pos - buffer[last].db_pos);
            if (gap >= min_gap && gap <= window)
                seeds.add({(uint32_t)last, (uint32_t)i});
            last = i++;
        }
    }
}

#define el '\n'

int main() {
    BLOSUM blosum = load_blosum("assets/BLOSUM62");
    LetterCodeMap index_map = load_letter_code("assets/LETTER_INDEX");

    DB database(index_map);
    database.load_file_to_database("assets/SEQUENCE_large");

    WordPermutation word_permutation(blosum, 11);
    //word_permutation.auto_get_permutations("permutation.cache");
    string query_string;
    cin >> query_string;

    Word initial = {0,
            (uint8_t)index_map[query_string[0]].value,
            (uint8_t)index_map[query_string[1]].value};
    WordIndex index = initial.toIndex();
    //map<AccessNumber, long long> appearence_map;
    HitBuffy er buffer;
    const uint32_t qLen =  sz(query_string);
    forsn(i, 2, qLen) {
        index.updateIndex(index_map[query_string[i]].value);
        for ( auto wordIndex: word_permutation[index.value])
            for (auto &wordMatch : database.wordReg[wordIndex] )
                buffer.add(Hit(wordMatch.fi, (Position)i, wordMatch.se, qLen));
    }
    cout << "hits " << buffer.count << el;

    HitBuffer::Radix radix;
    radix.sort(buffer);

    SeedBuffer seeds;
    two_hit(buffer, seeds, 40, 3); // window, min_gap = wordsize
    cout << "seeds " << seeds.count << el;
    int extensions;
    forn(s, seeds.count)
        if (s == 0 || buffer[seeds[s].second].key != buffer[seeds[s-1].second].key) extensions++;
    cout << "ext " << extensions << el;
}