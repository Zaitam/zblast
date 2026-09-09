#include <bits/stdc++.h>
#include <filesystem>
#include <cassert>
using namespace std;

#define forsn(i, s, n) for (int i = (s); i < (int)(n); i++)
#define forn(i, n) forsn(i, 0, n)
#define sz(x) (int)(x).size()
#define fi first
#define se second
#define pb push_back

const uint32_t WORDSIZE = 3;
using BLOSUM_WEIGHT = int8_t;
using BLOSUM = array<array<BLOSUM_WEIGHT, 25>, 25>;

BLOSUM load_blosum(const string &filename) {
    BLOSUM res{};
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

    void dfs(const Word &from, int k, int partial, uint16_t index, vector<uint16_t> &out) {
        if (k == 3){ out.pb(index); return;}
        int max_rest = 0;
        forsn(j, k+1, 3) max_rest += best[from.word[j]];
        forn(c, 25) {
            int next = partial + blosum[from.word[k]][c];
            if (next + max_rest < threshold) continue;
            dfs(from, k + 1, next, (uint16_t)(index * 25 + c), out);
        }
    }

    vector<uint16_t>& operator[](const uint16_t &value) {
        if (computed[value]) return permutations[value];
        Word from = {(uint8_t)(value / (25 * 25)), (uint8_t)(value / 25 % 25), (uint8_t)(value % 25)};
        dfs(from, 0, 0, 0, permutations[value]);
        computed[value] = true;
        return permutations[value];
    }

};
//const static string DATABASE_REG = "./DB.reg.db";
using AccessNumber = uint16_t;
using Position = uint16_t;
using WordAlign = pair<AccessNumber, Position>;

struct DB {
    array<vector<WordAlign>, 25 * 25 * 25> wordReg;
    LetterCodeMap &letter_code_map;

    vector<uint8_t>  residues;
    vector<uint32_t> offsets; // offsets[a] .. offsets[a+1] es la secuencia a
    Position max_sequence_len = 0;

    uint8_t* sequence_of(AccessNumber &access) {return residues.data() + offsets[access]; };
    uint32_t length_of(AccessNumber &access) { return offsets[access + 1] - offsets[access]; };
    size_t sequence_count() { return offsets.empty() ? 0 : offsets.size() - 1; };


    DB(LetterCodeMap & letter_map) : letter_code_map(letter_map) {};

    void add(const string &sequence, const AccessNumber access_number) {
        Word initial = {0,
            (uint8_t)letter_code_map[sequence[0]].value,
            (uint8_t)letter_code_map[sequence[1]].value};
        WordIndex index = initial.toIndex();

        if (offsets.empty()) offsets.pb(0);
        forn(i, sz(sequence)) residues.pb((uint8_t) letter_code_map[sequence[i]].value); // store sequence code
        offsets.pb(residues.size());
        max_sequence_len = max(max_sequence_len, (Position)sz(sequence));
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
            if (gap < min_gap) {i++; continue;} // Altschul 1997 p.3391: "Any hit that overlaps the most recent one is ignored."
            if (gap <= window) seeds.add({(uint32_t)last, (uint32_t)i});
            last = i++;
        }
    }
}

struct Extension {
    Position db_start, db_end;
    Position q_start, q_end;
    int32_t score;
    AccessNumber access;
};

struct Extender {
    BLOSUM &blosum;
    uint8_t *query;
    uint32_t query_len; // Variable por query
    uint8_t x_drop_ungapped, x_drop_gapped, gap_open, gap_extend;
    Extender(BLOSUM& blosum, uint8_t *query, uint32_t query_len, uint8_t x_drop_ungapped = 20, uint8_t x_drop_gapped = 40, uint8_t gap_open = 10, uint8_t gap_extend = 1) :
        blosum(blosum), query(query), query_len(query_len), x_drop_ungapped(x_drop_ungapped), x_drop_gapped(x_drop_gapped), gap_open(gap_open), gap_extend(gap_extend) {};

    // A partir del Seed -> Diagonal cte
    Extension ungapped(uint8_t *db, uint32_t db_len, uint32_t qp, uint32_t dp) {
        uint32_t back = min(WORDSIZE - 1, min(qp, dp));
        uint32_t q0 = qp - back;
        uint32_t d0 = dp - back;

        int32_t ref = 0;
        forn(k, back + 1) ref += blosum[query[q0 + k]][db[d0 + k]];

        // Stretch right
        int32_t score = 0;
        int32_t right_best = 0;
        uint32_t right_len = 0;
        for (uint32_t k = 1; qp + k < query_len && dp + k < db_len; ++k) {
            score += blosum[query[qp + k]][db[dp + k]];
            if (score > right_best) right_best = score, right_len = k;
            if (right_best - score > x_drop_ungapped) break;
        }
        // Stretch left
        score = 0;
        int32_t left_best = 0;
        uint32_t left_reach = 0;
        for (uint32_t k = 1; k <= q0 && k <= d0; ++k) {
            score += blosum[query[q0 - k]][db[d0 - k]];
            if (score > left_best) left_best = score, left_reach = k;
            if (left_best - score > x_drop_ungapped) break;
        }
        return {(Position)(d0 - left_reach), (Position)(dp + right_len + 1),
                (Position)(q0 - left_reach), (Position)(qp + right_len + 1),
                left_best + ref + right_best, 0};
    }

    static constexpr int32_t NEG = INT32_MIN / 4;
    int32_t half_gapped(uint8_t *q, int32_t qlen, uint8_t *d, int32_t dlen, int32_t &q_reach, int32_t &d_reach) {
        q_reach = d_reach = 0;
        if (qlen <= 0 || dlen <= 0) return 0;

        vector<int> Hprev(dlen + 1, NEG), Hcur(dlen + 1, NEG), Ecur(dlen + 1, NEG);
        int best = 0, lo = 0, hi = 0;

        Hprev[0] = 0; // fila 0 -> gaps
        for (int j = 1; j <= dlen; ++j) {
            int v = -(gap_open + j * gap_extend);
            if (best - v > x_drop_gapped) break;
            Hprev[j] = v;
            hi = j;
        }
        for (int i = 1; i <= qlen; ++i) {
            fill(Hcur.begin(), Hcur.end(), NEG);
            fill(Ecur.begin(), Ecur.end(), NEG);
            int f = NEG; // gap en db
            int new_lo = -1, new_hi = -1;
            const int jend = min(hi + 1, dlen);
            for (int j = lo; j <= jend; ++j) {
                int m = NEG; // diagonal: match/mismatch
                if (j >= 1 && Hprev[j-1] > NEG/2) m = Hprev[j-1] + blosum[q[i-1]][d[j-1]];
                int e = NEG; // gap en query
                if (j >= 1) {
                    int fh = (Hcur[j-1] > NEG/2) ? Hcur[j-1] - gap_open - gap_extend : NEG;
                    int fe = (Ecur[j-1] > NEG/2) ? Ecur[j-1] - gap_extend            : NEG;
                    e = max(fh, fe);
                }
                int fh = (Hprev[j] > NEG/2) ? Hprev[j] - gap_open - gap_extend : NEG;
                int ff = (f        > NEG/2) ? f        - gap_extend            : NEG;
                f = max(fh, ff);

                int h = max(m, max(e, f));
                if (h <= NEG/2 || best - h > x_drop_gapped) continue;
                Hcur[j] = h; Ecur[j] = e;
                if (h > best) best = h, q_reach = i, d_reach = j;
                if (new_lo < 0) new_lo = j;
                new_hi = j;
            }
            if (new_lo < 0) break;
            lo = new_lo; hi = new_hi;
            swap(Hprev, Hcur);
        }
        return best;
    }

    Extension gapped(uint8_t *db, uint32_t db_len, uint32_t qp, uint32_t dp) {
        int anchor = blosum[query[qp]][db[dp]];
        int rq, rd, lq, ld;
        int right = half_gapped(query + qp + 1, (int32_t)query_len - (int32_t)qp - 1,
                                (db + dp + 1), (int32_t)db_len - (int32_t)dp - 1, rq, rd);
        vector<uint8_t> qr(query, query + qp), dr(db, db + dp);
        reverse(qr.begin(), qr.end()); reverse(dr.begin(), dr.end());
        int left = half_gapped(qr.data(), (int32_t)qr.size(), dr.data(), (int32_t)dr.size(), lq, ld);
        return {(Position)(dp - ld), (Position)(dp + rd + 1),
                (Position)(qp - lq), (Position)(qp + rq + 1),
                anchor + left + right, 0};
    }
};

array<double, 25> load_background_freq(const string &filename, LetterCodeMap &map) {
    array<double, 25> res{};
    fstream f(filename.c_str(), ios::in);
    char c;
    double p;
    double total = 0;
    while(f >> c) {
        f >> p;
        res[map[c].value] = p, total += p;
    }
    return res;
}

struct GappedConfig { double lambda, K, H_bits; };
GappedConfig load_karlin_gapped(const string &filename, int gap_open, int gap_extend) {
    fstream f(filename.c_str(), ios::in);
    int go, ge;
    double lambda, k, h;
    while(f >> go) {
        f >> ge >> lambda >> k >> h;
        if (go == gap_open && ge == gap_extend) return {lambda, k, h};
    }
    throw runtime_error(filename + ": no hay fila para gap_open=" + to_string(gap_open) +
                        " gap_extend=" + to_string(gap_extend));
}

struct KarlinAltschul {
    array<double, 25> P{};
    double long lambda, K, H;

    KarlinAltschul(BLOSUM &blosum, bool with_gaps, array<double, 25> &fallback, GappedConfig &gapped, vector<uint8_t> *residues = nullptr) {
        double total = 0;
        if (residues) for (uint8_t r : *residues) if (r < 20) P[r] += 1, total += 1;

        if (total > 0) forn(i, 20) P[i] /= total;
        else forn(i, 20) P[i] = fallback[i];

        if (with_gaps) {
            lambda = gapped.lambda; K = gapped.K;
            H = gapped.H_bits * log(2.0); // la tabla esta en bits, aca van nats
            return;
        }
        // biseccion sobre F(l) = sum p_i p_j e^(l s_ij) - 1
        // F(0)=0, F'(0)<0 (score esperado al azar es negativo) -> 1 sola raiz positiva
        auto F = [&](double l) {
            double s = 0;
            forn(i, 20) forn(j, 20) s += P[i] * P[j] * exp(l * blosum[i][j]);
            return s - 1.0;
        };
        double lo = 1e-9, hi = 2.0;
        forn(it, 200) { double mid = (lo + hi) / 2; (F(mid) < 0 ? lo : hi) = mid; }
        lambda = (lo + hi) / 2;
        H = 0;
        forn(i, 20) forn(j, 20) H += P[i] * P[j] * exp(lambda * blosum[i][j]) * blosum[i][j];
        H *= lambda;
        K = gapped.K;
    }

    double bits(int score) const { return (lambda * score - log(K)) / log(2.0); }

    double evalue(int score, uint32_t m, double n, uint32_t num_seqs) const {
        // Mal borde -> alineamiento en ultimo residuio -> espacio menor que m*n -> iterar (en 5 converge) [Copiado de NCBI]
        double L = 0, me = m, ne = n;
        forn(it, 5) {
            me = max((double)(1.0 / K), (double)m - L);
            ne = max((double)(1.0 / K), n - (double)num_seqs * L);
            L = log(K * me * ne) / H;
            if (L < 0) L = 0;
        }
        return K * me * ne * exp(-lambda * score);
    }
};

#define el '\n'

int main() try {
    BLOSUM blosum = load_blosum("assets/BLOSUM62");
    LetterCodeMap index_map = load_letter_code("assets/LETTER_INDEX");
    array<double, 25> background = load_background_freq("assets/BACKGROUND_FREQ", index_map);
    GappedConfig ungapped_config = load_karlin_gapped("assets/KARLIN_GAPPED", 0, 0);
    GappedConfig gapped_config   = load_karlin_gapped("assets/KARLIN_GAPPED", 10, 1);


    DB database(index_map);
    database.load_file_to_database("assets/SEQUENCE_large");

    WordPermutation word_permutation(blosum, 11);
    //word_permutation.auto_get_permutations("permutation.cache");
    string query_string;
    cin >> query_string;

    uint32_t qLen = sz(query_string);
    vector<uint8_t> query;
    forn(i, qLen) query.pb((uint8_t)index_map[query_string[i]].value);

    Word initial = {0, query[0], query[1]};
    WordIndex index = initial.toIndex();
    HitBuffer buffer;
    forsn(i, 2, qLen) {
        index.updateIndex(query[i]);
        for ( auto wordIndex: word_permutation[index.value])
            for (auto &wordMatch : database.wordReg[wordIndex] )
                buffer.add(Hit(wordMatch.fi, (Position)i, wordMatch.se, qLen));
    }
    cout << "hits " << buffer.count << el;

    HitBuffer::Radix radix;
    radix.sort(buffer);

    SeedBuffer seeds;
    two_hit(buffer, seeds, 40, WORDSIZE); // window, min_gap = wordsize
    cout << "seeds " << seeds.count << el;

    Extender extender(blosum, query.data(), qLen);
    KarlinAltschul ku(blosum, false, background, ungapped_config, &database.residues);
    KarlinAltschul kg(blosum, true,  background, gapped_config,   &database.residues);

    vector<Extension> hsps;
    DiagonalKey covered_key = 0;
    Position covered_until = 0;   // arranca en 0: db_pos < 0 nunca es cierto sin signo
    forn(s, seeds.count) {
        Hit &hit = buffer[seeds[s].second];
        if (hit.key == covered_key && hit.db_pos < covered_until) continue;
        AccessNumber acc = hit.access();
        Extension e = extender.ungapped(database.sequence_of(acc), database.length_of(acc),
                                        hit.query_pos, hit.db_pos);
        e.access = acc;
        covered_key = hit.key, covered_until = e.db_end;
        if (ku.bits(e.score) >= 22) hsps.pb(e); // S_g del paper: ~22 bits
    }
    cout << "hsps " << sz(hsps) << el;

    for (auto &h : hsps) {
        Extension g = extender.gapped(database.sequence_of(h.access), database.length_of(h.access),
                                      (h.q_start + h.q_end) / 2, (h.db_start + h.db_end) / 2);
        g.access = h.access;
        if (g.score > h.score) h = g;
    }

    map<AccessNumber, Extension> best;
    for (auto &h : hsps) if (!best.count(h.access) || h.score > best[h.access].score) best[h.access] = h;
    vector<Extension> results;
    for (auto &kv : best) results.pb(kv.second);
    sort(results.begin(), results.end(), [](Extension &a,Extension &b) { return a.score > b.score; });

    double db_len = database.residues.size();
    uint32_t num_seqs = database.sequence_count();
    printf("\n%-8s %7s %7s %11s   %-13s %s\n", "access", "score", "bits", "E-value", "query", "base");
    int shown = 0;
    for (auto &h : results) {
        double e = kg.evalue(h.score, qLen, db_len, num_seqs);
        if (e > 0.001) break;
        if (++shown > 20) break;
        printf("%-8u %7d %7.1f %11.2g   %5u-%-7u %6u-%-8u\n", h.access, h.score,
               (double)kg.bits(h.score), e, h.q_start, h.q_end, h.db_start, h.db_end);
    }
    if (!shown) forn(i, min(5, sz(results))) {
        Extension &h = results[i];
        printf("%-8u %7d %7.1f %11.2g   %5u-%-7u %6u-%-8u  NS\n", h.access, h.score,
               (double)kg.bits(h.score), kg.evalue(h.score, qLen, db_len, num_seqs),
               h.q_start, h.q_end, h.db_start, h.db_end);
    }
    cout << el << "lambda " << (double)ku.lambda << "  K " << (double)ku.K << "  H " << (double)ku.H << " (from DB)" << el;
    return 0;
}
catch (exception &e) { cerr << "error: " << e.what() << el; return 1; }
