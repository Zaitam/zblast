#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <filesystem>
#include <cassert>

using namespace std;

#define forsn(i, s, n) for (int i = (s); i < (int) (n); i++)
#define forn(i, n) forsn(i, 0, n)
#define sz(x) (int) (x).size()
#define fi first
#define se second
#define pb push_back

const uint32_t WORDSIZE = 3;

// Source - https://stackoverflow.com/a/27270738
// Posted by ikh
// Retrieved 2026-09-09, License - CC BY-SA 3.0
template<int A, int B>
struct get_power {
    static const uint64_t value = A * get_power<A, B - 1>::value;
};

template<int A>
struct get_power<A, 0> {
    static const uint64_t value = 1;
};

const size_t WORDS = get_power<25, WORDSIZE>::value;
const size_t WORD_STRIDE = get_power<25, WORDSIZE - 1>::value; // Valor primera letra -> la tiro MOD y meter nueva
static_assert(WORDSIZE <= 6, "WORDSIZE demasiado grande: Se me rompe el indice");

using BLOSUM_WEIGHT = int8_t;
using BLOSUM = array<array<BLOSUM_WEIGHT, 25>, 25>;

// 0 .. 25^WORDSIZE -> uint16 o uint32 ... si muy chico salta el assert
using WordValue = conditional_t<(WORDS <= (size_t) numeric_limits<uint16_t>::max() + 1), uint16_t, uint32_t>;
static_assert(WORDS <= (size_t) numeric_limits<WordValue>::max() + 1, "WordValue no entra en WORDSIZE");

using AccessNumber = uint32_t;
using Position = uint32_t;
using Length = uint32_t;
using Score = int32_t; // un score de alineamiento
using StoredPosition = uint16_t;
using Diagonal = int32_t;
using DiagonalKey = uint64_t; // access << 32 | dbPos - queryPos + queryLen

const int WORD_THRESHOLD = 11; // T del vecindario
const Length TWO_HIT_WINDOW = 40; // A de Altschul 1997
const Score HSP_MIN_BITS = 22; //S_g para pasar al gapped
const Score X_DROP_UNGAPPED = 20, X_DROP_GAPPED = 40;
const Score GAP_OPEN = 11, GAP_EXTEND = 1; // existence 11, extension 1
const double MAX_EVALUE = 0.001; // Limite
const int MAX_RESULTS = 20;

BLOSUM load_blosum(const string &filename) {
    BLOSUM res{};
    fstream f(filename.c_str(), ios::in);
    int c;
    forn(i, 25)
        forn(j, 25) {
            f >> c;
            res[i][j] = (BLOSUM_WEIGHT) (c);
        }
    return res;
}

struct LetterCode {
    int8_t value;

    LetterCode &operator=(const int8_t &new_value) {
        value = new_value;
        return *this;
    }
};

struct LetterCodeMap {
    array<LetterCode, 256> map; // char -> codigo
    array<char, 25> letter; // codigo -> char

    LetterCodeMap() {
        forn(i, 256) map[i] = (int8_t) -1;
        letter.fill('X');
    }

    LetterCode &operator[](const char &letter_index) { return map[(uint8_t) letter_index]; }
};

LetterCodeMap load_letter_code(const string &filename) {
    LetterCodeMap res{};
    fstream f(filename.c_str(), ios::in);
    char c;
    unsigned int t;
    while (f >> c) {
        f >> t;
        res[c] = (int8_t) t;
        res.letter[t] = c;
        if (isalpha((uint8_t) c)) res[(char) (c | 32)] = (int8_t) t; // la misma en minuscula
    }
    const int8_t unknown = res['X'].value; // digits o simbolos o U o O
    forn(i, 256)
        if (res.map[i].value < 0) res.map[i] = unknown;
    return res;
}

struct WordIndex {
    WordValue value;

    WordIndex(WordValue index) : value(index) {};
    void updateIndex(const uint8_t &new_last) { value = (value % WORD_STRIDE) * 25 + new_last; }
};

struct Word {
    array<uint8_t, WORDSIZE> word;

    WordValue toValue() {
        WordValue value = 0;
        forn(i, WORDSIZE) value = value * 25 + word[i]; // inyectiva al menos
        return value;
    }

    WordIndex toIndex() { return WordIndex(toValue()); }
};
template<class F> // Recorre palabras de una secuencia -> ventana corrediza -> llama a body
void for_each_word(const uint8_t *sequence, Length len, F body) {
    if (len < WORDSIZE) return;
    WordValue prefix = 0;
    forn(i, WORDSIZE - 1) prefix = (WordValue) (prefix * 25 + sequence[i]);
    WordIndex word = prefix;
    forsn(last, WORDSIZE - 1, len) {
        word.updateIndex(sequence[last]);
        body((Position) last, word.value);
    }
}

struct WordPermutation {
    // WordValue -> [permutations]; Si WORDSIZE es grande, se toca un % muy bajo -> solo store los que usamos
    // TODO: Hacer un custom hash para WordValue -> vector<WordValue>
    // TODO: Cache en disco
    unordered_map<WordValue, vector<WordValue>> permutations;
    BLOSUM &blosum;
    int threshold;

    array<int8_t, 25> best; // best[x] = best core contra x
    WordPermutation(BLOSUM &blosum, int threshold) : blosum(blosum), threshold(threshold) {
        forn(x, 25) best[x] = *max_element(blosum[x].begin(), blosum[x].end());
    };

    void dfs(const Word &from, int k, int partial, WordValue index, vector<WordValue> &out) {
        if (k == (int) WORDSIZE) {
            out.pb(index);
            return;
        }
        int max_rest = 0;
        forsn(j, k + 1, WORDSIZE) max_rest += best[from.word[j]];
        forn(c, 25) {
            int next = partial + blosum[from.word[k]][c];
            if (next + max_rest < threshold) continue;
            dfs(from, k + 1, next, (WordValue) (index * 25 + c), out);
        }
    }

    vector<WordValue> &operator[](const WordValue &value) {
        auto found = permutations.find(value);
        if (found != permutations.end()) return found->second;
        Word from{};
        WordValue rest = value;
        for (int i = (int) WORDSIZE - 1; i >= 0; i--) from.word[i] = (uint8_t) (rest % 25), rest /= 25;
        vector<WordValue> &out = permutations[value];
        dfs(from, 0, 0, 0, out);
        return out;
    }
};

static const uint64_t ZBX_MAGIC = 0x3258424C415A55ULL; // "UZALBX2"
static const size_t ZBX_HEADER = 48; // 6 x uint64, 0 (mod 8)

struct DB {
    LetterCodeMap &letter_code_map;

    const uint8_t *mapped = nullptr; // Const pq PROT_READ -> SIGBUS
    size_t mapped_len = 0;
    const uint32_t *word_off = nullptr;
    const AccessNumber *post_access = nullptr;
    const StoredPosition *post_pos = nullptr;
    const uint32_t *offsets = nullptr; // offsets[a] .. offsets[a+1] es la secuencia a
    const uint8_t *residues = nullptr;
    uint64_t nres = 0, nseq = 0, npost = 0;

    Length max_sequence_len = 0;

    DB(LetterCodeMap &letter_map) : letter_code_map(letter_map) {};

    DB(const DB &) = delete; // SIngletonsih
    DB &operator=(const DB &) = delete;

    ~DB() {
        if (mapped) munmap((void *) mapped, mapped_len);
    }

    const uint8_t *sequence_of(const AccessNumber &access) const { return residues + offsets[access]; };
    Length length_of(const AccessNumber &access) const { return offsets[access + 1] - offsets[access]; };
    size_t sequence_count() const { return nseq ? nseq - 1 : 0; };

    struct Appearence {
        const AccessNumber *access;
        const StoredPosition *pos;
        uint32_t n;
        uint32_t size() const { return n; }
    };

    Appearence operator[](const WordValue &word) const {
        const uint32_t from = word_off[word];
        return {post_access + from, post_pos + from, word_off[word + 1] - from};
    }
    // .txt -> CSR (Compressed Sparse Row) en disco
    void build_index(const string &source, const string &out) {
        const filesystem::path parent = filesystem::path(out).parent_path();
        if (!parent.empty()) filesystem::create_directories(parent);
        vector<uint8_t> res;
        vector<uint32_t> off;
        Length maxlen = 0;
        {
            fstream f(source.c_str(), ios::in);
            int access;
            string sequence;
            off.pb(0);
            while (f >> access) {
                f >> sequence;
                forn(i, sz(sequence)) res.pb((uint8_t) letter_code_map[sequence[i]].value);
                off.pb(res.size());
                maxlen = max(maxlen, (Length) sz(sequence));
            }
        }

        // contar por palabra y prefix sum -> word_start queda con limits
        vector<uint32_t> word_start(WORDS + 1, 0);
        uint64_t total = 0;
        forn(a, sz(off) - 1)
            for_each_word(res.data() + off[a], off[a + 1] - off[a],
                          [&](Position, WordValue word) { word_start[word + 1]++, total++; });
        forn(w, WORDS) word_start[w + 1] += word_start[w];

        // write -> cursor[w] = next hole de w
        vector<AccessNumber> acc(total);
        vector<StoredPosition> pos(total);
        vector<uint32_t> cursor(word_start.begin(), word_start.end() - 1);
        forn(a, sz(off) - 1)
            for_each_word(res.data() + off[a], off[a + 1] - off[a], [&](Position last, WordValue word) {
                uint32_t at = cursor[word]++;
                acc[at] = (AccessNumber) a, pos[at] = (StoredPosition) last;
            });

        // Sort write por size -> array aligned -> primero 4bytes, 2bytes, 1byte
        // [Antes, mal ordenado, causaba que leer un uint32 este mal alineado...]
        fstream f(out.c_str(), ios::out | ios::binary);
        uint64_t header[6] = {ZBX_MAGIC, res.size(), off.size(), total, (uint64_t) maxlen, 0};
        f.write((char *) header, ZBX_HEADER);
        f.write((char *) word_start.data(), (WORDS + 1) * sizeof(uint32_t));
        f.write((char *) acc.data(), total * sizeof(AccessNumber));
        f.write((char *) off.data(), off.size() * sizeof(uint32_t));
        f.write((char *) pos.data(), total * sizeof(StoredPosition));
        f.write((char *) res.data(), res.size());
    }

    // mmap -> SO asocia archivo a rango de dirs -> el SO trae paginas del disco (page fault)
    // offload al disco automatico por el SO
    // corrimiento at += nres para cada pos
    bool map_index(const string &filename) {
        if (mapped) munmap((void *) mapped, mapped_len), mapped = nullptr;
        int fd = ::open(filename.c_str(), O_RDONLY);
        if (fd < 0) return false;
        struct stat st{};
        fstat(fd, &st);
        mapped_len = (size_t) st.st_size;
        void *p = mmap(nullptr, mapped_len, PROT_READ, MAP_PRIVATE, fd, 0);
        close(fd); // el mapeo sobrevive al fd cerrado
        if (p == MAP_FAILED) return false;
        mapped = (const uint8_t *) p;
        const uint64_t *header = (const uint64_t *) mapped;
        if (mapped_len < ZBX_HEADER || header[0] != ZBX_MAGIC) return false;
        nres = header[1], nseq = header[2], npost = header[3];
        max_sequence_len = (Length) header[4];

        size_t at = ZBX_HEADER;
        word_off = (const uint32_t *) (mapped + at), at += (WORDS + 1) * sizeof(uint32_t);
        post_access = (const AccessNumber *) (mapped + at), at += npost * sizeof(AccessNumber);
        offsets = (const uint32_t *) (mapped + at), at += nseq * sizeof(uint32_t);
        post_pos = (const StoredPosition *) (mapped + at), at += npost * sizeof(StoredPosition);
        residues = mapped + at, at += nres;
        return at == mapped_len;
    }

    bool open(const string &source, const string &index) {
        struct stat src{}, idx{};
        bool fresh = stat(source.c_str(), &src) == 0 && stat(index.c_str(), &idx) == 0 && idx.st_mtime >= src.st_mtime;
        if (fresh && map_index(index)) return true;
        build_index(source, index);
        return map_index(index);
    }
};

static const int ACCESS_KEY_SHIFT = 32;

const size_t INITIAL_CAPACITY = 1 << 16; // para los buffers
struct Hit {
    DiagonalKey key;
    Position query_pos;
    Position db_pos;

    Hit() = default;

    Hit(AccessNumber access, Position query_pos, Position db_pos, Length query_len) :
        query_pos(query_pos), db_pos(db_pos) {
        const Diagonal diagonal = (Diagonal) db_pos - (Diagonal) query_pos;
        key = ((DiagonalKey) access << ACCESS_KEY_SHIFT) | (DiagonalKey) (uint16_t) (diagonal + (Diagonal) query_len);
    }

    Diagonal diagonal() const { return (Diagonal) db_pos - (Diagonal) query_pos; }
    AccessNumber access() const { return (key >> ACCESS_KEY_SHIFT); }
};

struct HitBuffer {
    Hit *hits = nullptr;
    size_t count = 0, capacity = 0;

    // Singletonish
    HitBuffer() : hits(new Hit[INITIAL_CAPACITY]), capacity(INITIAL_CAPACITY) {};

    HitBuffer(const HitBuffer &) = delete; // no borrar doble
    HitBuffer &operator=(const HitBuffer &) = delete; // no borrar doble
    ~HitBuffer() { delete[] hits; };

    void reserve(size_t n) {
        if (n > capacity) {
            delete[] hits;
            hits = new Hit[n];
            capacity = n;
        }
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

    Hit &operator[](const size_t &i) { return hits[i]; };

    struct Radix {
        Hit *scratch = nullptr;
        size_t capacity = 0;

        Radix() = default;

        explicit Radix(size_t n) : scratch(new Hit[n]), capacity(n) {}

        ~Radix() { delete[] scratch; }

        Radix(const Radix &) = delete; // no borrar doble
        Radix &operator=(const Radix &) = delete; // no borrar doble

        void reserve(size_t n) {
            if (n > capacity) {
                delete[] scratch;
                scratch = new Hit[n];
                capacity = n;
            }
        }

        void sort(HitBuffer &buffer) {
            // TODO: NO GET KEY DIRECTO KEY
            // cuantas pasadas (ahora es 4 por el seize pero es algo que en compile time podemos saber  y calcular dinamicamente...)
            using Key = decltype(Hit::key);
            constexpr int PASSES = (int) sizeof(Key); // uint32 -> 4, uint64 -> 8
            if (buffer.count < 2) return; // abs ordenar con 2 o 1 elementos
            reserve(buffer.count);
            size_t histogram[PASSES][256] = {}; // 8kb -> 2kb => entra en L1

            for (size_t i = 0; i < buffer.count; ++i)
                forn(p, PASSES) ++histogram[p][(buffer.hits[i].key >> (p * 8)) & 0xFF];

            Hit *from = buffer.hits, *to = scratch;
            forn(p, PASSES) {
                if (histogram[p][(from[0].key >> (p * 8)) & 0xFF] == buffer.count) continue; // skip si bucket = total
                size_t *hist = histogram[p], offset = 0;
                size_t n;
                forn(bucket, 256) n = hist[bucket], hist[bucket] = offset, offset += n;
                for (size_t i = 0; i < buffer.count; i++) to[hist[(from[i].key >> (p * 8)) & 0xFF]++] = from[i];
                swap(from, to);
            }
            if (from == buffer.hits) return; // paridad -> cambiar el buffer a que quede en el otro array sino
            swap(buffer.hits, scratch);
            swap(buffer.capacity, capacity);
        }
    };
};

// basicamente un pair<uint32_t, uint32_t> pero se puede manejar por referencia o nullptr
struct Seed {
    uint32_t first, second;
};

struct SeedBuffer {
    Seed *seeds = nullptr;
    size_t count = 0, capacity = 0;

    SeedBuffer() :
        seeds(new Seed[INITIAL_CAPACITY]), capacity(INITIAL_CAPACITY) {}; // como HitBuffer -> hacer generico?
    SeedBuffer(const SeedBuffer &) = delete; // no borrar doble
    SeedBuffer &operator=(const SeedBuffer &) = delete; // no borrar doble
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

    Seed &operator[](const size_t &i) { return seeds[i]; }
};

void two_hit(HitBuffer &buffer, SeedBuffer &seeds, uint32_t window, uint32_t min_gap) {
    seeds.count = 0;
    size_t i = 0;
    while (i < buffer.count) {
        const DiagonalKey key = buffer[i].key;
        size_t last = i++;
        while (i < buffer.count && buffer[i].key == key) {
            const uint32_t gap = (uint32_t) (buffer[i].db_pos - buffer[last].db_pos);
            if (gap < min_gap) {
                i++;
                continue;
            } // Altschul 1997 p.3391: "Any hit that overlaps the most recent one is ignored."
            if (gap <= window) seeds.add({(uint32_t) last, (uint32_t) i});
            last = i++;
        }
    }
}

struct Extension {
    Position db_start, db_end;
    Position q_start, q_end;
    Score score;
    AccessNumber access;
};

struct Extender {
    BLOSUM &blosum;
    uint8_t *query;
    Length query_len; // Variable por query
    Score x_drop_ungapped, x_drop_gapped, gap_open, gap_extend;

    Extender(BLOSUM &blosum, uint8_t *query, Length query_len, Score x_drop_ungapped, Score x_drop_gapped,
             Score gap_open, Score gap_extend) :
        blosum(blosum), query(query), query_len(query_len), x_drop_ungapped(x_drop_ungapped),
        x_drop_gapped(x_drop_gapped), gap_open(gap_open), gap_extend(gap_extend) {};

    // A partir del Seed -> Diagonal cte
    Extension ungapped(const uint8_t *db, Length db_len, Position qp, Position dp) {
        Length back = min(WORDSIZE - 1, min(qp, dp));
        Position q0 = qp - back;
        Position d0 = dp - back;

        Score ref = 0;
        forn(k, back + 1) ref += blosum[query[q0 + k]][db[d0 + k]];

        // Stretch right
        Score score = 0;
        Score right_best = 0;
        Length right_len = 0;
        for (Length k = 1; qp + k < query_len && dp + k < db_len; ++k) {
            score += blosum[query[qp + k]][db[dp + k]];
            if (score > right_best) right_best = score, right_len = k;
            if (right_best - score > x_drop_ungapped) break;
        }
        // Stretch left
        score = 0;
        Score left_best = 0;
        Length left_reach = 0;
        for (Length k = 1; k <= q0 && k <= d0; ++k) {
            score += blosum[query[q0 - k]][db[d0 - k]];
            if (score > left_best) left_best = score, left_reach = k;
            if (left_best - score > x_drop_ungapped) break;
        }
        return {(Position) (d0 - left_reach),    (Position) (dp + right_len + 1), (Position) (q0 - left_reach),
                (Position) (qp + right_len + 1), left_best + ref + right_best,    0};
    }

    static constexpr Score NEG = numeric_limits<Score>::min();

    Score half_gapped(const uint8_t *q, Length qlen, const uint8_t *d, Length dlen, Length &q_reach, Length &d_reach) {
        q_reach = d_reach = 0;
        if (qlen == 0 || dlen == 0) return 0;
        const int32_t qn = (int32_t) qlen, dn = (int32_t) dlen;

        vector<Score> Hprev(dn + 1, NEG), Hcur(dn + 1, NEG), Ecur(dn + 1, NEG);
        Score best = 0;
        int lo = 0, hi = 0;

        Hprev[0] = 0; // fila 0 -> gaps
        for (int j = 1; j <= dn; ++j) {
            Score v = -(gap_open + j * gap_extend);
            if (best - v > x_drop_gapped) break;
            Hprev[j] = v;
            hi = j;
        }
        for (int i = 1; i <= qn; ++i) {
            fill(Hcur.begin(), Hcur.end(), NEG);
            fill(Ecur.begin(), Ecur.end(), NEG);
            Score f = NEG; // gap en db
            int new_lo = -1, new_hi = -1;
            const int jend = min(hi + 1, dn);
            for (int j = lo; j <= jend; ++j) {
                Score m = NEG; // diagonal: match/mismatch
                if (j >= 1 && Hprev[j - 1] > NEG / 2) m = Hprev[j - 1] + blosum[q[i - 1]][d[j - 1]];
                Score e = NEG; // gap en query
                if (j >= 1) {
                    Score fh = (Hcur[j - 1] > NEG / 2) ? Hcur[j - 1] - gap_open - gap_extend : NEG;
                    Score fe = (Ecur[j - 1] > NEG / 2) ? Ecur[j - 1] - gap_extend : NEG;
                    e = max(fh, fe);
                }
                Score fh = (Hprev[j] > NEG / 2) ? Hprev[j] - gap_open - gap_extend : NEG;
                Score ff = (f > NEG / 2) ? f - gap_extend : NEG;
                f = max(fh, ff);

                Score h = max(m, max(e, f));
                if (h <= NEG / 2 || best - h > x_drop_gapped) continue;
                Hcur[j] = h;
                Ecur[j] = e;
                if (h > best) best = h, q_reach = (Length) i, d_reach = (Length) j;
                if (new_lo < 0) new_lo = j;
                new_hi = j;
            }
            if (new_lo < 0) break;
            lo = new_lo;
            hi = new_hi;
            swap(Hprev, Hcur);
        }
        return best;
    }

    Extension gapped(const uint8_t *db, Length db_len, Position qp, Position dp) {
        Score anchor = blosum[query[qp]][db[dp]];
        Length rq, rd, lq, ld; // cuanto llego cada mitad desde el ancla
        Score right = half_gapped(query + qp + 1, query_len - qp - 1, db + dp + 1, db_len - dp - 1, rq, rd);
        vector<uint8_t> qr(query, query + qp), dr(db, db + dp);
        reverse(qr.begin(), qr.end());
        reverse(dr.begin(), dr.end());
        Score left = half_gapped(qr.data(), (Length) qr.size(), dr.data(), (Length) dr.size(), lq, ld);
        return {(Position) (dp - ld),     (Position) (dp + rd + 1), (Position) (qp - lq),
                (Position) (qp + rq + 1), anchor + left + right,    0};
    }
};

array<double, 25> load_background_freq(const string &filename, LetterCodeMap &map) {
    array<double, 25> res{};
    fstream f(filename.c_str(), ios::in);
    char c;
    double p;
    double total = 0;
    while (f >> c) {
        f >> p;
        res[map[c].value] = p, total += p;
    }
    return res;
}

struct GappedConfig {
    double lambda, K, H_bits;
};

GappedConfig load_karlin_gapped(const string &filename, int gap_open, int gap_extend) {
    fstream f(filename.c_str(), ios::in);
    GappedConfig res{0, 0, 0};
    int go, ge;
    double lambda, k, h;
    while (f >> go) {
        f >> ge >> lambda >> k >> h;
        if (go == gap_open && ge == gap_extend) res = {lambda, k, h};
    }
    return res;
}

struct KarlinAltschul {
    array<double, 25> P{};
    double long lambda, K, H;

    KarlinAltschul(BLOSUM &blosum, bool with_gaps, array<double, 25> &fallback, GappedConfig &gapped,
                   const uint8_t *residues, uint64_t nres) {
        if (with_gaps) {
            // Si hay gaps -> no importa frequencia
            lambda = gapped.lambda;
            K = gapped.K;
            H = gapped.H_bits * log(2.0); // la tabla esta en bits, aca van nats
            return;
        }

        vector<uint8_t> standard;
        array<bool, 25> is_standard{};
        forn(i, 25)
            if (fallback[i] > 0) standard.pb((uint8_t) i), is_standard[i] = true;

        double total = 0;
        if (residues)
            for (uint64_t i = 0; i < nres; i++)
                if (is_standard[residues[i]]) P[residues[i]] += 1, total += 1;

        if (total > 0)
            for (uint8_t i: standard) P[i] /= total;
        else
            for (uint8_t i: standard) P[i] = fallback[i];
        // biseccion sobre F(l) = sum p_i p_j e^(l s_ij) - 1
        // F(0)=0, F'(0)<0 (score esperado al azar es negativo) -> 1 sola raiz positiva
        auto F = [&](double l) {
            double s = 0;
            for (uint8_t i: standard)
                for (uint8_t j: standard) s += P[i] * P[j] * exp(l * blosum[i][j]);
            return s - 1.0;
        };
        double lo = 1e-9, hi = 2.0;
        forn(it, 200) {
            double mid = (lo + hi) / 2;
            (F(mid) < 0 ? lo : hi) = mid;
        }
        lambda = (lo + hi) / 2;
        H = 0;
        for (uint8_t i: standard)
            for (uint8_t j: standard) H += P[i] * P[j] * exp(lambda * blosum[i][j]) * blosum[i][j];
        H *= lambda;
        K = gapped.K;
    }

    double bits(Score score) const { return (lambda * score - log(K)) / log(2.0); }

    double evalue(Score score, uint32_t m, double n, uint32_t num_seqs) const {
        // Mal borde -> alineamiento en ultimo residuio -> espacio menor que m*n -> iterar (en 5 converge) [Copiado de NCBI]
        double L = 0, me = m, ne = n;
        forn(it, 5) {
            me = max((double) (1.0 / K), (double) m - L);
            ne = max((double) (1.0 / K), n - (double) num_seqs * L);
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
    GappedConfig gapped_config = load_karlin_gapped("assets/KARLIN_GAPPED", GAP_OPEN, GAP_EXTEND);

    DB database(index_map);
    const string source_file = "assets/SEQUENCE_large";
    const string index_file = "cache/SEQUENCE_large.zbx";
    if (!database.open(source_file, index_file)) {
        cerr << "no pude abrir el indice " << index_file << el;
        return 1;
    }

    WordPermutation word_permutation(blosum, WORD_THRESHOLD);
    //word_permutation.auto_get_permutations("permutation.cache");
    string query_string;
    cin >> query_string;

    Length qLen = sz(query_string);
    vector<uint8_t> query;
    forn(i, qLen) query.pb((uint8_t) index_map[query_string[i]].value);

    HitBuffer buffer;

    // Pre-reserve buffer size to total sum of postings -> vecindarios memoized pq no se recomputan
    size_t total = 0;
    for_each_word(query.data(), qLen, [&](Position, WordValue word) {
        for (WordValue neighbour: word_permutation[word]) total += database[neighbour].size();
    });
    buffer.reserve(total);

    for_each_word(query.data(), qLen, [&](Position i, WordValue word) {
        for (WordValue neighbour: word_permutation[word]) {
            DB::Appearence matches = database[neighbour]; // dos punteros al mmap
            forn(k, matches.n) buffer.add(Hit(matches.access[k], i, matches.pos[k], qLen));
        }
    });
    cout << "hits " << buffer.count << el;

    HitBuffer::Radix radix;
    radix.sort(buffer);

    SeedBuffer seeds;
    two_hit(buffer, seeds, TWO_HIT_WINDOW, WORDSIZE); // min_gap = wordsize
    cout << "seeds " << seeds.count << el;

    Extender extender(blosum, query.data(), qLen, X_DROP_UNGAPPED, X_DROP_GAPPED, GAP_OPEN, GAP_EXTEND);
    KarlinAltschul ku(blosum, false, background, ungapped_config, database.residues, database.nres);
    KarlinAltschul kg(blosum, true, background, gapped_config, database.residues, database.nres);

    vector<Extension> hsps;
    DiagonalKey covered_key = 0;
    Position covered_until = 0; // arranca en 0: db_pos < 0 nunca es cierto sin signo
    forn(s, seeds.count) {
        Hit &hit = buffer[seeds[s].second];
        if (hit.key == covered_key && hit.db_pos < covered_until) continue;
        AccessNumber acc = hit.access();
        Extension e = extender.ungapped(database.sequence_of(acc), database.length_of(acc), hit.query_pos, hit.db_pos);
        e.access = acc;
        covered_key = hit.key, covered_until = e.db_end;
        if (ku.bits(e.score) >= HSP_MIN_BITS) hsps.pb(e);
    }
    cout << "hsps " << sz(hsps) << el;

    for (auto &h: hsps) {
        Extension g = extender.gapped(database.sequence_of(h.access), database.length_of(h.access),
                                      (h.q_start + h.q_end) / 2, (h.db_start + h.db_end) / 2);
        g.access = h.access;
        if (g.score > h.score) h = g;
    }

    map<AccessNumber, Extension> best;
    for (auto &h: hsps)
        if (!best.count(h.access) || h.score > best[h.access].score) best[h.access] = h;
    vector<Extension> results;
    for (auto &kv: best) results.pb(kv.second);
    sort(results.begin(), results.end(), [](Extension &a, Extension &b) { return a.score > b.score; });

    double db_len = (double) database.nres;
    uint32_t num_seqs = database.sequence_count();
    printf("\n%-8s %7s %7s %11s   %-13s %s\n", "access", "score", "bits", "E-value", "query", "base");
    int shown = 0;
    for (auto &h: results) {
        double e = kg.evalue(h.score, qLen, db_len, num_seqs);
        if (e > MAX_EVALUE) break;
        if (++shown > MAX_RESULTS) break;
        printf("%-8u %7d %7.1f %11.2g   %5u-%-7u %6u-%-8u\n", h.access, h.score, (double) kg.bits(h.score), e,
               h.q_start, h.q_end, h.db_start, h.db_end);
    }
    if (!shown)
        forn(i, min(5, sz(results))) {
            Extension &h = results[i];
            printf("%-8u %7d %7.1f %11.2g   %5u-%-7u %6u-%-8u  NS\n", h.access, h.score, (double) kg.bits(h.score),
                   kg.evalue(h.score, qLen, db_len, num_seqs), h.q_start, h.q_end, h.db_start, h.db_end);
        }
    cout << el << "lambda " << (double) ku.lambda << "  K " << (double) ku.K << "  H " << (double) ku.H << " (from DB)"
         << el;
    return 0;
} catch (exception &e) {
    cerr << "error: " << e.what() << el;
    return 1;
}
