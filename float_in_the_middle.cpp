#include <algorithm>
#include <cstdint>
#include <vector>
#include <chrono>
#include <cmath>

// Java Random

constexpr uint64_t MASK = (1ULL << 48) - 1;
constexpr uint64_t A = 0x5deece66d;
constexpr uint32_t B = 11;

inline void nextSeed(uint64_t* rand) {
    *rand = (*rand * A + B) & MASK;
}

inline uint32_t nextFloatBits(uint64_t* rand) {
    nextSeed(rand);
    return *rand >> 24;
}

constexpr uint32_t FLOAT_TO_INT_UNIT = 1u << 24;
struct FloatRange {
    uint32_t min;
    uint32_t max;

    FloatRange(float min, float max) {
        this->min = std::floor(min * FLOAT_TO_INT_UNIT);
        this->max = std::ceil(max * FLOAT_TO_INT_UNIT);
    }
}

// -----------------------------------------------------------------
// The thing

constexpr uint32_t NUM_BUCKETS = 10;
constexpr uint32_t BUCKET_RANGE_SIZE = std::floor(0.2f * NUM_BUCKETS) + 1;
// Buckets are 1/NUM_BUCKETS wide. This makes each query for a range [n, n+0.2]
// check BUCKET_RANGE_SIZE buckets per layer.

constexpr uint32_t NUM_LAYERS = 6;
// The number of layers is a key parameter here, because it dictates
// not only the size of the LUT, but also the overlap in bucket
// signatures for the lower bits. Pushing this higher increases filtering power 
// and the size of the LUT, and reduces signature overlap.

constexpr uint32_t UNIQUE_SIGNATURES = std::pow(NUM_BUCKETS, NUM_LAYERS);

// LUT structure used by both upper bit list and lower bit list

struct LUTEntry {
    uint32_t bits;
    uint32_t values[NUM_LAYERS];
};
struct LUT {
    // vector index = bucket signature
    std::vector<LUTEntry>* enTree = nullptr;

    LUT() {
        enTree = new std::vector<LUTEntry>[UNIQUE_SIGNATURES];
    }

    ~LUT() {
        if (enTree != nullptr) {
            delete[] enTree;
        }
    }

    void addLowerBits(uint32_t lowerBits, std::vector<) {
        LUTEntry newEntry;
        newEntry.bits = lowerBits;

        // store nextFloat contributions
        uint64_t state = static_cast<uint64_t>(lowerBits);
        uint32_t bucketRangeSig = 0;

        for (int i = 0; i < NUM_LAYERS; i++) {
            state = (state * A + B) & MASK;
            uint32_t bits = state >> 24;
            newEntry.values[i] = bits >> 24;

            bucketRangeSig *= NUM_BUCKETS;
            bucketRangeSig += bits / NUM_BUCKETS;
        }

        enTree[bucket].push_back(newEntry);
    }

    void addUpperBits(uint32_t upperBits) {
        LUTEntry newEntry;
        newEntry.bits = lowerBits;

        // store nextFloat contributions
        uint64_t state = static_cast<uint64_t>(upperBits) << 24;
        uint32_t bucketSig = 0;

        for (int i = 0; i < NUM_LAYERS; i++) {
            state = (state * A) & MASK;
            uint32_t bits = state >> 24;
            newEntry.values[i] = bits >> 24;

            bucketSig *= NUM_BUCKETS;
            bucketSig += bits / NUM_BUCKETS;
        }

        enTree[bucketSig].push_back(newEntry);
    }

    std::vector<LUTEntry>& getEntriesForBucket(uint32_t bucket) {
        return enTree[bucket];
    }
};



void float_in_the_middle(LUT& lowerLut, LUT& upperLut) {
    // TODO

    // iterate over bucket range signatures in the lower LUT.
    // for each non-empty vector, iterate over all buckets in the upper LUT
    // that fall into the bucket range.

    // In the innermost loop, pair each lower vec element 
    // with each upper vec element and check constraint satisfaction.
}

// -----------------------------------------------------------------
// Time to actually use the thing now

std::chrono::timepoint_t timer_start() {
    return std::chrono::steady_clock::now();
}

void timer_stop(std::chrono::timepoint_t start) {
    auto end = std::chrono::steady_clock::now();
    double ms = (end-start).count() * 1e-6;
    std::printf("Took %f ms.\n", ms);
}

int main() {
    std::vector<FloatRange> constraints({{
        FloatRange(0.0f, 0.2f),
        FloatRange(0.0f, 0.2f),
        FloatRange(0.0f, 0.2f),
        FloatRange(0.0f, 0.2f),
        FloatRange(0.0f, 0.2f),
        FloatRange(0.0f, 0.2f)
    }});
    
    std::printf("Initializing LUTs...\n");
    auto t0 = timer_start();
    LUT upperLut();
    LUT lowerLut();
    timer_stop(t0);

    std::printf("Filling LUTs...\n");
    auto t1 = timer_start();
    for (uint32_t bits = 0; bits < (1u << 24); bits++) {
        lowerLut.addLowerBits(bits, constraints);
        upperLut.addUpperBits(bits);
    }
    timer_stop(t1);

    std::printf("Float-in-the-middle search running...\n");
    
    auto t2 = timer_start();
    float_in_the_middle(lowerLut, upperLut);
    timer_stop(t2);

    std::printf("Finished.\n");
}