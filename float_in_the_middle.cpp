#include <cstdint>
#include <cstdio>
#include <vector>
#include <chrono>
#include <cmath>

// Java Random

constexpr uint64_t MASK_48 = (1ULL << 48) - 1;
constexpr uint32_t MASK_24 = (1U << 24) - 1;
constexpr uint64_t A = 0x5deece66d;
constexpr uint32_t B = 11;

inline void nextSeed(uint64_t* rand) {
    *rand = (*rand * A + B) & MASK_48;
}

inline uint32_t nextFloatBits(uint64_t* rand) {
    nextSeed(rand);
    return *rand >> 24;
}

constexpr uint32_t FLOAT_TO_INT_UNIT = 1u << 24;
struct FloatRange {
    uint32_t min;
    uint32_t max;
    uint64_t lcgA;
    uint64_t lcgB;

    FloatRange(float min, float max, int64_t lcgSteps) {
        this->min = std::floor(min * FLOAT_TO_INT_UNIT);
        this->max = std::ceil(max * FLOAT_TO_INT_UNIT);
        this->lcgA = A;
        this->lcgB = B;
        lcgCombine(lcgSteps);
    }

private:
    void lcgCombine(int64_t steps) {
        uint64_t real_steps = static_cast<uint64_t>(steps) & MASK_48;
        uint64_t combA = lcgA;
        uint64_t combB = lcgB;
        lcgA = 1;
        lcgB = 0;

        // A2(A1*x + B1) + B2 = A1A2 * x + A2B1 + B2
        for (int i = 0; i < 48; i++) {
            if (real_steps & 1) {
                lcgB = (lcgB * combA + combB) & MASK_48;
                lcgA = (lcgA * combA) & MASK_48;
            }
            combB = (combB * combA + combB) & MASK_48;
            combA = (combA * combA) & MASK_48;
            real_steps >>= 1;
        }
    }
};

// -----------------------------------------------------------------
// The thing

constexpr uint32_t LOWER_BIT_COUNT = 24;
// Idk what will be best here yet, will find out empirically.

constexpr uint32_t NUM_BUCKETS = 10; // 0.0-0.1, 0.1-0.2, 0.2-0.3
constexpr uint32_t BUCKET_WIDTH = (FLOAT_TO_INT_UNIT + NUM_BUCKETS - 1) / NUM_BUCKETS;
constexpr uint32_t BUCKET_RANGE_SIZE = 3;
// Buckets are 1/NUM_BUCKETS wide. This makes each query for a range [n, n+0.2]
// check BUCKET_RANGE_SIZE = floor(0.2 * NUM_BUCKETS) + 1 buckets per layer.

constexpr uint32_t NUM_LAYERS = 6;
// The number of layers is a key parameter here, because it dictates
// not only the size of the LUT, but also the overlap in bucket
// signatures for the lower bits. Pushing this higher increases filtering power 
// and the size of the LUT, and reduces signature overlap.

constexpr uint32_t pow(uint32_t a, uint32_t n) {
    uint32_t r = 1;
    for (uint32_t i = 0; i < n; i++) {
        r *= a;
    }
    return r;
}
constexpr uint32_t UNIQUE_SIGNATURES = pow(NUM_BUCKETS, NUM_LAYERS);

// LUT structure used by both upper bit list and lower bit list

struct LUTEntry {
    uint32_t bits;
    uint32_t values[NUM_LAYERS];
};
struct LUT {
    // vector index = bucket signature
    std::vector<LUTEntry>* enTree = nullptr;
    const std::vector<FloatRange>& constraints;

    LUT(const std::vector<FloatRange>& constraints) : constraints(constraints) {
        enTree = new std::vector<LUTEntry>[UNIQUE_SIGNATURES];
    }

    ~LUT() {
        if (enTree != nullptr) {
            delete[] enTree;
        }
    }

    void addLowerBits(uint32_t lowerBits) {
        LUTEntry newEntry;
        newEntry.bits = lowerBits;

        // store nextFloat contributions
        uint32_t bucketRangeSig = 0;

        for (int i = 0; i < NUM_LAYERS; i++) {
            uint64_t state = (lowerBits * constraints[i].lcgA + constraints[i].lcgB) & MASK_48;
            uint32_t bits = state >> 24;
            newEntry.values[i] = bits;

            bucketRangeSig *= NUM_BUCKETS;
            uint32_t rangeMin = FLOAT_TO_INT_UNIT + constraints[i].min - bits;
            if (rangeMin > FLOAT_TO_INT_UNIT) rangeMin -= FLOAT_TO_INT_UNIT;

            bucketRangeSig += rangeMin / BUCKET_WIDTH; // FIXME this could be inaccurate
        }

        if (bucketRangeSig > UNIQUE_SIGNATURES) {
            printf("BAD LOWER INDEX: %u\n", bucketRangeSig);
            return;
        }
        enTree[bucketRangeSig].push_back(newEntry);
    }

    void addUpperBits(uint32_t upperBits) {
        LUTEntry newEntry;
        newEntry.bits = upperBits;

        // store nextFloat contributions
        uint32_t bucketSig = 0;

        for (int i = 0; i < NUM_LAYERS; i++) {
            uint32_t bits = static_cast<uint32_t>((upperBits * constraints[i].lcgA) & MASK_24);
            newEntry.values[i] = bits;

            bucketSig *= NUM_BUCKETS;
            bucketSig += bits / BUCKET_WIDTH; // FIXME this could be inaccurate
        }

        if (bucketSig > UNIQUE_SIGNATURES) {
            printf("BAD UPPER INDEX: %u\n", bucketSig);
            return;
        }
        enTree[bucketSig].push_back(newEntry);
    }

    std::vector<LUTEntry>& getEntriesForBucket(uint32_t bucket) {
        if (bucket >= UNIQUE_SIGNATURES) {
            std::printf("Out of bounds access - %u\n", bucket);
            return enTree[0];
        }
        return enTree[bucket];
    }
};

uint64_t counter = 0;

void lookup_bucket_range(std::vector<LUTEntry>& lowerEntries, LUT& upperLut, uint32_t remaining_sig, uint32_t partial_bucket, uint32_t depth) {
    if (depth == 6) {
        for (int i = 0; i < BUCKET_RANGE_SIZE; i++) {
            uint32_t new_partial = partial_bucket * NUM_BUCKETS;
            uint32_t digit = (remaining_sig + i) % NUM_BUCKETS;
            std::vector<LUTEntry>& upperEntries = upperLut.getEntriesForBucket(new_partial);
            
            for (auto& low : lowerEntries) {
                for (auto& high : upperEntries) {
                    // Full check against stored FloatRange constraints
                    bool passed = true;
                    const auto& constr = upperLut.constraints;
                    for (int c = 0; c < NUM_LAYERS; c++) {
                        uint32_t combinedValue = (low.values[c] + high.values[c]) & MASK_24;
                        passed &= constr[i].min <= combinedValue && combinedValue <= constr[i].max;
                        if (!passed) break; 
                    }

                    if (passed) {
                        counter++;
                    }
                }
            }
        }
        return;
    }
    
    for (int i = 0; i < BUCKET_RANGE_SIZE; i++) {
        uint32_t new_partial = partial_bucket * NUM_BUCKETS;
        uint32_t digit = (remaining_sig + i) % NUM_BUCKETS;
        //std::printf("partial: %u  digit: %u\n", new_partial, digit);
        lookup_bucket_range(lowerEntries, upperLut, remaining_sig/NUM_BUCKETS, new_partial+digit, depth+1);
    }
}

void float_in_the_middle(LUT& lowerLut, LUT& upperLut) {
    for (uint32_t startBucketSig = 0; startBucketSig < UNIQUE_SIGNATURES; startBucketSig++) {
        //std::printf("bucket: %u\n", startBucketSig);
        std::vector<LUTEntry>& lowerEntries = lowerLut.getEntriesForBucket(startBucketSig);
        if (lowerEntries.empty())
            continue;

        lookup_bucket_range(lowerEntries, upperLut, startBucketSig, 0, 1);
    }
}

// -----------------------------------------------------------------
// Time to actually use the thing now

typedef std::chrono::time_point<std::chrono::steady_clock> timepoint;

timepoint timer_now() {
    return std::chrono::steady_clock::now();
}

void timer_get_elapsed(timepoint start) {
    auto end = timer_now();
    double ms = (end-start).count() * 1e-6;
    std::printf("Took %f ms.\n", ms);
}

int main() {
    timepoint t3;
    {
        std::vector<FloatRange> constraints({{
            FloatRange(0.0f, 0.2f, 1),
            FloatRange(0.0f, 0.2f, 2),
            FloatRange(0.0f, 0.2f, 3),
            FloatRange(0.0f, 0.2f, 4),
            FloatRange(0.0f, 0.2f, 5),
            FloatRange(0.0f, 0.2f, 6)
        }});
        
        std::printf("Initializing LUTs...\n");
        auto t0 = timer_now();
        LUT upperLut(constraints);
        LUT lowerLut(constraints);
        timer_get_elapsed(t0);

        std::printf("Filling LUTs...\n");
        auto t1 = timer_now();
        for (uint32_t bits = 0; bits < (1u << LOWER_BIT_COUNT); bits++) {
            lowerLut.addLowerBits(bits);
        }
        for (uint32_t bits = 0; bits < (1u << 48-LOWER_BIT_COUNT); bits++) {
            upperLut.addUpperBits(bits);
        }
        timer_get_elapsed(t1);

        std::printf("Float-in-the-middle search running...\n");
        
        auto t2 = timer_now();
        float_in_the_middle(lowerLut, upperLut);
        timer_get_elapsed(t2);

        std::printf("%llu states passed full check,\n", counter);
        std::printf("%llu states were expected.\n", (1ULL << 48) / pow(5, NUM_LAYERS));

        t3 = timer_now();
        std::printf("Deallocating memory...\n");
    }
    timer_get_elapsed(t3);
    std::printf("Done.\n");
}
