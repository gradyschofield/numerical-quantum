
class BigSum {
    float positiveAccumulators[25];
    float negativeAccumulators[25];
    float orderLimits[25];
    int binOffset;
    int orders;

public:

    __device__ BigSum() {
        float minBin = 1E-5;
        orders = 25;
        binOffset = round(log10(minBin));
        double lowestLimit = pow(10, binOffset);
        for (int i = 0; i < orders; ++i) {
            orderLimits[i] = lowestLimit * pow(10, i);
            positiveAccumulators[i] = 0;
            negativeAccumulators[i] = 0;
        }
    }
    __device__ void add(float x) {
        int bin = min(orders-1, (int)max(0.0, round(log10(fabs(x))) - binOffset));
        float positiveInc = max(0.0f, x);
        float negativeInc = min(0.0f, x);
        positiveAccumulators[bin] += positiveInc;
        negativeAccumulators[bin] += negativeInc;
    }

    __device__ operator double() const {
        double positiveSum = 0;
        double negativeSum = 0;
        for(int i = 0; i < orders; ++i) {
            positiveSum += positiveAccumulators[i];
        }
        for(int i = 0; i < orders; ++i) {
            negativeSum += negativeAccumulators[i];
        }
        return positiveSum + negativeSum;
    }

    __device__ void operator+=(double x) {
        add(x);
    }

    __device__ void operator-=(double x) {
        add(-x);
    }

    __device__ void carryAccumulators() {
        for (int bin = 0; bin < orders - 1; ++bin) {
            float setValue = positiveAccumulators[bin] > orderLimits[bin] ? 0 : positiveAccumulators[bin];
            float carryValue = positiveAccumulators[bin] > orderLimits[bin] ? positiveAccumulators[bin] : 0;
            positiveAccumulators[bin + 1] += carryValue;
            positiveAccumulators[bin] = setValue;
        }
        for (int bin = 0; bin < orders - 1; ++bin) {
            float setValue = negativeAccumulators[bin] < -orderLimits[bin] ? 0 : negativeAccumulators[bin];
            float carryValue = negativeAccumulators[bin] < -orderLimits[bin] ? negativeAccumulators[bin] : 0;
            negativeAccumulators[bin + 1] += carryValue;
            negativeAccumulators[bin] = setValue;
        }
    }

};

