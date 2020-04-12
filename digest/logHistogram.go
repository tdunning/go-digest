package tdigest

import (
	"errors"
	"fmt"
	"math"
)

// LogHistogram has bins that are very nearly logarithmically spaced from a min to a max value.
type LogHistogram struct {
	min, max             float64
	logFactor, logOffset float64
	count                []uint32
}

func NewLogHistogram(min, max, epsilonFactor float64) (*LogHistogram, error) {
	if epsilonFactor == 0 {
		epsilonFactor = 0.1
	}
	// slight guard band to account for inaccuracy in approxLog2
	epsilonFactor *= 0.98

	if max <= 2*min {
		return nil, errors.New(fmt.Sprintf("illegal/nonsensical min, max (%.2f, %.2f)", min, max))
	}
	if min <= 0 || max <= 0 {
		return nil, errors.New("min and max must be positive")
	}
	if epsilonFactor < 1e-6 || epsilonFactor > 0.5 {
		return nil, errors.New(
			fmt.Sprintf("Unreasonable number of bins per decade %.2g. Expected value in range [1e-6,0.5]",
				epsilonFactor))
	}

	tmp := math.Log(2) / math.Log(1+epsilonFactor)
	r := &LogHistogram{
		min:       min,
		max:       max,
		logFactor: tmp,
		logOffset: approxLog2(min) * tmp,
	}
	binCount := r.BucketIndex(max) + 1
	if binCount > 10000 {
		return nil, errors.New(
			fmt.Sprintf("Excessive number of bins %d resulting from min,max = %.2g, %.2g",
				binCount, min, max))
	}
	r.count = make([]uint32, binCount)
	return r, nil
}

/**
 * Approximates log_2(value) by abusing floating point hardware. The floating point exponent
 * is used to get the integer part of the log. The mantissa is then adjusted with a second order
 * polynomial to get a better approximation. The error is bounded to be less than ±0.01 and is
 * zero at every power of two (which also implies the approximation is continuous).
 *
 * @param value The argument of the log
 * @return log_2(value) (within an error of about ± 0.01)
 */
func approxLog2(value float64) float64 {
	valueBits := math.Float64bits(value)
	exponent := int64((valueBits&0x7ff0000000000000)>>52) - 1024
	m := math.Float64frombits((valueBits & 0x800fffffffffffff) + 0x3ff0000000000000)
	return m*(2-(1.0/3)*m) + float64(exponent) - (2.0 / 3.0)
}

/**
 * Computes an approximate value of 2^x. This is done as an exact inverse of #approxLog2 so
 * that bin boundaries can be computed exactly.
 *
 * @param x The power of 2 desired.
 * @return 2^x approximately.
 */
func pow2(x float64) float64 {
	exponent := math.Floor(x) - 1
	x = x - exponent
	m := 3 - math.Sqrt(7-3*x)
	return math.Pow(2, exponent+1) * m
}

func (hist LogHistogram) BucketIndex(x float64) int {
	return int(approxLog2(x)*hist.logFactor - hist.logOffset)
}

func (hist LogHistogram) LowerBound(k int) float64 {
	return pow2((float64(k) + hist.logOffset) / hist.logFactor)
}

func (hist LogHistogram) Add(v float64) {
	hist.count[hist.bucket(v)]++
}

func (hist LogHistogram) GetBounds() []float64 {
	r := make([]float64, len(hist.count))
	for i := range r {
		r[i] = hist.LowerBound(i)
	}
	return r
}

func (hist LogHistogram) GetCounts() []uint32 {
	return hist.count
}

func (hist *LogHistogram) AddHistograms(others ...*LogHistogram) error {
	for _, other := range others {
		if other.min != hist.min || other.max != hist.max || len(other.count) != len(hist.count) {
			return errors.New("can only merge histograms with identical bounds and precision")
		}
		for i, k := range other.count {
			hist.count[i] += k
		}
	}
	return nil
}

// exposed for testing
func (hist LogHistogram) bucket(x float64) int {
	if x <= hist.min {
		return 0
	} else if x >= hist.max {
		return len(hist.count) - 1
	} else {
		return hist.BucketIndex(x)
	}
}
