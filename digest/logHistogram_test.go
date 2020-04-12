package tdigest

import (
	"math"
	"testing"
)

func TestBoundaries(t *testing.T) {
	hist, _ := NewLogHistogram(0.1, 1000, 0.1)

	bounds := hist.GetBounds()
	oldV := 0.0
	for i, v := range bounds {
		if oldV == 0 {
			oldV = v
		}
		if v/oldV > math.Exp(0.1) {
			t.Errorf("large jump of %.4f at index %d", v/oldV, i)
		}
		oldV = v
	}
}

func TestInverse(t *testing.T) {
	for x := 0.001; x <= 100; x += 1e-3 {
		log := approxLog2(x)
		system := math.Log2(x)
		if math.Abs(log-system) > 1e-2 {
			t.Errorf("pow2(approxLog, %.3f, %.3f", x, log)
		}

		roundTrip := pow2(log)
		if math.Abs(x-roundTrip) > 1e-13 {
			t.Errorf("pow2(approxLog, %.3f, %.3f, %.3f", x, log, roundTrip)
		}
	}
}

func TestLogHistogram_Add(t *testing.T) {
	hist, _ := NewLogHistogram(0.1, 1000, 0.1)
	hist2, _ := NewLogHistogram(0.1, 1000, 0.1)
	hist3, _ := NewLogHistogram(0.1, 1000, 0.1)
	hist4, _ := NewLogHistogram(0.1, 1000, 0.1)

	for x := 0.05; x < 2000; x *= math.Sqrt(1.1) {
		hist.Add(x)
		if x < 20 {
			hist2.Add(x)
		} else if x < 300 {
			hist3.Add(x)
		} else {
			hist4.Add(x)
		}
	}
	counts := make([]int, 20)
	kx := hist.GetCounts()
	for i, k := range kx {
		if i > 0 && i < len(kx)-1 {
			counts[k]++
		}
	}
	if kx[0] != 17 {
		t.Error("Should have had large underflow count")
	}
	if kx[len(kx)-1] != 16 {
		t.Error("Should have had large overflow count")
	}
	if counts[1]+counts[3] > 6 || counts[2] < int(0.9*float64(len(kx))) || counts[1]+counts[2]+counts[3] != len(kx)-2 {
		t.Error("counts should be 2 ± 1 (and almost always 2")
	}

	err := hist2.AddHistograms(hist3, hist4)
	if err != nil {
		t.Error("AddHistograms threw an error")
	}
	counts2 := hist2.GetCounts()
	errors := 0
	for i, count := range hist.GetCounts() {
		if count != counts2[i] {
			errors++
			t.Errorf("bad merge count[%d] should be %d, was %d", i, count, counts2[i])
			if errors > 4 {
				break
			}
		}
	}
}
