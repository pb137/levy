package levy

import (
	"fmt"
	"math"
)

// PdfBergstrom is a generator for the tail of a Levy-stable distribution that uses a
// Bergstrom representation. The representation is only suitable for large values of x.
// Convergence improves with increasing values of alpha.
type PdfBergstrom struct {
	eps   float64
	limit int
}

// Creates a new generator for tail distribution.
func NewPdfBergstrom() *PdfBergstrom {
	return &PdfBergstrom{1e-12, 30}
}

func (pdf *PdfBergstrom) ScaledValue(x, alpha, beta float64) (float64, error) {

	var err error
	zeta := beta * math.Tan(0.5*math.Pi*alpha)
	eps := pdf.eps / x * math.Pi
	done := false
	n := 1
	sum := 0.0
	for !done {
		a := 1.0
		if n%2 == 0 {
			a = -1.0
		}
		a *= math.Gamma(float64(n)*alpha+1) / math.Gamma(float64(n)+1)
		a *= math.Pow(1+zeta*zeta, 0.5*float64(n))
		a *= math.Sin(float64(n) * (0.5*math.Pi*alpha + math.Atan(zeta)))

		delta := a * math.Pow(x, -alpha*float64(n))
		sum += delta

		if math.Abs(delta) < eps {
			done = true
		}
		if n >= pdf.limit {
			done = true
			err = fmt.Errorf("Iteration limit in tail approximation exceeded (%d)", pdf.limit)
		}
		n++
	}
	sum /= x * math.Pi
	return sum, err
}
