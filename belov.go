package levy

import (
	"Pb/quad"
	"math"
)

// PdfBelov is a generator for a Levy-stable distribution that uses Nolan's
// representation and integration scheme of I. A. Belov
// This splits the integral into two pieces and uses fixed abiscca
// methods to do the integral. The method does not work well when the
// integrand gets oscilatory - ie for small values of alpha ~ 0.5.
// Probably better to stick to PdfZ methods.
type PdfBelov struct{}

func NewPdfBelov() *PdfBelov {
	return &PdfBelov{}
}

func (pdf *PdfBelov) ScaledValue(x, alpha, beta float64) (float64, error) {

	if alpha == 2 {
		return math.Exp(-0.25*x*x) / math.Sqrt(4.0*math.Pi), nil
	} else if alpha == 1 && beta == 0 {
		return 1.0 / ((1.0 + x*x) * math.Pi), nil
	} else {
		tval := 8.0
		v1 := quad.NewGaussLaguerre().Integrate(
			func(t float64) float64 {
				return componentBelov(x, t+tval, alpha, beta)
			}, 0.0)
		v2 := quad.NewGaussLegendre1024().Integrate(
			func(t float64) float64 {
				return componentBelov(x, t, alpha, beta)
			}, 0.0, tval)
		return v1 + v2, nil
	}
	return 0, nil
}

func componentBelov(x float64, t float64, alpha float64, beta float64) float64 {
	h := x*t + beta*(t-math.Pow(t, alpha))*math.Tan(0.5*math.Pi*alpha)
	return math.Cos(h) * math.Exp(-math.Pow(t, alpha)) / math.Pi
}
