package levy

import (
	"math"
	"math/rand"
)

// Generates sample from general Levy-stable distibution.
func Sample(alpha float64, beta float64, mu float64, sigma float64) float64 {
	if beta == 0.0 {
		return symmetric(alpha, mu, sigma)
	}

	v := math.Pi * (rand.Float64() - 0.5)

	w := 0.0
	for w == 0.0 {
		w = rand.ExpFloat64()
	}

	if alpha == 1.0 {
		x := ((0.5*math.Pi+beta*v)*math.Tan(v) -
			beta*math.Log(0.5*math.Pi*w*math.Cos(v)/(0.5*math.Pi+beta*v))) / (0.5 * math.Pi)
		return sigma*x + beta*sigma*math.Log(sigma)/(0.5*math.Pi) + mu
	}

	t := beta * math.Tan(0.5*math.Pi*alpha)
	s := math.Pow(1.0+t*t, 1.0/(2.0*alpha))
	b := math.Atan(t) / alpha
	x := s * math.Sin(alpha*(v+b)) *
		math.Pow(math.Cos(v-alpha*(v+b))/w, (1.0-alpha)/alpha) /
		math.Pow(math.Cos(v), 1.0/alpha)
	return sigma*x + mu
}

// Generates a sample from a Guassian distribution (alpha=2, beta=0).
func GaussSample(mu float64, sigma float64) float64 {
	return symmetric(2.0, mu, sigma/math.Sqrt(2.0))
}

// Generates a sample from a Cauchy distribution (alpha=1, beta=0).
func CauchySample(mu float64, sigma float64) float64 {
	return symmetric(1.0, mu, sigma)
}

// Generates a sample from a Levy distribution (alpha=0.5, beta=1).
func LevySample(mu float64, sigma float64) float64 {
	return Sample(0.5, 1.0, mu, sigma)
}

func symmetric(alpha float64, mu float64, sigma float64) float64 {

	u := math.Pi * (rand.Float64() - 0.5)

	if alpha == 1.0 {
		return sigma*math.Tan(u) + mu
	}

	v := 0.0
	for v == 0.0 {
		v = rand.ExpFloat64()
	}

	if alpha == 2.0 {
		return 2.0*math.Sin(u)*math.Sqrt(v)*sigma + mu
	}

	t := math.Sin(alpha*u) / math.Pow(math.Cos(u), 1.0/alpha)
	s := math.Pow(math.Cos((1.0-alpha)*u)/v, (1.0-alpha)/alpha)
	return sigma*t*s + mu
}
