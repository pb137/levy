package levy

import (
	"fmt"
	"github.com/pb137/quad"
	"github.com/pb137/solve"
	"math"
)

// PdfZ is a generator for a Levy-stable distribution that uses Zolatarev
// representation and a custom integrator. In this representation, the integrand
// contains a peaked distribution. The peak is located using bisection and the integrand is
// split into two pieces, above and below the peak. An adaptive integration routine used for
// each part of the integrand.
type PdfZ struct {
	eps_quad, eps_bisect     float64
	alpha_tol, beta_tol      float64
	limit_quad, limit_bisect int
	quad                     quad.Integrator
}

// Creates a new generator for Levy-stable distribution, using an adaptive
// Gauss-Kronrod Integrator from the quad package.
//
// Default parameters are:
// 		eps_quad = integration target accuracy = 1e-14
//		eps_bisect = accuracy for locating max of integrand = 1e-10
//		limit_quad = integration iteration limit = 45
//		limit_bisect = bisection iteration limit = 50
func NewPdfZ() *PdfZ {

	eps_quad := 1e-14   // quadrature precision
	eps_bisect := 1e-10 // bisection precision
	alpha_tol := 1e-6   // tolerance for alpha to be close to special value (1 or 2)
	beta_tol := 1e-6    // tolerance for beta to be close to special value (0)
	limit_quad := 45    // quadrature iteration limit
	limit_bisect := 50  // bisection iteration limit
	q := quad.NewCustomQAG(quad.NewGaussKronrod61())
	return &PdfZ{eps_quad, eps_bisect, alpha_tol, beta_tol, limit_quad, limit_bisect, q}
}

// Creates a new generator for Levy-stable distribution, using a supplied quad.Integrator
// and parameters.
//
// Parameters are:
// 		eps_quad = integration target accuracy
//		eps_bisect = accuracy for locating max of integrand
//		limit_quad = integration iteration limit
//		limit_bisect = bisection iteration limit
func NewCustomPdfZ(q quad.Integrator, eps_quad, eps_bisect, alpha_tol, beta_tol float64, limit_quad, limit_bisect int) *PdfZ {
	return &PdfZ{eps_quad, eps_bisect, alpha_tol, beta_tol, limit_quad, limit_bisect, q}
}

func (pdf *PdfZ) ScaledValue(x, alpha, beta float64) (float64, error) {

	if closeTo(alpha, 2, pdf.alpha_tol) {

		// Gaussian case, for appropriately normalised levy distribution
		return math.Exp(-0.25*x*x) / math.Sqrt(4.0*math.Pi), nil

	} else if closeTo(alpha, 1, pdf.alpha_tol) && !closeTo(beta, 0, pdf.beta_tol) {

		// This tends to suffer from small oscillations in the integrated distribution
		// Need futher integration to sort this out

		gamma := math.Exp(-0.5 * math.Pi * x / beta)
		a := -0.5 * math.Pi
		b := 0.5 * math.Pi

		p, err := pdf.integrate(
			func(theta float64) float64 {
				return componentEq1(theta, beta)*gamma - 1.0
			},
			func(theta float64) float64 {
				return integrandEq1(theta, beta, gamma)
			},
			a, b)

		p *= 0.5 * gamma / math.Abs(beta)
		return p, err

	} else if closeTo(alpha, 1, pdf.alpha_tol) && closeTo(beta, 0, pdf.beta_tol) {

		// Cauchy distribution
		return 1.0 / ((1.0 + x*x) * math.Pi), nil

	} else if !closeTo(alpha, 1, pdf.alpha_tol) {

		zeta := -beta * math.Tan(0.5*math.Pi*alpha)

		if x == zeta {

			eps := math.Atan(-zeta) / alpha
			p := math.Gamma(1+1/alpha) * math.Cos(eps) /
				(math.Pi * math.Pow(1+zeta*zeta, 0.5/alpha))
			return p, nil
		} else if x > zeta {

			eps := math.Atan(-zeta) / alpha
			gamma := math.Pow(x-zeta, alpha/(alpha-1.0))
			a := -eps
			b := 0.5 * math.Pi

			p, err := pdf.integrate(
				func(theta float64) float64 {
					return componentNeq1(theta, alpha, beta, eps)*gamma - 1.0
				},
				func(theta float64) float64 {
					return integrandNeq1(theta, alpha, beta, eps, gamma)
				},
				a, b)

			p *= alpha * math.Pow(x-zeta, 1.0/(alpha-1.0)) / (math.Pi * math.Abs(alpha-1))
			return p, err
		} else if x < zeta {
			return pdf.ScaledValue(-x, alpha, -beta)
		}
	}
	return 0, nil
}

func (pdf *PdfZ) integrate(f func(theta float64) float64, g func(theta float64) float64, a, b float64) (float64, error) {

	// Find max in integrand, so that integration can capture peak in integrand
	// f is the function that provides location of maximum
	// g is the integrand

	// Find max of integrand by solving f = 0
	max, _, err := solve.Bisection(f, a, b, pdf.eps_bisect, pdf.limit_bisect)

	//	fmt.Printf("Got max at=%g val=%g\n",max,fmax)

	// Integration up to the peak
	I_1, _, e := pdf.quad.Integrate(g, a, max, pdf.eps_quad, 0.0, pdf.limit_quad)

	if e != nil {
		if err == nil {
			err = e
		} else {
			err = fmt.Errorf("%s; %s ", err, e)
		}
	}

	// Integration from the peak
	I_2, _, e := pdf.quad.Integrate(g, max, b, pdf.eps_quad, 0.0, pdf.limit_quad)

	if e != nil {
		if err == nil {
			err = e
		} else {
			err = fmt.Errorf("%s; %s ", err, e)
		}
	}

	return (I_1 + I_2), err
}

func componentNeq1(theta, alpha, beta, eps float64) float64 {

	// alpha != 1, beta != 0 case
	// These values are checked by calling function
	return math.Pow(math.Cos(alpha*eps), 1.0/(alpha-1.0)) *
		math.Pow(math.Cos(theta)/math.Sin(alpha*(theta+eps)), alpha/(alpha-1.0)) *
		math.Cos(alpha*eps+(alpha-1.0)*theta) / math.Cos(theta)
}

func integrandNeq1(theta, alpha, beta, eps, gamma float64) float64 {

	val := componentNeq1(theta, alpha, beta, eps)

	// If val = +Inf integrand =  0
	// If gamma = NaN (actually + Inf) then integrand = 0
	if math.IsNaN(gamma) || math.IsInf(val, 1) {
		return 0
	}

	//	return val*math.Exp(-val*gamma)
	return math.Exp(math.Log(val) - val*gamma)
}

func componentEq1(theta, beta float64) float64 {

	// alpha == 1 case, beta != 0 case
	// These values are checked by calling function
	return (1.0 + 2.0*beta*theta/math.Pi) *
		math.Exp((0.5*math.Pi/beta+theta)*math.Tan(theta)) /
		math.Cos(theta)
}

func integrandEq1(theta, beta, gamma float64) float64 {

	val := componentEq1(theta, beta)

	// If val = +Inf integrand =  0
	// If gamma = NaN (actually + Inf) then integrand = 0
	if math.IsNaN(gamma) || math.IsInf(val, 1) {
		return 0
	}

	//	return val*math.Exp(-val*gamma)
	return math.Exp(math.Log(val) - val*gamma)
}

func closeTo(val1, val2, tol float64) bool {
	return math.Abs(val1-val2) < math.Abs(tol)
}
