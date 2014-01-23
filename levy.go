// Package levy implements routines for generating samples and pdfs from levy stable distributions.
// The pdf is calculated via direct numerical integration using adaptive quadrature as described in
// 		Borak, Szymon; Härdle, Wolfgang Karl; Weron, Rafał (2005) :
//		Stable distributions, SFB 649 discussion paper, No. 2005,008, http://hdl.handle.net/10419/25027.
package levy

import (
	"fmt"
)

// Interface Pdf generates value from a scaled Levy-stable distribution (ie mu=0, sigma=1.0).
type Pdf interface {
	ScaledValue(x, alpha, beta float64) (float64, error)
}

// Value returns the pdf of Levy-stable distribution.
func Value(pdf Pdf, x, alpha, beta, mu, sigma float64) (float64, error) {

	if alpha <= 0 || alpha > 2 {
		return 0, fmt.Errorf("alpha (%g) outside allowed range [0,2)", alpha)
	}

	if beta < -1 || beta > 1 {
		return 0, fmt.Errorf("beta (%g) outside allowed range [-1,1]", beta)
	}

	x = (x - mu) / sigma
	val, err := pdf.ScaledValue(x, alpha, beta)
	val /= sigma
	return val, err
}
