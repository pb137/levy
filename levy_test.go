package levy

import (
	"code.google.com/p/plotinum/plotter"
	"fmt"
	"github.com/pb137/plotutil"
	"github.com/pb137/quad"
	"testing"
)

func TestPrintPdfMethods(t *testing.T) {

	fs := make([]plotter.XYs, 4)
	labels := make([]string, 4)

	pdf1 := NewCustomPdfZ(quad.NewSimpsonsRule(), 1e-10, 1e-10, 1e-6, 1e-6, 15, 50)
	fs[0] = plotutil.CreatePlotData(func(x float64) float64 {
		p, _ := Value(pdf1, x, 0.5, 0.5, 0.0, 1.0)
		return p
	}, -5.0, 5.0, 100)
	labels[0] = "Simpson"

	//	pdf2 := NewCustomPdfZ(quad.NewQAG(), 1e-14, 1e-10, 45, 50)
	pdf2 := NewPdfZ()
	fs[1] = plotutil.CreatePlotData(func(x float64) float64 {
		p, _ := Value(pdf2, x, 0.5, 0.5, 0.0, 1.0)
		return p
	}, -5.0, 5.0, 100)
	labels[1] = "Default QAG - GK21"

	pdf3 := NewCustomPdfZ(quad.NewCustomQAG(quad.NewGaussKronrod61()), 1e-14, 1e-10, 1e-6, 1e-6, 45, 50)
	fs[2] = plotutil.CreatePlotData(func(x float64) float64 {
		p, _ := Value(pdf3, x, 0.5, 0.5, 0.0, 1.0)
		return p
	}, -5.0, 5.0, 100)
	labels[2] = "Custom QAG - GK61"

	pdf4 := NewPdfBelov()
	fs[3] = plotutil.CreatePlotData(func(x float64) float64 {
		p, _ := Value(pdf4, x, 0.5, 0.5, 0.0, 1.0)
		return p
	}, -5.0, 5.0, 100)
	labels[3] = "Two integrands"

	plotutil.CreatePlot(fs, labels, map[string]string{"Title": "alpha_0.5_integrand_comp"})
}

func TestRangeAlpha(t *testing.T) {

	fs := make([]plotter.XYs, 6)
	labels := make([]string, 6)

	alpha_min := 0.5
	dalpha := 0.25

	pdf := NewPdfZ()
	// pdf := NewCustomPdf(quad.NewCustomQAG(quad.NewGaussKronrod61()), 1e-14, 1e-10, 45, 50)
	//pdf := NewCustomPdf(quad.NewSimpsonsRule(), 1e-14, 1e-10, 15, 50)
	for i := 0; i < 6; i++ {

		alpha := alpha_min + float64(i)*dalpha
		fs[i] = plotutil.CreatePlotData(func(x float64) float64 {
			p, _ := Value(pdf, x, alpha, 0.5, 0.0, 1.0)
			return p
		}, -15.0, 15.0, 200)
		labels[i] = fmt.Sprintf("alpha = %g", alpha)
	}
	plotutil.CreatePlot(fs, labels, map[string]string{"Title": "alpha_plot_beta=0.5"})
}

func TestRangeAlpha2(t *testing.T) {

	fs := make([]plotter.XYs, 3)
	labels := make([]string, 3)

	pdf := NewPdfZ()
	alphas := [...]float64{1.0, 1.01, 1.01}
	for i, alpha := range alphas {
		fs[i] = plotutil.CreatePlotData(func(x float64) float64 {
			p, _ := Value(pdf, x, alpha, 0.5, 0.0, 1.0)
			return p
		}, -15.0, 15.0, 200)
		labels[i] = fmt.Sprintf("alpha = %g", alpha)
	}
	plotutil.CreatePlot(fs, labels, map[string]string{"Title": "alpha_1ish_plot_beta=0.5"})
}

func TestRangeBeta(t *testing.T) {

	fs := make([]plotter.XYs, 6)
	labels := make([]string, 6)

	alpha := 0.5
	beta_min := -0.75
	dbeta := 0.25

	pdf := NewPdfZ()
	for i := 0; i < 6; i++ {

		beta := beta_min + float64(i)*dbeta
		fs[i] = plotutil.CreatePlotData(func(x float64) float64 {
			p, _ := Value(pdf, x, alpha, beta, 0.0, 1.0)
			return p
		}, -15.0, 15.0, 200)
		labels[i] = fmt.Sprintf("alpha = %g", alpha)
	}
	plotutil.CreatePlot(fs, labels, map[string]string{"Title": "alpha_plot_alpha=0.5"})
}

func TestBergstrom(t *testing.T) {

	alpha := 0.75
	beta := 0.25
	pdfZ := NewPdfZ()
	pdfB := NewPdfBergstrom()
	for i := 1; i <= 200; i += 10 {
		x := float64(i)
		p, _ := Value(pdfZ, x, alpha, beta, 0.0, 1.0)
		pb, _ := Value(pdfB, x, alpha, beta, 0.0, 1.0)
		t.Logf("x=%g p=%g pb=%g\n", x, p, pb)
	}
}

func TestPrintBergstrom(t *testing.T) {

	fs := make([]plotter.XYs, 2)
	labels := make([]string, 2)

	alpha := 0.5
	beta := 0.25
	pdfZ := NewPdfZ()
	pdfB := NewPdfBergstrom()
	fs[0] = plotutil.CreateLogPlotData(func(x float64) float64 {
		p, _ := Value(pdfZ, x, alpha, beta, 0.0, 1.0)
		return p
	}, 10, 10000, 100)
	labels[0] = "Default QAG - GK21"

	fs[1] = plotutil.CreateLogPlotData(func(x float64) float64 {
		p, _ := Value(pdfB, x, alpha, beta, 0.0, 1.0)
		return p
	}, 10.0, 10000, 100)
	labels[1] = "Bergstrom"

	plotutil.CreatePlot(fs, labels, map[string]string{"Title": "bergstrom_alpha_05_beta_025"})
}

func TestRelError(t *testing.T) {

	fs := make([]plotter.XYs, 1)
	labels := make([]string, 1)

	alpha := 1.9
	beta := 0.25
	pdfZ := NewPdfZ()
	pdfB := NewPdfBergstrom()
	fs[0] = plotutil.CreatePlotData(func(x float64) float64 {
		p, _ := Value(pdfZ, x, alpha, beta, 0.0, 1.0)
		pb, _ := Value(pdfB, x, alpha, beta, 0.0, 1.0)
		return (p - pb) / p
	}, 10, 100, 20)
	labels[0] = "Relative Error"

	plotutil.CreatePlot(fs, labels, map[string]string{"Title": "bergstrom_rel_error_alpha_0.5"})
}
