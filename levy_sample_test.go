package levy

import (
	"code.google.com/p/plotinum/plotter"
	"github.com/pb137/plotutil"
	"math"
	"math/rand"
	"testing"
)

func TestPrintSamples(t *testing.T) {

	fs := make([]plotter.XYs, 4)
	labels := make([]string, 4)

	fs[0] = GetSamples(func() float64 {
		return Sample(1.1, 0.0, 0.0, 1.0)
	}, 100)
	labels[0] = "alpha = 1.1"

	fs[1] = GetSamples(func() float64 {
		return Sample(1.5, 0.0, 0.0, 1.0)
	}, 100)
	labels[1] = "alpha = 1.5"

	fs[2] = GetSamples(func() float64 {
		return Sample(1.9, 0.0, 0.0, 1.0)
	}, 100)
	labels[2] = "alpha = 1.9"

	fs[3] = GetSamples(func() float64 {
		return Sample(2.0, 0.0, 0.0, 1.0)
	}, 100)
	labels[3] = "alpha = 2.0"

	plotutil.CreatePlot(fs, labels, map[string]string{"Title": "Samples"})
}

func TestPrintSamples2D(t *testing.T) {

	fs := make([]plotter.XYs, 4)
	labels := make([]string, 4)

	fs[0] = GetSamples2D(func() float64 {
		return Sample(1.1, 0.0, 0.0, 1.0)
	}, 500)
	labels[0] = "alpha = 1.1"

	fs[1] = GetSamples2D(func() float64 {
		return Sample(1.5, 0.0, 0.0, 1.0)
	}, 500)
	labels[1] = "alpha = 1.5"

	fs[2] = GetSamples2D(func() float64 {
		return Sample(1.9, 0.0, 0.0, 1.0)
	}, 500)
	labels[2] = "alpha = 1.9"

	fs[3] = GetSamples2D(func() float64 {
		return Sample(2.0, 0.0, 0.0, 1.0)
	}, 500)
	labels[3] = "alpha = 2.0"
	plotutil.CreatePlot(fs, labels, map[string]string{"Title": "Samples2D", "FixAspectRatio": "True"})
}

func GetSamples(f func() float64, N int) plotter.XYs {

	pts := make(plotter.XYs, N)
	pts[0].X = float64(0)
	pts[0].Y = f()
	for i := 1; i < N; i++ {
		pts[i].X = float64(i)
		pts[i].Y = pts[i-1].Y + f()
	}
	return pts
}

func GetSamples2D(f func() float64, N int) plotter.XYs {

	pts := make(plotter.XYs, N)
	pts[0].X = 0.0
	pts[0].Y = 0.0
	for i := 1; i < N; i++ {
		r := f()
		theta := 2.0 * math.Pi * rand.Float64()
		pts[i].X = pts[i-1].X + r*math.Cos(theta)
		pts[i].Y = pts[i-1].Y + r*math.Sin(theta)
	}
	return pts
}
