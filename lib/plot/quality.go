//
// Copyright (C) 2017-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//

package plot

import (
	"image/color"
	"math"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/font"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

func BoxplotQuality(quals [][]uint64, path string, title string) error {
	// Change default font
	plot.DefaultFont = font.Font{Typeface: "Liberation", Variant: "Sans"}
	// Max read length
	maxLengthIdx := 0
	onlyZero := false
	for i := 0; i < len(quals); i++ {
		onlyZero = true
		for j := 0; j < len(quals[i]); j++ {
			if quals[i][j] > 0 {
				onlyZero = false
				break
			}
		}
		if onlyZero {
			break
		} else {
			maxLengthIdx = i
		}
	}

	// Start plot
	p := plot.New()
	// Quality thresholds: 1%
	l1, err := plotter.NewLine(hLine(0.5, float64(maxLengthIdx)+1.5, -10.*math.Log10(0.01)))
	if err != nil {
		return err
	}
	l1.LineStyle.Width = vg.Points(1.)
	l1.LineStyle.Color = color.RGBA{R: 0, G: 170, B: 0}
	p.Add(l1)
	// Quality thresholds: 5%
	l2, err := plotter.NewLine(hLine(0.5, float64(maxLengthIdx)+1.5, -10.*math.Log10(0.05)))
	if err != nil {
		return err
	}
	l2.LineStyle.Width = vg.Points(1.)
	l2.LineStyle.Color = color.RGBA{R: 255, G: 140, B: 0}
	p.Add(l2)
	// Add bars
	var yMax float64
	w := vg.Points(200 / float64(maxLengthIdx+1))
	for i := 0; i <= maxLengthIdx; i++ {
		b := new(plotter.BoxPlot)
		b.Location = float64(i + 1)
		b.Width = w
		b.CapWidth = 3 * w / 4
		// Trick: Ignoring the 0 quality (as it can also means no sequence)
		b.Median = float64(Quantile(quals[i][1:], 0.5))
		b.Quartile1 = float64(Quantile(quals[i][1:], 0.25))
		b.Quartile3 = float64(Quantile(quals[i][1:], 0.75))
		b.AdjLow = float64(Quantile(quals[i][1:], 0.035))
		b.AdjHigh = float64(Quantile(quals[i][1:], 0.965))
		b.Min = b.AdjLow
		b.Max = b.AdjHigh
		b.MedianStyle = draw.LineStyle{Color: color.RGBA{R: 255, G: 43, B: 15}, Width: vg.Points(0.7)}
		b.GlyphStyle = plotter.DefaultGlyphStyle
		b.BoxStyle = draw.LineStyle{Color: color.Gray{90}, Width: vg.Points(0.7)}
		b.WhiskerStyle = draw.LineStyle{Color: color.Gray{180}, Width: vg.Points(0.4)}
		p.Add(b)
		if yMax < b.Max {
			yMax = b.Max
		}
	}
	// Legends
	p.X.Max = float64(maxLengthIdx) + 2.
	p.Y.Min = 0.
	p.Y.Max = yMax + 1.
	p.X.Label.Text = "Read position"
	p.Y.Label.Text = "Quality score"
	p.Title.Text = title

	// Save
	if err := p.Save(5*vg.Inch, 3*vg.Inch, path); err != nil {
		return err
	}
	return nil
}

func Quantile(data []uint64, q float64) int {
	// Factor
	var factor float64
	for i := 0; i < len(data); i++ {
		factor += float64(data[i])
	}
	factor *= q
	// Cumulative sum
	cumsum := make([]float64, len(data))
	for i := 1; i < len(data); i++ {
		cumsum[i] = cumsum[i-1] + float64(data[i])
	}
	// Apply factor
	for i := 0; i < len(data); i++ {
		cumsum[i] = math.Abs(cumsum[i] - factor)
	}
	// Find quantile
	vmin := cumsum[0]
	imin := 0
	for i := 0; i < len(cumsum); i++ {
		if cumsum[i] < vmin {
			vmin = cumsum[i]
			imin = i
		}
	}
	return imin + 1
}

func hLine(start, end, y float64) plotter.XYs {
	pts := make(plotter.XYs, 2)
	pts[0].X = start
	pts[0].Y = y
	pts[1].X = end
	pts[1].Y = y
	return pts
}
