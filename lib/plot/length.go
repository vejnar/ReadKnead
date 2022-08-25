//
// Copyright (C) 2017-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//

package plot

import (
	"strconv"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/font"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

func BarplotLength(lengths map[int]uint64, path string, title string) error {
	// Change default font
	plot.DefaultFont = font.Font{Typeface: "Liberation", Variant: "Sans"}
	// Max read length
	maxLength := 0
	for l, _ := range lengths {
		if l > maxLength {
			maxLength = l
		}
	}
	// Prepare data
	data := make(plotter.Values, maxLength+1)
	for i := 0; i < len(data); i++ {
		if v, ok := lengths[i]; ok {
			data[i] = float64(v)
		}
	}
	// Lowest and highest lengths
	var xMin, xMax int
	for i := 0; i < len(data); i++ {
		if data[i] != 0 {
			xMin = i
			break
		}
	}
	xMin -= 2
	for i := len(data) - 1; i >= 0; i-- {
		if data[i] != 0 {
			xMax = i
			break
		}
	}
	xMax += 2

	// Start plot
	p := plot.New()
	// Add bars
	w := vg.Points(200 / float64(xMax-xMin))
	bars, err := plotter.NewBarChart(data, w)
	if err != nil {
		return err
	}
	bars.LineStyle.Width = vg.Length(0)
	bars.Color = plotutil.Color(0)
	p.Add(bars)
	// Legends
	p.X.Label.Text = "Read length"
	p.X.Min = float64(xMin)
	p.X.Max = float64(xMax)
	p.Y.Tick.Marker = commaTicks{}
	p.Title.Text = title

	// Save
	if err := p.Save(5*vg.Inch, 3*vg.Inch, path); err != nil {
		return err
	}
	return nil
}

type commaTicks struct{}

// Ticks computes the default tick marks, but inserts commas
// into the labels for the major tick marks.
func (commaTicks) Ticks(min, max float64) []plot.Tick {
	tks := plot.DefaultTicks{}.Ticks(min, max)
	for i, t := range tks {
		if t.Label == "" { // Skip minor ticks, they are fine.
			continue
		}
		tks[i].Label = addCommas(strconv.FormatFloat(t.Value, 'f', -1, 64))
	}
	return tks
}

// AddCommas adds commas after every 3 characters.
func addCommas(s string) string {
	if len(s) <= 3 {
		return s
	} else {
		return addCommas(s[0:len(s)-3]) + "," + s[len(s)-3:]
	}
}
