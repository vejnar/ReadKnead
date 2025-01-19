//
// Copyright © 2017 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package operations

import (
	"encoding/json"
	"fmt"
	"os"

	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"
	"git.sr.ht/~vejnar/ReadKnead/lib/plot"
)

type OpStat struct {
	OpsR1, OpsR2                                                 map[string]map[string]uint64
	KeptPair, TotalPair                                          uint64
	qualsInR1, qualsInR2, qualsOutR1, qualsOutR2                 [][]uint64
	lengthsInR1, lengthsInR2, lengthsOutR1, lengthsOutR2         map[int]uint64
	maxLengthInR1, maxLengthInR2, maxLengthOutR1, maxLengthOutR2 int
	statsInPath, statsOutPath, ReportPath, Label                 string
	paired                                                       bool
	asciiMin                                                     int
}

func initQL(maxReadLength int, maxQual int) (quals [][]uint64, lengths map[int]uint64) {
	// Quality per base
	// 2D array (maxReadLength x maxQual)
	quals = make([][]uint64, maxReadLength)
	cols := make([]uint64, maxReadLength*maxQual)
	for i := range quals {
		quals[i] = cols[i*maxQual : (i+1)*maxQual]
	}
	// Read length
	lengths = make(map[int]uint64)
	return
}

func updateQL(quals [][]uint64, lengths map[int]uint64, qualsAdd [][]uint64, lengthsAdd map[int]uint64) {
	// Quality
	for i := 0; i < len(quals); i++ {
		for j := 0; j < len(quals[i]); j++ {
			quals[i][j] += qualsAdd[i][j]
		}
	}
	// Length
	for k, v := range lengthsAdd {
		lengths[k] += v
	}
}

func writeQL(statsPath string, pairName string, quals [][]uint64, lengths map[int]uint64) error {
	// Quality
	if f, err := os.Create(statsPath + pairName + "_qual.txt"); err != nil {
		return err
	} else {
		var outString string
		for _, l := range quals {
			outString = fmt.Sprint(l)
			f.WriteString(outString[1:len(outString)-1] + "\n")
		}
		f.Close()
	}
	// Length
	if f, err := os.Create(statsPath + pairName + "_length.txt"); err != nil {
		return err
	} else {
		kmin, kmax := minMaxKeys(lengths)
		for k := kmin; k <= kmax; k++ {
			if v, ok := lengths[k]; ok {
				fmt.Fprintf(f, "%d\t%d\n", k, v)
			}
		}
		f.Close()
	}
	return nil
}

func minMaxKeys(m map[int]uint64) (min, max int) {
	for k := range m {
		if k < min {
			min = k
		}
		if k > max {
			max = k
		}
	}
	return
}

func NewOpStat(statsInPath string, statsOutPath string, reportPath string, label string, maxReadLength int, maxQual int, asciiMin int, paired bool, opsR1 []Operation, opsR2 []Operation) *OpStat {
	ot := OpStat{paired: paired, statsInPath: statsInPath, statsOutPath: statsOutPath, ReportPath: reportPath, Label: label, asciiMin: asciiMin}
	// Init. for statistics: In
	if ot.statsInPath != "" {
		// Read1
		ot.qualsInR1, ot.lengthsInR1 = initQL(maxReadLength, maxQual)
		// Read2
		if ot.paired {
			ot.qualsInR2, ot.lengthsInR2 = initQL(maxReadLength, maxQual)
		}
	}
	// Init. for statistics: Out
	if ot.statsOutPath != "" {
		ot.qualsOutR1, ot.lengthsOutR1 = initQL(maxReadLength, maxQual)
		// Read2
		if ot.paired {
			ot.qualsOutR2, ot.lengthsOutR2 = initQL(maxReadLength, maxQual)
		}
	}
	// Init. operations stats
	ot.OpsR1 = make(map[string]map[string]uint64)
	for _, op := range opsR1 {
		ot.OpsR1[op.Label()] = make(map[string]uint64)
	}
	ot.OpsR2 = make(map[string]map[string]uint64)
	for _, op := range opsR2 {
		ot.OpsR2[op.Label()] = make(map[string]uint64)
	}
	return &ot
}

func (ot *OpStat) CountIn(p *fastq.ExtPair) {
	if ot.statsInPath != "" {
		// Statistics: In, Quality
		for i, q := range p.R1.Qual {
			ot.qualsInR1[i][int(q)-ot.asciiMin] += 1
		}
		// Statistics: In, Length
		ot.lengthsInR1[len(p.R1.Seq)]++
		// Max read length
		if len(p.R1.Seq) > ot.maxLengthInR1 {
			ot.maxLengthInR1 = len(p.R1.Seq)
		}
		if ot.paired {
			// Statistics: In, Quality
			for i, q := range p.R2.Qual {
				ot.qualsInR2[i][int(q)-ot.asciiMin] += 1
			}
			// Statistics: In, Length
			ot.lengthsInR2[len(p.R2.Seq)]++
			// Max read length
			if len(p.R2.Seq) > ot.maxLengthInR2 {
				ot.maxLengthInR2 = len(p.R2.Seq)
			}
		}
	}
}

func (ot *OpStat) CountOut(p *fastq.ExtPair) {
	if ot.statsOutPath != "" {
		// Statistics: Out, Quality
		for i, q := range p.R1.Qual {
			ot.qualsOutR1[i][int(q)-ot.asciiMin] += 1
		}
		// Statistics: Out, Length
		ot.lengthsOutR1[len(p.R1.Seq)]++
		// Max read length
		if len(p.R1.Seq) > ot.maxLengthOutR1 {
			ot.maxLengthOutR1 = len(p.R1.Seq)
		}
		if ot.paired {
			// Statistics: Out, Quality
			for i, q := range p.R2.Qual {
				ot.qualsOutR2[i][int(q)-ot.asciiMin] += 1
			}
			// Statistics: Out, Length
			ot.lengthsOutR2[len(p.R2.Seq)]++
			// Max read length
			if len(p.R2.Seq) > ot.maxLengthOutR2 {
				ot.maxLengthOutR2 = len(p.R2.Seq)
			}
		}
	}
}

func (ot *OpStat) Update(otn *OpStat) {
	// Quality & Length: In
	if ot.statsInPath != "" {
		// Read1
		updateQL(ot.qualsInR1, ot.lengthsInR1, otn.qualsInR1, otn.lengthsInR1)
		// Read2
		if ot.paired {
			updateQL(ot.qualsInR2, ot.lengthsInR2, otn.qualsInR2, otn.lengthsInR2)
		}
	}
	// Quality & Length: Out
	if ot.statsOutPath != "" {
		updateQL(ot.qualsOutR1, ot.lengthsOutR1, otn.qualsOutR1, otn.lengthsOutR1)
		// Read2
		if ot.paired {
			updateQL(ot.qualsOutR2, ot.lengthsOutR2, otn.qualsOutR2, otn.lengthsOutR2)
		}
	}
	// Operations stats
	for op := range otn.OpsR1 {
		for name, v := range otn.OpsR1[op] {
			ot.OpsR1[op][name] += v
		}
	}
	for op := range otn.OpsR2 {
		for name, v := range otn.OpsR2[op] {
			ot.OpsR2[op][name] += v
		}
	}
	// Counts
	ot.KeptPair += otn.KeptPair
	ot.TotalPair += otn.TotalPair
}

func (ot *OpStat) Write() error {
	// Label separator
	var sep string
	var err error
	if ot.Label != "" {
		sep = " "
	}
	// Statistics: In
	if ot.statsInPath != "" {
		// Raw & Plot
		err = writeQL(ot.statsInPath, "_r1", ot.qualsInR1[:ot.maxLengthInR1], ot.lengthsInR1)
		if err != nil {
			return err
		}
		err = plot.BoxplotQuality(ot.qualsInR1[:ot.maxLengthInR1], ot.statsInPath+"_r1_qual.pdf", ot.Label+sep+"Read1 input")
		if err != nil {
			return err
		}
		err = plot.BarplotLength(ot.lengthsInR1, ot.statsInPath+"_r1_length.pdf", ot.Label+sep+"Read1 input")
		if err != nil {
			return err
		}
		if ot.paired {
			err = writeQL(ot.statsInPath, "_r2", ot.qualsInR2[:ot.maxLengthInR2], ot.lengthsInR2)
			if err != nil {
				return err
			}
			err = plot.BoxplotQuality(ot.qualsInR2[:ot.maxLengthInR2], ot.statsInPath+"_r2_qual.pdf", ot.Label+sep+"Read2 input")
			if err != nil {
				return err
			}
			err = plot.BarplotLength(ot.lengthsInR2, ot.statsInPath+"_r2_length.pdf", ot.Label+sep+"Read2 input")
			if err != nil {
				return err
			}
		}
	}
	// Statistics: Out
	if ot.statsOutPath != "" {
		if ot.maxLengthOutR1 > 0 {
			err = writeQL(ot.statsOutPath, "_r1", ot.qualsOutR1[:ot.maxLengthOutR1], ot.lengthsOutR1)
			if err != nil {
				return err
			}
			err = plot.BoxplotQuality(ot.qualsOutR1[:ot.maxLengthOutR1], ot.statsOutPath+"_r1_qual.pdf", ot.Label+sep+"Read1 output")
			if err != nil {
				return err
			}
			err = plot.BarplotLength(ot.lengthsOutR1, ot.statsOutPath+"_r1_length.pdf", ot.Label+sep+"Read1 output")
			if err != nil {
				return err
			}
		}
		if ot.paired {
			if ot.maxLengthOutR2 > 0 {
				err = writeQL(ot.statsOutPath, "_r2", ot.qualsOutR2[:ot.maxLengthOutR2], ot.lengthsOutR2)
				if err != nil {
					return err
				}
				err = plot.BoxplotQuality(ot.qualsOutR2[:ot.maxLengthOutR2], ot.statsOutPath+"_r2_qual.pdf", ot.Label+sep+"Read2 output")
				if err != nil {
					return err
				}
				err = plot.BarplotLength(ot.lengthsOutR2, ot.statsOutPath+"_r2_length.pdf", ot.Label+sep+"Read2 output")
				if err != nil {
					return err
				}
			}
		}
	}
	// Report
	if ot.ReportPath != "" {
		// Remove step(s) without any stats
		for k, v := range ot.OpsR1 {
			if len(v) == 0 {
				delete(ot.OpsR1, k)
			}
		}
		for k, v := range ot.OpsR2 {
			if len(v) == 0 {
				delete(ot.OpsR2, k)
			}
		}
		// Merge read1, read2 and global stats to final report
		statFinal := make(map[string]map[string]map[string]uint64)
		statFinal["read1"] = ot.OpsR1
		if ot.paired {
			statFinal["read2"] = ot.OpsR2
		}
		statFinal["pair"] = make(map[string]map[string]uint64)
		statFinal["pair"]["all"] = make(map[string]uint64)
		statFinal["pair"]["all"]["output"] = ot.KeptPair
		statFinal["pair"]["all"]["input"] = ot.TotalPair
		// JSON
		report, _ := json.MarshalIndent(statFinal, "", "  ")
		if ot.ReportPath != "-" {
			if f, err := os.Create(ot.ReportPath); err != nil {
				return err
			} else {
				f.Write(report)
				f.Close()
			}
		} else {
			fmt.Println(string(report))
		}
	}
	return nil
}
