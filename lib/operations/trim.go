//
// Copyright © 2017 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package operations

import (
	"errors"
	"fmt"

	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"
	"git.sr.ht/~vejnar/ReadKnead/lib/param"
	"git.sr.ht/~vejnar/ReadKnead/lib/trim"

	"git.sr.ht/~vejnar/bktrim"
	"github.com/buger/jsonparser"
)

const (
	TrimAlignNW = iota
	TrimBKTrim
	TrimBKTrimPaired
	TrimSearch
	TrimMatch
	TrimQuality
)

type Trim struct {
	name               string
	end                int
	sequences          [][]byte
	sequencesPaired    [][]byte
	addTrimmed         bool
	addTrimmedRef      bool
	addSeparator       bool
	algo               int
	algoName           string
	minSequence        int
	minScore           float32
	position           int
	keep               []bool
	lengthLigand       int
	addLigand          bool
	addLigandSeparator bool
	applyTrimSeq       bool
	bkMatrices         []*bktrim.Matrix
	window             int
	unqualifiedPropMax float32
	minQuality         int
	param              param.Parameters
}

func NewTrim(data []byte, param param.Parameters) (*Trim, error) {
	t := Trim{name: "trim", param: param}
	sequence, err := jsonparser.GetUnsafeString(data, "sequence")
	if err != nil && err != jsonparser.KeyPathNotFoundError {
		return &t, err
	}
	if sequence != "" {
		t.sequences = append(t.sequences, []byte(sequence))
	}
	err = nil
	jsonparser.ArrayEach(data, func(value []byte, dataType jsonparser.ValueType, offset int, err2 error) {
		if err == nil {
			var s string
			s, err = jsonparser.ParseString(value)
			if err != nil {
				return
			}
			t.sequences = append(t.sequences, []byte(s))
		}
	}, "sequences")
	if err != nil {
		return &t, err
	}
	sequencePaired, err := jsonparser.GetUnsafeString(data, "sequence_paired")
	if err != nil && err != jsonparser.KeyPathNotFoundError {
		return &t, err
	}
	if sequencePaired != "" {
		t.sequencesPaired = append(t.sequencesPaired, []byte(sequencePaired))
	}
	err = nil
	jsonparser.ArrayEach(data, func(value []byte, dataType jsonparser.ValueType, offset int, err2 error) {
		if err == nil {
			var s string
			s, err = jsonparser.ParseString(value)
			if err != nil {
				return
			}
			t.sequencesPaired = append(t.sequencesPaired, []byte(s))
		}
	}, "sequences_paired")
	if err != nil {
		return &t, err
	}
	addTrimmed, err := jsonparser.GetBoolean(data, "add_trimmed")
	if err == jsonparser.KeyPathNotFoundError {
		t.addTrimmed = false
	} else if err != nil {
		return &t, err
	} else {
		t.addTrimmed = addTrimmed
	}
	addTrimmedRef, err := jsonparser.GetBoolean(data, "add_trimmed_ref")
	if err == jsonparser.KeyPathNotFoundError {
		addTrimmedRef = false
	} else if err != nil {
		return &t, err
	} else {
		t.addTrimmedRef = addTrimmedRef
	}
	addSeparator, err := jsonparser.GetBoolean(data, "add_separator")
	if err == jsonparser.KeyPathNotFoundError {
		addSeparator = true
	} else if err != nil {
		return &t, err
	} else {
		t.addSeparator = addSeparator
	}
	algoRaw, err := jsonparser.GetString(data, "algo")
	if err == jsonparser.KeyPathNotFoundError {
		if param.Paired {
			algoRaw = "bktrim_paired"
		} else {
			algoRaw = "bktrim"
		}
	} else if err != nil {
		return &t, err
	}
	if !(algoRaw == "bktrim" || algoRaw == "bktrim_paired" || algoRaw == "align" || algoRaw == "search" || algoRaw == "match" || algoRaw == "trimqual") {
		return &t, fmt.Errorf("Unknown trimming algorithm: %s", algoRaw)
	}
	if algoRaw != "trimqual" && len(t.sequences) == 0 {
		return &t, fmt.Errorf("Sequence to trim not found")
	}
	if algoRaw != "bktrim_paired" {
		end, err := jsonparser.GetInt(data, "end")
		if err != nil {
			if errors.Is(err, jsonparser.KeyPathNotFoundError) {
				return &t, fmt.Errorf("%w: %v", err, "end")
			}
			return &t, err
		} else {
			t.end = int(end)
		}
	}
	if algoRaw == "align" {
		t.algo = TrimAlignNW
		t.algoName = algoRaw
	} else if algoRaw == "bktrim" || algoRaw == "bktrim_paired" {
		epsilon, err := jsonparser.GetFloat(data, "epsilon")
		if err == jsonparser.KeyPathNotFoundError {
			epsilon = 0.1
		} else if err != nil {
			return &t, err
		}
		epsilonIndel, err := jsonparser.GetFloat(data, "epsilon_indel")
		if err == jsonparser.KeyPathNotFoundError {
			epsilonIndel = 0.03
		} else if err != nil {
			return &t, err
		}
		minOverlap, err := jsonparser.GetInt(data, "min_overlap")
		if err == jsonparser.KeyPathNotFoundError {
			minOverlap = 3
		} else if err != nil {
			return &t, err
		}
		if algoRaw == "bktrim" {
			t.algo = TrimBKTrim
			t.algoName = algoRaw
			t.bkMatrices = trim.NewMatrixAdapter(t.sequences, t.end, epsilon, epsilonIndel, int(minOverlap), param.AsciiMin)
		} else if algoRaw == "bktrim_paired" {
			if len(t.sequencesPaired) == 0 {
				return &t, fmt.Errorf("Sequence to trim on paired read not found")
			}
			t.algo = TrimBKTrimPaired
			t.algoName = algoRaw
			t.bkMatrices = trim.NewMatrixAdapterPaired(t.sequences, t.sequencesPaired, epsilon, epsilonIndel, int(minOverlap), param.AsciiMin)
		}
	} else if algoRaw == "search" {
		t.algo = TrimSearch
		t.algoName = algoRaw
	} else if algoRaw == "match" {
		t.algo = TrimMatch
		t.algoName = algoRaw
	} else if algoRaw == "trimqual" {
		t.algo = TrimQuality
		t.algoName = algoRaw
	}
	minSequence, err := jsonparser.GetInt(data, "min_sequence")
	if err == jsonparser.KeyPathNotFoundError {
		t.minSequence = 0
	} else if err != nil {
		return &t, err
	} else {
		t.minSequence = int(minSequence)
	}
	minScore, err := jsonparser.GetFloat(data, "min_score")
	if err == jsonparser.KeyPathNotFoundError {
		if algoRaw == "align" {
			t.minScore = float32(len(t.sequences[0]) * 5.)
		} else {
			t.minScore = 0.8
		}
	} else if err != nil {
		return &t, err
	} else {
		t.minScore = float32(minScore)
	}
	position, err := jsonparser.GetInt(data, "position")
	if err != nil && err != jsonparser.KeyPathNotFoundError {
		return &t, err
	} else {
		t.position = int(position)
	}
	err = nil
	t.keep = make([]bool, len(trim.TrimTypes))
	found := false
	jsonparser.ArrayEach(data, func(value []byte, dataType jsonparser.ValueType, offset int, err2 error) {
		if err == nil {
			var k string
			k, err = jsonparser.ParseString(value)
			if err != nil {
				return
			}
			for itt, tt := range trim.TrimTypes {
				if k == tt.String() {
					t.keep[itt] = true
				}
			}
			found = true
		}
	}, "keep")
	if !found {
		for i := range t.keep {
			t.keep[i] = true
		}
	}
	if err != nil {
		return &t, err
	}
	lengthLigand, err := jsonparser.GetInt(data, "length_ligand")
	if err == jsonparser.KeyPathNotFoundError {
		t.lengthLigand = 0
	} else if err != nil {
		return &t, err
	} else {
		t.lengthLigand = int(lengthLigand)
	}
	addLigand, err := jsonparser.GetBoolean(data, "add_ligand")
	if err == jsonparser.KeyPathNotFoundError {
		t.addLigand = false
	} else if err != nil {
		return &t, err
	} else {
		t.addLigand = addLigand
	}
	addLigandSeparator, err := jsonparser.GetBoolean(data, "add_ligand_separator")
	if err == jsonparser.KeyPathNotFoundError {
		t.addLigandSeparator = true
	} else if err != nil {
		return &t, err
	} else {
		t.addLigandSeparator = addLigandSeparator
	}
	applyTrimSeq, err := jsonparser.GetBoolean(data, "apply_trim_seq")
	if err == jsonparser.KeyPathNotFoundError {
		t.applyTrimSeq = true
	} else if err != nil {
		return &t, err
	} else {
		t.applyTrimSeq = applyTrimSeq
	}
	window, err := jsonparser.GetInt(data, "window")
	if err == jsonparser.KeyPathNotFoundError {
		t.window = 5
	} else if err != nil {
		return &t, err
	} else {
		t.window = int(window)
	}
	unqualifiedPropMax, err := jsonparser.GetFloat(data, "unqualified_prop_max")
	if err == jsonparser.KeyPathNotFoundError {
		t.unqualifiedPropMax = 0.6
	} else if err != nil {
		return &t, err
	} else {
		t.unqualifiedPropMax = float32(unqualifiedPropMax)
	}
	minQuality, err := jsonparser.GetInt(data, "min_quality")
	if err == jsonparser.KeyPathNotFoundError {
		t.minQuality = 15
	} else if err != nil {
		return &t, err
	} else {
		t.minQuality = int(minQuality)
	}
	return &t, nil
}

func (op *Trim) Name() string {
	return op.name
}

func (op *Trim) IsThreadSafe() bool {
	return true
}

func (op *Trim) GetDpx(idx int) ([][]byte, int) {
	return [][]byte{}, idx
}

func (op *Trim) Transform(p *fastq.ExtPair, r int, ot *OpStat, verboseLevel int) int {
	var trimIdx int
	var trimType trim.TrimType
	var trimScore float32
	var trimSeq []byte
	if r == 1 {
		if verboseLevel > 2 {
			fmt.Printf("%s %s %s r%d\n%s\n", op.name, op.algoName, p.R1.Name, r, p.R1.Seq)
		}
		switch op.algo {
		case TrimAlignNW:
			trimType, trimIdx, trimScore, trimSeq = trim.TrimAlign(&p.R1, op.sequences, op.minSequence, op.minScore, op.end, op.applyTrimSeq, verboseLevel)
		case TrimBKTrim:
			trimType, trimIdx, trimScore, trimSeq = trim.TrimBKTrim(&p.R1, op.bkMatrices, op.minSequence, op.end, op.applyTrimSeq, verboseLevel)
		case TrimBKTrimPaired:
			trimType, trimIdx, trimScore, trimSeq = trim.TrimBKTrimPaired(p, op.bkMatrices, op.applyTrimSeq, verboseLevel)
		case TrimSearch:
			trimType, trimIdx, trimScore, trimSeq = trim.TrimSearch(&p.R1, op.sequences, op.minSequence, op.minScore, op.end, op.applyTrimSeq, verboseLevel)
		case TrimMatch:
			trimType, trimIdx, trimScore, trimSeq = trim.TrimMatch(&p.R1, op.sequences, op.position, op.minSequence, op.minScore, op.end, op.applyTrimSeq, verboseLevel)
		case TrimQuality:
			trimType, trimIdx, trimScore, trimSeq = trim.TrimQuality(&p.R1, op.window, op.unqualifiedPropMax, op.minQuality, op.param.AsciiMin, op.end, op.applyTrimSeq, verboseLevel)
		}
		if verboseLevel > 2 {
			fmt.Printf("%s %s length:%d score:%.2f\n", p.R1.Seq, trimType, len(p.R1.Seq), trimScore)
		}
		ot.OpsR1[op.name][trimType.String()]++
		if op.keep[trimType] {
			// Add trimmed sequence
			if len(trimSeq) > 0 {
				if op.addTrimmedRef {
					if op.addSeparator {
						p.R1.Name = joinThree(p.R1.Name, []byte("#"), op.sequences[trimIdx])
						p.R2.Name = joinThree(p.R2.Name, []byte("#"), op.sequences[trimIdx])
					} else {
						p.R1.Name = joinTwo(p.R1.Name, op.sequences[trimIdx])
						p.R2.Name = joinTwo(p.R2.Name, op.sequences[trimIdx])
					}
				}
				if op.addTrimmed {
					if op.addSeparator {
						p.R1.Name = joinThree(p.R1.Name, []byte("#"), trimSeq)
						p.R2.Name = joinThree(p.R2.Name, []byte("#"), trimSeq)
					} else {
						p.R1.Name = joinTwo(p.R1.Name, trimSeq)
						p.R2.Name = joinTwo(p.R2.Name, trimSeq)
					}
				}
			}
			// Ligand
			if trimType != trim.NoTrimType && op.lengthLigand > 0 {
				pc := Clip{name: op.name, end: op.end, length: op.lengthLigand, addClipped: op.addLigand, addSeparator: op.addLigandSeparator}
				return pc.Transform(p, 1, ot, verboseLevel)
			}
			return 0
		} else {
			return 1
		}
	} else {
		if verboseLevel > 2 {
			fmt.Printf("%s %d %s r%d\n%s\n", op.name, op.algo, p.R2.Name, r, p.R2.Seq)
		}
		switch op.algo {
		case TrimAlignNW:
			trimType, trimIdx, trimScore, trimSeq = trim.TrimAlign(&p.R2, op.sequences, op.minSequence, op.minScore, op.end, op.applyTrimSeq, verboseLevel)
		case TrimBKTrim:
			trimType, trimIdx, trimScore, trimSeq = trim.TrimBKTrim(&p.R2, op.bkMatrices, op.minSequence, op.end, op.applyTrimSeq, verboseLevel)
		case TrimSearch:
			trimType, trimIdx, trimScore, trimSeq = trim.TrimSearch(&p.R2, op.sequences, op.minSequence, op.minScore, op.end, op.applyTrimSeq, verboseLevel)
		case TrimMatch:
			trimType, trimIdx, trimScore, trimSeq = trim.TrimMatch(&p.R2, op.sequences, op.position, op.minSequence, op.minScore, op.end, op.applyTrimSeq, verboseLevel)
		case TrimQuality:
			trimType, trimIdx, trimScore, trimSeq = trim.TrimQuality(&p.R2, op.window, op.unqualifiedPropMax, op.minQuality, op.param.AsciiMin, op.end, op.applyTrimSeq, verboseLevel)
		}
		if verboseLevel > 2 {
			fmt.Printf("%s %s length:%d score:%.2f\n", p.R2.Seq, trimType, len(p.R2.Seq), trimScore)
		}
		ot.OpsR2[op.name][trimType.String()]++
		if op.keep[trimType] {
			// Add trimmed sequence
			if len(trimSeq) > 0 {
				if op.addTrimmedRef {
					if op.addSeparator {
						p.R1.Name = joinThree(p.R1.Name, []byte("#"), op.sequences[trimIdx])
						p.R2.Name = joinThree(p.R2.Name, []byte("#"), op.sequences[trimIdx])
					} else {
						p.R1.Name = joinTwo(p.R1.Name, op.sequences[trimIdx])
						p.R2.Name = joinTwo(p.R2.Name, op.sequences[trimIdx])
					}
				}
				if op.addTrimmed {
					if op.addSeparator {
						p.R1.Name = joinThree(p.R1.Name, []byte("#"), trimSeq)
						p.R2.Name = joinThree(p.R2.Name, []byte("#"), trimSeq)
					} else {
						p.R1.Name = joinTwo(p.R1.Name, trimSeq)
						p.R2.Name = joinTwo(p.R2.Name, trimSeq)
					}
				}
			}
			// Ligand
			if trimType != trim.NoTrimType && op.lengthLigand > 0 {
				pc := Clip{name: op.name, end: op.end, length: op.lengthLigand, addClipped: op.addLigand, addSeparator: op.addLigandSeparator}
				return pc.Transform(p, 2, ot, verboseLevel)
			}
			return 0
		} else {
			return 1
		}
	}
	return 0
}
