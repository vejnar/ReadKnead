//
// Copyright © 2025 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package operations

import (
	"fmt"

	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"
	"git.sr.ht/~vejnar/ReadKnead/lib/param"

	"github.com/buger/jsonparser"
)

type Quality struct {
	name       string
	minQuality float32
	function   string
	param      param.Parameters
}

func NewQuality(data []byte, param param.Parameters) (*Quality, error) {
	q := Quality{name: "quality", param: param}
	// minQuality
	minQuality, err := jsonparser.GetFloat(data, "min_quality")
	if err == jsonparser.KeyPathNotFoundError {
		q.minQuality = -1
	} else if err != nil {
		return &q, err
	} else {
		q.minQuality = float32(minQuality)
	}
	// function
	function, err := jsonparser.GetString(data, "function")
	if err == jsonparser.KeyPathNotFoundError {
		q.function = "average"
	} else if err != nil {
		return &q, err
	} else {
		q.function = function
	}
	if !(q.function == "average") {
		return &q, fmt.Errorf("unknown function: %s", q.function)
	}
	return &q, nil
}

func (op *Quality) Name() string {
	return op.name
}

func (op *Quality) IsThreadSafe() bool {
	return true
}

func (op *Quality) GetDpx(idx int) ([][]byte, int) {
	return [][]byte{}, idx
}

func (op *Quality) Transform(p *fastq.ExtPair, r int, ot *OpStat, verboseLevel int) int {
	var read *fastq.Record
	if r == 1 {
		read = &p.R1
	} else {
		read = &p.R2
	}
	if verboseLevel > 2 {
		fmt.Printf("%s %s r%d\n%s\n%s\n", op.name, read.Name, r, read.Seq, read.Qual)
	}
	if op.function == "average" {
		sumQual := 0
		for i := 0; i < len(read.Qual); i++ {
			sumQual += int(read.Qual[i]) - op.param.AsciiMin
		}
		avgQual := float32(sumQual) / float32(len(read.Qual))
		if verboseLevel > 3 {
			fmt.Printf("> avg:%.2f discard:%t\n", avgQual, avgQual < op.minQuality)
		}
		if avgQual < op.minQuality {
			return 1
		} else {
			return 0
		}
	}
	return 0
}
