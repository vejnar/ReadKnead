//
// Copyright (C) 2017-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package operations

import (
	"fmt"

	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"

	"github.com/buger/jsonparser"
)

type Length struct {
	name      string
	minLength int
	maxLength int
}

func NewLength(data []byte) (*Length, error) {
	l := Length{name: "length"}
	// minLength
	minLength, err := jsonparser.GetInt(data, "min_length")
	if err == jsonparser.KeyPathNotFoundError {
		l.minLength = -1
	} else if err != nil {
		return &l, err
	} else {
		l.minLength = int(minLength)
	}
	// maxLength
	maxLength, err := jsonparser.GetInt(data, "max_length")
	if err == jsonparser.KeyPathNotFoundError {
		l.maxLength = -1
	} else if err != nil {
		return &l, err
	} else {
		l.maxLength = int(maxLength)
	}
	return &l, nil
}

func (op *Length) Name() string {
	return op.name
}

func (op *Length) IsThreadSafe() bool {
	return true
}

func (op *Length) GetDpx(idx int) ([][]byte, int) {
	return [][]byte{}, idx
}

func (op *Length) Transform(p *fastq.ExtPair, r int, ot *OpStat, verboseLevel int) int {
	if r == 1 {
		if op.minLength != -1 && len(p.R1.Seq) < op.minLength {
			if verboseLevel > 2 {
				fmt.Printf("length r1 too_short %d\n", len(p.R1.Seq))
			}
			ot.OpsR1[op.name]["too_short"]++
			return 1
		} else {
			return 0
		}
		if op.maxLength != -1 && len(p.R1.Seq) > op.maxLength {
			if verboseLevel > 2 {
				fmt.Printf("length r1 too_long %d\n", len(p.R1.Seq))
			}
			ot.OpsR1[op.name]["too_long"]++
			return 1
		} else {
			return 0
		}
	} else {
		if op.minLength != -1 && len(p.R2.Seq) < op.minLength {
			if verboseLevel > 2 {
				fmt.Printf("length r2 too_short %d\n", len(p.R2.Seq))
			}
			ot.OpsR2[op.name]["too_short"]++
			return 1
		} else {
			return 0
		}
		if op.maxLength != -1 && len(p.R2.Seq) > op.maxLength {
			if verboseLevel > 2 {
				fmt.Printf("length r2 too_long %d\n", len(p.R2.Seq))
			}
			ot.OpsR2[op.name]["too_long"]++
			return 1
		} else {
			return 0
		}
	}
	return 0
}
