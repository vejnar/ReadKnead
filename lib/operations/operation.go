//
// Copyright © 2017 Charles E. Vejnar
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

type Operation interface {
	Name() string
	IsThreadSafe() bool
	GetDpx(int) ([][]byte, int)
	Transform(*fastq.ExtPair, int, *OpStat, int) int
}

func ReadOps(data []byte, param param.Parameters) ([]Operation, error) {
	var ops []Operation
	var err error
	jsonparser.ArrayEach(data, func(value []byte, dataType jsonparser.ValueType, offset int, err2 error) {
		if err == nil {
			var op Operation
			var opName string
			opName, err = jsonparser.GetUnsafeString(value, "name")
			if err != nil {
				err = fmt.Errorf("Operation \"name\" missing")
				return
			}
			switch opName {
			case "clip":
				op, err = NewClip(value)
			case "demultiplex":
				op, err = NewDemultiplex(value)
			case "length":
				op, err = NewLength(value)
			case "quality":
				op, err = NewQuality(value, param)
			case "random":
				op, err = NewRandom(value)
			case "rename":
				op, err = NewRename(value)
			case "trim":
				op, err = NewTrim(value, param)
			default:
				err = fmt.Errorf("unknown operation: %s", opName)
			}
			if err != nil {
				return
			}
			ops = append(ops, op)
		}
	})
	return ops, err
}

func joinTwo(a []byte, b []byte) []byte {
	n := len(a) + len(b)
	o := make([]byte, n)
	copy(o, a)
	copy(o[len(a):], b)
	return o
}

func joinThree(a []byte, b []byte, c []byte) []byte {
	n := len(a) + len(b) + len(c)
	o := make([]byte, n)
	copy(o, a)
	copy(o[len(a):], b)
	copy(o[len(a)+len(b):], c)
	return o
}
