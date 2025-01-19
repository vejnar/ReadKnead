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
	"strconv"

	"git.sr.ht/~vejnar/ReadKnead/lib/bio"
	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"

	"github.com/buger/jsonparser"
)

type Rename struct {
	name         string
	label        string
	newName      []byte
	base36       bool
	keepBarcode  bool
	mergeBarcode bool
	allReads     bool
}

func NewRename(data []byte) (*Rename, error) {
	r := Rename{name: "rename"}
	// label
	label, err := jsonparser.GetUnsafeString(data, "label")
	if err != nil && err != jsonparser.KeyPathNotFoundError {
		return &r, err
	}
	if label == "" {
		r.label = r.name
	} else {
		r.label = label
	}
	newName, err := jsonparser.GetUnsafeString(data, "new_name")
	if err != nil {
		if errors.Is(err, jsonparser.KeyPathNotFoundError) {
			return &r, fmt.Errorf("%w: %v", err, "new_name")
		}
		return &r, err
	} else {
		r.newName = []byte(newName)
	}
	base36, err := jsonparser.GetBoolean(data, "base36")
	if err == jsonparser.KeyPathNotFoundError {
		r.base36 = false
	} else if err != nil {
		return &r, err
	} else {
		r.base36 = base36
	}
	keepBarcode, err := jsonparser.GetBoolean(data, "keep_barcode")
	if err == jsonparser.KeyPathNotFoundError {
		r.keepBarcode = false
	} else if err != nil {
		return &r, err
	} else {
		r.keepBarcode = keepBarcode
	}
	mergeBarcode, err := jsonparser.GetBoolean(data, "merge_barcode")
	if err == jsonparser.KeyPathNotFoundError {
		r.mergeBarcode = false
	} else if err != nil {
		return &r, err
	} else {
		r.mergeBarcode = mergeBarcode
	}
	allReads, err := jsonparser.GetBoolean(data, "all_reads")
	if err == jsonparser.KeyPathNotFoundError {
		r.allReads = true
	} else if err != nil {
		return &r, err
	} else {
		r.allReads = allReads
	}
	return &r, nil
}

func (op *Rename) Name() string {
	return op.name
}

func (op *Rename) Label() string {
	return op.label
}

func (op *Rename) IsThreadSafe() bool {
	return false
}

func (op *Rename) GetDpx(idx int) ([][]byte, int) {
	return [][]byte{}, idx
}

func (op *Rename) Transform(p *fastq.ExtPair, r int, ot *OpStat, verboseLevel int) int {
	var barcode []byte
	if r == 1 || op.allReads {
		if op.keepBarcode {
			barcode = getBarcode(p.R1.Name)
			if op.mergeBarcode {
				barcode = mergeBarcode(barcode)
			}
		}
		if op.base36 {
			p.R1.Name = joinThree(op.newName, []byte(strconv.FormatUint(p.ID, 36)), barcode)
		} else {
			p.R1.Name = joinThree(op.newName, []byte(strconv.FormatUint(p.ID, 10)), barcode)
		}
	}
	if r == 2 || op.allReads {
		if op.keepBarcode {
			barcode = getBarcode(p.R2.Name)
			if op.mergeBarcode {
				barcode = mergeBarcode(barcode)
			}
		}
		if op.base36 {
			p.R2.Name = joinThree(op.newName, []byte(strconv.FormatUint(p.ID, 36)), barcode)
		} else {
			p.R2.Name = joinThree(op.newName, []byte(strconv.FormatUint(p.ID, 10)), barcode)
		}
	}
	return 0
}

func getBarcode(name []byte) []byte {
	cutStart := len(name) - 1
	for i := cutStart; i >= 0; i-- {
		if !bio.IsDNA(name[i]) && name[i] != '#' {
			break
		}
		if name[i] == '#' {
			cutStart = i
		}
	}
	if len(name)-cutStart > 1 {
		return name[cutStart:]
	} else {
		return []byte{}
	}
}

func mergeBarcode(barcode []byte) []byte {
	var newBarcode []byte
	// Keep barcode delimiter
	newBarcode = append(newBarcode, barcode[0])
	// Add all non-delimiter letters
	for i := 1; i < len(barcode); i++ {
		if barcode[i] != '#' {
			newBarcode = append(newBarcode, barcode[i])
		}
	}
	return newBarcode
}
