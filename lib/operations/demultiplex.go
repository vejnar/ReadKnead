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

	"git.sr.ht/~vejnar/ReadKnead/lib/bio"
	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"

	"github.com/buger/jsonparser"
)

type Demultiplex struct {
	Barcodes     [][]byte
	BarcodesID   []int
	name         string
	end          int
	barcodeIdx   int
	lengthLigand int
	maxMismatch  int
	useSeq       bool
}

func NewDemultiplex(data []byte) (*Demultiplex, error) {
	d := Demultiplex{name: "demultiplex"}
	var err error
	jsonparser.ArrayEach(data, func(value []byte, dataType jsonparser.ValueType, offset int, err2 error) {
		if err == nil {
			var b string
			b, err = jsonparser.ParseString(value)
			if err != nil {
				return
			}
			d.Barcodes = append(d.Barcodes, []byte(b))
		}
	}, "barcodes")
	if err != nil {
		return &d, err
	}
	end, err := jsonparser.GetInt(data, "end")
	if err == jsonparser.KeyPathNotFoundError {
		d.useSeq = false
	} else if err != nil {
		return &d, err
	} else {
		d.end = int(end)
		d.useSeq = true
	}
	bcidx, err := jsonparser.GetInt(data, "barcode_idx")
	if err == jsonparser.KeyPathNotFoundError {
		d.useSeq = true
	} else if err != nil {
		return &d, err
	} else {
		d.barcodeIdx = int(bcidx)
		d.useSeq = false
	}
	maxMismatch, err := jsonparser.GetInt(data, "max_mismatch")
	if err == jsonparser.KeyPathNotFoundError {
		d.maxMismatch = 0
	} else if err != nil {
		return &d, err
	} else {
		d.maxMismatch = int(maxMismatch)
	}
	lengthLigand, err := jsonparser.GetInt(data, "length_ligand")
	if err == jsonparser.KeyPathNotFoundError {
		d.lengthLigand = 0
	} else if err != nil {
		return &d, err
	} else {
		d.lengthLigand = int(lengthLigand)
	}
	return &d, nil
}

func (op *Demultiplex) Name() string {
	return op.name
}

func (op *Demultiplex) IsThreadSafe() bool {
	return true
}

func (op *Demultiplex) GetDpx(idx int) ([][]byte, int) {
	var names [][]byte
	names = append(names, []byte("undetermined"))
	idx++
	for _, b := range op.Barcodes {
		names = append(names, b)
		op.BarcodesID = append(op.BarcodesID, idx)
		idx++
	}
	return names, idx
}

func (op *Demultiplex) Transform(p *fastq.ExtPair, r int, ot *OpStat, verboseLevel int) int {
	var seq []byte
	var bcs [][]byte
	var okSeq bool
	bestBarcodeSeq := []byte("undetermined")
	bestBarcode := -1
	if r == 1 {
		if verboseLevel > 2 {
			fmt.Printf("%s %s r%d\n%s\n", op.name, p.R1.Name, r, p.R1.Seq)
		}
		for ibc, bc := range op.Barcodes {
			okSeq = false
			if op.useSeq {
				if len(p.R1.Seq) > len(bc) {
					if op.end == 5 {
						seq = p.R1.Seq[:len(bc)]
					} else {
						seq = p.R1.Seq[len(p.R1.Seq)-len(bc):]
					}
					okSeq = true
				}
			} else {
				bcs = GetBarcodes(p.R1.Name)
				seq = bcs[op.barcodeIdx]
				if len(seq) == len(bc) {
					okSeq = true
				}
			}
			if okSeq {
				// Count mismatch(es) with barcode
				nmismatch := 0
				for i := 0; i < len(bc); i++ {
					if bc[i] != seq[i] {
						nmismatch++
					}
				}
				// Demultiplex
				if nmismatch <= op.maxMismatch {
					bestBarcode = op.BarcodesID[ibc]
					bestBarcodeSeq = bc
				}
				if verboseLevel > 3 {
					fmt.Printf("+ %s barcode:%s nmismatch:%d\n", bc, seq, nmismatch)
				}
				if nmismatch == 0 {
					break
				}
			}
		}
	} else {
		if verboseLevel > 2 {
			fmt.Printf("%s %s r%d\n%s\n", op.name, p.R2.Name, r, p.R2.Seq)
		}
		for ibc, bc := range op.Barcodes {
			okSeq = false
			if op.useSeq {
				if len(p.R2.Seq) > len(bc) {
					if op.end == 5 {
						seq = p.R2.Seq[:len(bc)]
					} else {
						seq = p.R2.Seq[len(p.R2.Seq)-len(bc):]
					}
					okSeq = true
				}
			} else {
				bcs = GetBarcodes(p.R2.Name)
				seq = bcs[op.barcodeIdx]
				if len(seq) == len(bc) {
					okSeq = true
				}
			}
			if okSeq {
				// Count mismatch(es) with barcode
				nmismatch := 0
				for i := 0; i < len(bc); i++ {
					if bc[i] != seq[i] {
						nmismatch++
					}
				}
				// Demultiplex
				if nmismatch <= op.maxMismatch {
					bestBarcode = op.BarcodesID[ibc]
					bestBarcodeSeq = bc
				}
				if verboseLevel > 3 {
					fmt.Printf("+ %s barcode:%s nmismatch:%d\n", bc, seq, nmismatch)
				}
				if nmismatch == 0 {
					break
				}
			}
		}
	}
	// Demultiplex if barcode found
	if bestBarcode != -1 {
		p.WID = bestBarcode
		// Clip barcode and ligand
		pc := Clip{name: op.name + "-clip", end: op.end, length: len(bestBarcodeSeq) + op.lengthLigand}
		pc.Transform(p, r, ot, verboseLevel)
		if verboseLevel > 2 {
			fmt.Printf("barcode found:%s\n", string(bestBarcodeSeq))
		}
	} else {
		if verboseLevel > 2 {
			fmt.Println("No barcode found")
		}
	}
	// Stats
	if r == 1 {
		ot.OpsR1[op.name][string(bestBarcodeSeq)]++
	} else {
		ot.OpsR2[op.name][string(bestBarcodeSeq)]++
	}
	return 0
}

func GetBarcodes(name []byte) [][]byte {
	var barcodes [][]byte
	var bc []byte
	lastCut := len(name)
	for i := len(name) - 1; i >= 0; i-- {
		if !bio.IsDNA(name[i]) && name[i] != '#' {
			break
		}
		if name[i] == '#' {
			bc = name[i+1 : lastCut]
			if len(bc) > 0 {
				barcodes = append([][]byte{bc}, barcodes...)
			}
			lastCut = i
		}
	}
	return barcodes
}
