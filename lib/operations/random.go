//
// Copyright © 2017 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package operations

import (
	"math/rand"

	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"

	"github.com/buger/jsonparser"
)

type Random struct {
	name        string
	probability float32
}

func NewRandom(data []byte) (*Random, error) {
	r := Random{name: "random"}
	probability, err := jsonparser.GetFloat(data, "probability")
	if err == jsonparser.KeyPathNotFoundError {
		r.probability = 1.
	} else if err != nil {
		return &r, err
	} else {
		r.probability = float32(probability)
	}
	return &r, nil
}

func (op *Random) Name() string {
	return op.name
}

func (op *Random) IsThreadSafe() bool {
	return true
}

func (op *Random) GetDpx(idx int) ([][]byte, int) {
	return [][]byte{}, idx
}

func (op *Random) Transform(p *fastq.ExtPair, r int, ot *OpStat, verboseLevel int) int {
	if rand.Float32() > op.probability {
		return 1
	}
	return 0
}
