//
// Copyright © 2017 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//

package trim

import (
	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"
)

// Trim by matching
func TrimMatch(r *fastq.Record, adaptors [][]byte, position int, minAdaptor int, minScore float32, trimSide int, applyTrimSeq bool, verboseLevel int) (TrimType, int, float32, []byte) {
	var trimIdx, trimStart, trimEnd, adLength, nMatch int
	var trimType TrimType
	var trimScore, tmpScore float32
	var trimSeq []byte
	trimEnd = len(r.Seq)
	// Find matching adaptor
	for iad, ad := range adaptors {
		adLength = len(ad)
		// Count mismatch(es) with adaptor
		nMatch = 0
		for i := 0; i < min(adLength, len(r.Seq)-position); i++ {
			if ad[i] == r.Seq[position+i] {
				nMatch++
			}
		}
		// Trim or not trim
		tmpScore = float32(nMatch) / float32(adLength)
		if nMatch == adLength {
			trimIdx = iad
			trimType = TrimExactType
			trimScore = 1.
			break
		} else if tmpScore >= minScore && tmpScore > trimScore {
			trimIdx = iad
			trimType = TrimAlignType
			trimScore = tmpScore
		}
	}
	if trimType != NoTrimType {
		if trimSide == 5 {
			trimStart = position + len(adaptors[trimIdx])
			trimSeq = r.Seq[:trimStart]
		} else if trimSide == 3 {
			trimEnd = position
			trimSeq = r.Seq[trimEnd:]
		}
		// Apply trimming
		if applyTrimSeq {
			r.Seq = r.Seq[trimStart:trimEnd]
			r.Qual = r.Qual[trimStart:trimEnd]
		}
	}
	return trimType, trimIdx, trimScore, trimSeq
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
