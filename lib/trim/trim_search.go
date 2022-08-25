//
// Copyright (C) 2017-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//

package trim

import (
	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"
)

// Trim by search
func TrimSearch(r *fastq.Record, adaptors [][]byte, minAdaptor int, minScore float32, trimSide int, applyTrimSeq bool, verboseLevel int) (TrimType, int, float32, []byte) {
	var trimIdx, trimStart, trimEnd, nMatch, nMismatch int
	var trimType TrimType
	var trimScore, tmpScore float32
	var trimSeq []byte
	trimEnd = len(r.Seq)
	// Find matching adaptor
	for iadaptor, adaptor := range adaptors {
		if trimSide == 5 {
			for i := len(r.Seq); i >= minAdaptor; i-- {
				nMatch = 0
				nMismatch = 0
				subseq := r.Seq[:i]
				for j, k := len(subseq)-1, len(adaptor)-1; j >= 0 && k >= 0; j, k = j-1, k-1 {
					if subseq[j] == adaptor[k] {
						nMatch++
					} else {
						nMismatch++
					}
				}
				tmpScore = float32(nMatch) / float32(nMatch+nMismatch)
				if tmpScore > trimScore {
					if nMismatch == 0 {
						trimType = TrimExactType
						trimScore = 1.
					} else if tmpScore >= minScore {
						trimType = TrimAlignType
						trimScore = tmpScore
					}
					if trimType != NoTrimType {
						trimIdx = iadaptor
						trimStart = i
						trimSeq = r.Seq[:trimStart]
						break
					}
				}
			}
		} else if trimSide == 3 {
			for i := 0; i <= len(r.Seq)-minAdaptor; i++ {
				nMatch = 0
				nMismatch = 0
				subseq := r.Seq[i:]
				for j := 0; j < len(subseq) && j < len(adaptor); j++ {
					if subseq[j] == adaptor[j] {
						nMatch++
					} else {
						nMismatch++
					}
				}
				tmpScore = float32(nMatch) / float32(nMatch+nMismatch)
				if tmpScore > trimScore {
					if nMismatch == 0 {
						trimType = TrimExactType
						trimScore = 1.
					} else if tmpScore >= minScore {
						trimType = TrimAlignType
						trimScore = tmpScore
					}
					if trimType != NoTrimType {
						trimIdx = iadaptor
						trimEnd = i
						trimSeq = r.Seq[trimEnd:]
						break
					}
				}
			}
		}
	}
	// Apply trimming
	if applyTrimSeq && trimType != NoTrimType {
		r.Seq = r.Seq[trimStart:trimEnd]
		r.Qual = r.Qual[trimStart:trimEnd]
	}
	return trimType, trimIdx, trimScore, trimSeq
}
