//
// Copyright © 2025 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//

package trim

import (
	"fmt"

	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"
)

// Trim by quality
func TrimQuality(r *fastq.Record, window int, unqualifiedPropMax float32, minQuality int, asciiMin int, trimSide int, applyTrimSeq bool, verboseLevel int) (TrimType, int, float32, []byte) {
	var nLow, trimIdx, trimStart, trimEnd int
	var trimType TrimType
	var trimScore, tmpScore float32
	var trimSeq []byte
	trimEnd = len(r.Seq)
	if len(r.Seq) >= window {
		if trimSide == 5 {
			for i := 0; i <= len(r.Seq)-window; i++ {
				nLow = 0
				for j := i; j < i+window; j++ {
					if int(r.Qual[j])-asciiMin < minQuality {
						nLow++
					}
				}
				tmpScore = float32(nLow) / float32(window)
				if verboseLevel > 3 {
					qual := make([]int, window)
					for k := range qual {
						qual[k] = int(r.Qual[i:i+window][k]) - asciiMin
					}
					fmt.Printf("> %v [%s] [%s]\n> #-unqualified:%d prop:%.2f", qual, r.Qual[i:i+window], r.Seq[i:i+window], nLow, tmpScore)
				}
				if tmpScore >= unqualifiedPropMax {
					trimType = TrimExactType
					trimStart = i + window
					trimScore = tmpScore
					if verboseLevel > 3 {
						fmt.Printf(" trimType:%s trimStart:%d\n", trimType, trimStart)
					}
				} else {
					if verboseLevel > 3 {
						fmt.Printf(" trimType:%s trimStart:%d stop-sliding\n", trimType, trimStart)
					}
					break
				}
			}
			if trimType != NoTrimType {
				trimSeq = r.Seq[:trimStart]
			}
		} else if trimSide == 3 {
			for i := len(r.Seq) - 1; i >= window-1; i-- {
				nLow = 0
				for j := i; j > i-window; j-- {
					if int(r.Qual[j])-asciiMin < minQuality {
						nLow++
					}
				}
				tmpScore = float32(nLow) / float32(window)
				if verboseLevel > 3 {
					qual := make([]int, window)
					for k := range qual {
						qual[k] = int(r.Qual[i-window+1:i+1][k]) - asciiMin
					}
					fmt.Printf("> %v [%s] [%s]\n> #-unqualified:%d prop:%.2f", qual, r.Qual[i-window+1:i+1], r.Seq[i-window+1:i+1], nLow, tmpScore)
				}
				if tmpScore >= unqualifiedPropMax {
					trimType = TrimExactType
					trimEnd = i - window + 1
					trimScore = tmpScore
					if verboseLevel > 3 {
						fmt.Printf(" trimType:%s trimEnd:%d\n", trimType, trimEnd)
					}
				} else {
					if verboseLevel > 3 {
						fmt.Printf(" trimType:%s trimEnd:%d stop-sliding\n", trimType, trimEnd)
					}
					break
				}
			}
			if trimType != NoTrimType {
				trimSeq = r.Seq[trimEnd:]
			}
		}
		// Apply trimming
		if applyTrimSeq && trimType != NoTrimType {
			r.Seq = r.Seq[trimStart:trimEnd]
			r.Qual = r.Qual[trimStart:trimEnd]
		}
	}
	return trimType, trimIdx, trimScore, trimSeq
}
