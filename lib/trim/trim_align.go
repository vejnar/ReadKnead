//
// Copyright (C) 2017-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//

package trim

import (
	"bytes"
	"fmt"
	"strings"

	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"

	"github.com/vejnar/nwalgo"
)

// Trim by aligning
func TrimAlign(r *fastq.Record, adaptors [][]byte, minAdaptor int, minScore float32, trimSide int, applyTrimSeq bool, verboseLevel int) (TrimType, int, float32, []byte) {
	var trimIdx, trimStart, trimEnd, alignScore int
	var trimType TrimType
	var trimScore, tmpScore float32
	var trimSeq []byte
	trimEnd = len(r.Seq)
	// Find matching adaptor
	for iadaptor, adaptor := range adaptors {
		// Find the adaptor by search first
		adaptorIndex := bytes.Index(r.Seq, adaptor)
		if adaptorIndex == -1 {
			var alignSeq, alignAdaptor []byte
			var adaptorIndexAln int
			var minAdaptorOk bool
			// Align
			alignSeq, alignAdaptor, alignScore = nwalgo.AlignBytes(r.Seq, adaptor, 5, -10, -10, true)
			tmpScore = float32(alignScore)
			if verboseLevel > 3 {
				fmt.Printf("> %s\n> %s %.0f\n", alignSeq, alignAdaptor, tmpScore)
			}
			// Check alignment quality
			if minAdaptor == 0 {
				if trimSide == 5 {
					adaptorIndexAln = lastIndexHomo(alignAdaptor, byte('-'))
				} else if trimSide == 3 {
					adaptorIndexAln = indexHomo(alignAdaptor, byte('-'))
				}
				minAdaptorOk = true
			} else {
				var requiredAdaptorPart []byte
				if trimSide == 5 {
					requiredAdaptorPart = adaptor[len(adaptor)-minAdaptor:]
					adaptorIndexAln = bytes.LastIndex(alignAdaptor, requiredAdaptorPart)
				} else if trimSide == 3 {
					requiredAdaptorPart = adaptor[:minAdaptor]
					adaptorIndexAln = bytes.Index(alignAdaptor, requiredAdaptorPart)
				}
				minAdaptorOk = bytes.Equal(alignSeq[adaptorIndexAln:adaptorIndexAln+minAdaptor], requiredAdaptorPart)
			}
			if adaptorIndexAln != -1 && minAdaptorOk && minScore <= tmpScore && tmpScore > trimScore {
				trimType = TrimAlignType
				trimIdx = iadaptor
				trimScore = tmpScore
				if trimSide == 5 {
					trimStart = adaptorIndexAln + minAdaptor - bytes.Count(alignSeq[:adaptorIndexAln], []byte("-"))
					trimSeq = r.Seq[:trimStart]
				} else if trimSide == 3 {
					trimEnd = adaptorIndexAln - bytes.Count(alignSeq[:adaptorIndexAln], []byte("-"))
					trimSeq = r.Seq[trimEnd:]
				}
			}
		} else {
			trimType = TrimExactType
			trimIdx = iadaptor
			// Set score to 10x per match
			trimScore = 10. * float32(len(adaptor))
			if trimSide == 5 {
				trimStart = adaptorIndex + len(adaptor)
				trimSeq = r.Seq[:trimStart]
			} else if trimSide == 3 {
				trimEnd = adaptorIndex
				trimSeq = r.Seq[trimEnd:]
			}
			if verboseLevel > 3 {
				fmt.Printf("> %s\n> %s%s %.0f\n", r.Seq, strings.Repeat("-", adaptorIndex), adaptor, trimScore)
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

func indexHomo(s []byte, sep byte) int {
	if s[0] != sep {
		return -1
	}
	for i := 0; i < len(s)-1; i++ {
		if s[i] != sep {
			return i
		}
	}
	return -1
}

func lastIndexHomo(s []byte, sep byte) int {
	last := len(s) - 1
	if s[last] != sep {
		return -1
	}
	for i := last; i >= 0; i-- {
		if s[i] != sep {
			return i
		}
	}
	return -1
}
