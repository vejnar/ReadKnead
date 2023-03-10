//
// Copyright © 2017 Charles E. Vejnar
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

	"git.sr.ht/~vejnar/bktrim"
)

func NewMatrixAdapter(adaptors [][]byte, trimSide int, epsilon float64, epsilonIndel float64, minOverlap int, baseQual int) []*bktrim.Matrix {
	var ms []*bktrim.Matrix
	for _, adaptor := range adaptors {
		m := bktrim.NewMatrix(epsilon, epsilonIndel, minOverlap, baseQual)
		if trimSide == 5 {
			m.AddAdapter(adaptor, bktrim.TRIM_HEAD, 0)
		} else if trimSide == 3 {
			m.AddAdapter(adaptor, bktrim.TRIM_TAIL, 0)
		}
		ms = append(ms, m)
	}
	return ms
}

func NewMatrixAdapterPaired(adaptors [][]byte, adaptorsPaired [][]byte, epsilon float64, epsilonIndel float64, minOverlap int, baseQual int) []*bktrim.Matrix {
	var ms []*bktrim.Matrix
	for iadaptor, adaptor := range adaptors {
		m := bktrim.NewMatrix(epsilon, epsilonIndel, minOverlap, baseQual)
		m.AddAdapter(adaptor, bktrim.TRIM_TAIL, 0)
		m.AddAdapter(adaptorsPaired[iadaptor], bktrim.TRIM_TAIL, 1)
		ms = append(ms, m)
	}
	return ms
}

// Trim by aligning
func TrimBKTrim(r *fastq.Record, bkMatrices []*bktrim.Matrix, minAdaptor int, trimSide int, applyTrimSeq bool, verboseLevel int) (TrimType, int, float32, []byte) {
	var trimIdx, trimStart, trimEnd int
	var trimType TrimType
	var trimScore, tmpScore float32
	var trimSeq []byte
	var alignQual bool
	trimEnd = len(r.Seq)
	// Find matching adaptor
	for im, m := range bkMatrices {
		// Find the adaptor by search first
		adaptorIndex := bytes.Index(r.Seq, m.Adapter1.Seq)
		if adaptorIndex == -1 {
			determined, sol := m.FindAdapter(r.Seq, r.Qual)
			if determined {
				tmpScore = float32(sol.Score)
				if verboseLevel > 3 {
					n := 0
					if sol.Pos > 0 {
						n = sol.Pos
					}
					fmt.Printf("> %s\n> %s%d %.2f\n", r.Seq, strings.Repeat("-", n), sol.Pos, tmpScore)
				}
				if sol.Pos <= 0 {
					trimType = TrimTooShortType
					trimIdx = im
					trimScore = tmpScore
					trimSeq = r.Seq
					trimStart = 0
					trimEnd = 0
				} else {
					// Check alignment quality
					alignQual = true
					if minAdaptor > 0 {
						if trimSide == 5 {
							if sol.Pos <= minAdaptor {
								alignQual = false
							} else {
								alignQual = bytes.Equal(m.Adapter1.Seq[len(m.Adapter1.Seq)-minAdaptor:], r.Seq[sol.Pos-minAdaptor:sol.Pos])
							}
						} else if trimSide == 3 {
							if len(r.Seq) < sol.Pos+minAdaptor {
								alignQual = false
							} else {
								alignQual = bytes.Equal(m.Adapter1.Seq[:minAdaptor], r.Seq[sol.Pos:sol.Pos+minAdaptor])
							}
						}
					}
					if alignQual && trimScore < tmpScore {
						trimType = TrimAlignType
						trimIdx = im
						trimScore = tmpScore
						if trimSide == 5 {
							trimStart = sol.Pos
							trimSeq = r.Seq[:trimStart]
						} else if trimSide == 3 {
							trimEnd = sol.Pos
							trimSeq = r.Seq[trimEnd:]
						}
					}
				}
			}
		} else {
			trimType = TrimExactType
			trimIdx = im
			trimScore = 0.
			if trimSide == 5 {
				trimStart = adaptorIndex + len(m.Adapter1.Seq)
				trimSeq = r.Seq[:trimStart]
			} else if trimSide == 3 {
				trimEnd = adaptorIndex
				trimSeq = r.Seq[trimEnd:]
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

// Trim by aligning
func TrimBKTrimPaired(p *fastq.ExtPair, bkMatrices []*bktrim.Matrix, applyTrimSeq bool, verboseLevel int) (TrimType, int, float32, []byte) {
	var trimIdx int
	var trimType TrimType
	var trimScore float64
	var trimSeq []byte
	var trimMatrix *bktrim.Matrix
	var sol1, sol2 bktrim.Solution
	// Find matching adaptor
	for im, m := range bkMatrices {
		found, tmpSol1, tmpSol2 := m.FindAdapterWithPE(p.R1.Seq, p.R1.Qual, p.R2.Seq, p.R2.Qual)
		if found && trimScore <= tmpSol1.Score+tmpSol2.Score {
			trimType = TrimAlignType
			trimIdx = im
			trimScore = tmpSol1.Score + tmpSol2.Score
			trimMatrix = m
			sol1 = tmpSol1
			sol2 = tmpSol2
		}
	}
	// Apply trimming
	if applyTrimSeq && trimType != NoTrimType {
		if trimMatrix.CombinePairSeqs(p.R1.Seq, p.R1.Qual, p.R2.Seq, p.R2.Qual, sol1.Pos, sol2.Pos) {
			p.R1.Seq = p.R1.Seq[:sol1.Pos]
			p.R1.Qual = p.R1.Qual[:sol1.Pos]
			p.R2.Seq = p.R2.Seq[:sol2.Pos]
			p.R2.Qual = p.R2.Qual[:sol2.Pos]
			return trimType, trimIdx, float32(trimScore), p.R1.Seq[sol1.Pos:]
		}
	}
	return trimType, trimIdx, float32(trimScore), trimSeq
}
