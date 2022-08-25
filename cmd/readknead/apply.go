//
// Copyright (C) 2017-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package main

import (
	"context"
	"fmt"
	"path/filepath"
	"strings"

	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"
	"git.sr.ht/~vejnar/ReadKnead/lib/operations"
	"git.sr.ht/~vejnar/ReadKnead/lib/param"

	"golang.org/x/sync/errgroup"
)

func ApplyOperations(fastqsR1 []string, fastqsR2 []string, fqPathOut string, fqFnameOutR1 string, fqFnameOutR2 string, fqCmdIn []string, fqCmdOut []string, opsR1 []operations.Operation, opsR2 []operations.Operation, param param.Parameters, statsInPath string, statsOutPath string, maxReadLength int, reportPath string, label string, bufSize int, nWorker int, verboseLevel int) (nPair uint64, err error) {
	// Check there is a single non thread-safe operation
	var nNotThreadSafe int
	for _, op := range opsR1 {
		if !op.IsThreadSafe() {
			nNotThreadSafe++
		}
	}
	for _, op := range opsR2 {
		if !op.IsThreadSafe() {
			nNotThreadSafe++
		}
	}
	if nNotThreadSafe > 1 {
		return nPair, fmt.Errorf("Too many non thread-safe operation (maximum 1)")
	}

	// Demultiplex name(s)
	var dpxNames, names [][]byte
	var dpxID int
	for _, op := range opsR1 {
		names, dpxID = op.GetDpx(dpxID)
		dpxNames = append(dpxNames, names...)
	}
	for _, op := range opsR2 {
		names, dpxID = op.GetDpx(dpxID)
		dpxNames = append(dpxNames, names...)
	}
	if len(dpxNames) == 0 {
		dpxNames = append(dpxNames, []byte("all"))
	}
	if verboseLevel > 2 {
		fmt.Printf("Barcodes: ")
		for _, n := range dpxNames {
			fmt.Printf("%s ", string(n))
		}
		fmt.Printf("\n")
	}

	// Open output FASTQ files
	var fqws1, fqws2 []*fastq.FqWriter
	var fqw *fastq.FqWriter
	var writeFq bool
	if fqPathOut != "" || fqFnameOutR1 != "" || fqFnameOutR2 != "" {
		for _, n := range dpxNames {
			fqf := filepath.Base(fastqsR1[0])
			if fqFnameOutR1 != "" {
				fqf = strings.Replace(fqFnameOutR1, "[DPX]", string(n), 1)
			}
			if verboseLevel > 2 {
				fmt.Println("Opening", filepath.Join(fqPathOut, fqf))
			}
			fqw, err = fastq.Wopen(filepath.Join(fqPathOut, fqf), fqCmdOut, bufSize)
			if err != nil {
				return nPair, err
			}
			fqws1 = append(fqws1, fqw)
			defer func(fqw *fastq.FqWriter) {
				if ferr := fqw.Close(); ferr != nil {
					if err != nil {
						err = fmt.Errorf("%w Then %s", err, ferr)
					} else {
						err = ferr
					}
				}
			}(fqw)
			if param.Paired {
				fqf = filepath.Base(fastqsR2[0])
				if fqFnameOutR2 != "" {
					fqf = strings.Replace(fqFnameOutR2, "[DPX]", string(n), 1)
				}
				if verboseLevel > 2 {
					fmt.Println("Opening", filepath.Join(fqPathOut, fqf))
				}
				fqw, err = fastq.Wopen(filepath.Join(fqPathOut, fqf), fqCmdOut, bufSize)
				if err != nil {
					return nPair, err
				}
				fqws2 = append(fqws2, fqw)
				defer func(fqw *fastq.FqWriter) {
					if ferr := fqw.Close(); ferr != nil {
						if err != nil {
							err = fmt.Errorf("%w Then %s", err, ferr)
						} else {
							err = ferr
						}
					}
				}(fqw)
			}
		}
		writeFq = true
	}

	// Init context
	ctx, _ := context.WithCancel(context.Background()) // cancelFunc not used
	// Start sync errgroup
	g, gctx := errgroup.WithContext(ctx)

	// Start receiving channel
	chFinal := make(chan fastq.ExtPair, nWorker*10000)
	chTransit := make(chan fastq.ExtPair, nWorker)

	g.Go(func() error {
		defer close(chFinal)
		var id, maxid uint64
		pairCache := make(map[uint64]fastq.ExtPair)
		for p := range chTransit {
			if p.ID == id {
				chFinal <- p
				id++
			} else {
				pairCache[p.ID] = p
				for {
					if pc, ok := pairCache[id]; ok {
						chFinal <- pc
						delete(pairCache, pc.ID)
						id++
					} else {
						break
					}
				}
			}
		}
		// Send back what remains in the cache
		if len(pairCache) > 0 {
			for k := range pairCache {
				if k > maxid {
					maxid = k
				}
			}
			for i := id; i <= maxid; i++ {
				pc := pairCache[i]
				chFinal <- pc
			}
		}
		return nil
	})

	// Start read channel
	chPair := make(chan fastq.ExtPair, nWorker*2)

	g.Go(func() error {
		defer close(chPair)
		var id uint64
		var fqr1, fqr2 *fastq.FqReader
		var r1, r2 fastq.Record
		var err error
		for iFq := 0; iFq < len(fastqsR1); iFq++ {
			// Open FASTQ files
			fqr1, err = fastq.Ropen(fastqsR1[iFq], fqCmdIn, bufSize)
			if err != nil {
				return err
			}
			defer fqr1.Close()
			if param.Paired {
				fqr2, err = fastq.Ropen(fastqsR2[iFq], fqCmdIn, bufSize)
				if err != nil {
					return err
				}
				defer fqr2.Close()
				r2, err = fqr2.Iter()
				if err != nil {
					return err
				}
			} else {
				fqr2 = new(fastq.FqReader)
			}
			// Iter reads
			for r1, err = fqr1.Iter(); !fqr1.Done && !fqr2.Done; r1, err = fqr1.Iter() {
				if err != nil {
					return err
				}
				select {
				case <-gctx.Done():
					return gctx.Err()
				case chPair <- fastq.ExtPair{ID: id, Ok: true, R1: r1, R2: r2}:
				}
				id++
				if param.Paired {
					r2, err = fqr2.Iter()
					if err != nil {
						return err
					}
				}
			}
		}
		return nil
	})

	// Init. OpStat
	ots := make([]*operations.OpStat, nWorker)
	for i := 0; i < nWorker; i++ {
		ots[i] = operations.NewOpStat(statsInPath, statsOutPath, reportPath, label, maxReadLength, param.MaxQual, param.AsciiMin, param.Paired, opsR1, opsR2)
	}

	// Spawn worker goroutine(s)
	g.Go(func() error {
		defer close(chTransit)
		// Get thread-safe operation(s)
		var opsR1ts []operations.Operation
		var opsR2ts []operations.Operation
		for iop := 0; iop < len(opsR1); iop++ {
			if opsR1[iop].IsThreadSafe() {
				opsR1ts = append(opsR1ts, opsR1[iop])
			}
		}
		for iop := 0; iop < len(opsR2); iop++ {
			if opsR2[iop].IsThreadSafe() {
				opsR2ts = append(opsR2ts, opsR2[iop])
			}
		}
		// Start worker(s)
		wg, wgctx := errgroup.WithContext(gctx)
		chWorker := make(chan int, nWorker)
		for i := 0; i < nWorker; i++ {
			chWorker <- i
			wg.Go(func() error {
				nw := <-chWorker
				// Loop over data
				var iop, ipair int
				var ok1, ok2 bool
				for p := range chPair {
					// Stat In
					ots[nw].CountIn(&p)
					// Apply operation(s)
					ok1 = true
					ok2 = true
					if verboseLevel > 2 {
						fmt.Println("\n******** Read/pair", ipair, "********")
					}
					for iop = 0; iop < len(opsR1ts); iop++ {
						if opsR1ts[iop].Transform(&p, 1, ots[nw], verboseLevel) != 0 {
							ok1 = false
							break
						}
						if verboseLevel > 2 {
							fmt.Println()
						}
					}
					for iop = 0; iop < len(opsR2ts); iop++ {
						if opsR2ts[iop].Transform(&p, 2, ots[nw], verboseLevel) != 0 {
							ok2 = false
							break
						}
						if verboseLevel > 2 {
							fmt.Println()
						}
					}
					// Send back
					if ok1 && ok2 {
						// Stat Out
						ots[nw].CountOut(&p)
						ots[nw].KeptPair++
					} else {
						p.Ok = false
					}
					ots[nw].TotalPair++
					select {
					case <-wgctx.Done():
						return wgctx.Err()
					case chTransit <- p:
					}
					ipair++
				}
				return nil
			})
		}
		// Wait for the workers to finish
		err := wg.Wait()
		if err != nil {
			return err
		}
		return nil
	})

	// Get the non thread-safe operation (one per read)
	var opNt1, opNt2 operations.Operation
	var opNtStat1, opNtStat2 *operations.OpStat
	var doNtr1, doNtr2 int
	for iop := 0; iop < len(opsR1); iop++ {
		if !opsR1[iop].IsThreadSafe() {
			opNt1 = opsR1[iop]
			doNtr1 = 1
			break
		}
	}
	for iop := 0; iop < len(opsR2); iop++ {
		if !opsR2[iop].IsThreadSafe() {
			opNt2 = opsR2[iop]
			doNtr2 = 2
			break
		}
	}

	// Write FASTQ
	var id uint64
	id = 1
	for p := range chFinal {
		// Output
		if writeFq && p.Ok {
			// Non thread-safe operation
			if doNtr1 != 0 || doNtr2 != 0 {
				p.ID = id
				id++
				if doNtr1 != 0 {
					opNt1.Transform(&p, doNtr1, opNtStat1, verboseLevel)
				}
				if doNtr2 != 0 {
					opNt2.Transform(&p, doNtr2, opNtStat2, verboseLevel)
				}
			}
			// Write
			err = fqws1[p.WID].WriteRecord(p.R1)
			if err != nil {
				return nPair, err
			}
			if param.Paired {
				err = fqws2[p.WID].WriteRecord(p.R2)
				if err != nil {
					return nPair, err
				}
			}
		}
	}

	// Wait for all errgroup goroutines
	err = g.Wait()
	if err != nil {
		return nPair, err
	}

	// Write Report
	for i := 1; i < nWorker; i++ {
		ots[0].Update(ots[i])
	}
	err = ots[0].Write()
	if err != nil {
		return nPair, err
	}

	return ots[0].TotalPair, err
}
