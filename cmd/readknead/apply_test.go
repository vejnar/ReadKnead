//
// Copyright © 2017 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package main

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
	"time"

	"git.sr.ht/~vejnar/ReadKnead/lib/operations"
	"git.sr.ht/~vejnar/ReadKnead/lib/param"
)

func TestApplyOperations(t *testing.T) {
	tmp := t.TempDir()

	tests := [][]string{
		{
			"sample1_R1.fastq",
			"",
			"sample1_R1.fastq",
			"",
			"clip_trim.json",
			"sample1_R1.fastq.golden",
		},
		{
			"sample2_R1.fastq",
			"sample2_R2.fastq",
			"sample2_R1.fastq",
			"sample2_R2.fastq",
			"paired_end_trim.json",
			"sample2_R*.fastq.golden",
		},
		{
			"sample2_R1.fastq",
			"sample2_R2.fastq",
			"sample2_[DPX]_R1.fastq",
			"sample2_[DPX]_R2.fastq",
			"demultiplex.json",
			"sample2_*_R*.fastq.golden",
		},
		{
			"sample3_R1.fastq",
			"",
			"sample3_R1.fastq",
			"",
			"trimqual.json",
			"sample3_R1.fastq.golden",
		},
		{
			"sample4_R1.fastq",
			"",
			"sample4_R1.fastq",
			"",
			"trim_length.json",
			"sample4_R1.fastq.golden",
		},
	}

	for _, test := range tests {
		var fastqsR1, fastqsR2 []string
		fastqsR1 = []string{filepath.Join("testdata", test[0])}
		if len(test[1]) > 0 {
			fastqsR2 = []string{filepath.Join("testdata", test[1])}
		}
		fqPathOut := tmp
		fqFnameOutR1 := test[2]
		fqFnameOutR2 := test[3]
		fqCmdIn := []string{}
		fqCmdOut := []string{}
		opsR1Path := test[4]
		statsInPath := ""
		statsOutPath := ""
		maxReadLength := 1000
		reportPath := ""
		label := ""
		bufSize := 41943040
		nWorker := 1
		verboseLevel := 10

		param := param.Parameters{AsciiMin: 33, MaxQual: 43, Paired: false}
		if len(fastqsR2) > 0 {
			param.Paired = true
		}

		// Load read operations from JSON
		var opsR1, opsR2 []operations.Operation
		var err error
		opsR1, err = operations.ReadOps(readAll(filepath.Join("testdata", opsR1Path)), param)
		if err != nil {
			t.Fatalf("failed reading json: %s", err)
		}

		// Run
		nPair, err := ApplyOperations(fastqsR1, fastqsR2, fqPathOut, fqFnameOutR1, fqFnameOutR2, fqCmdIn, fqCmdOut, opsR1, opsR2, param, statsInPath, statsOutPath, maxReadLength, reportPath, label, bufSize, nWorker, verboseLevel)
		if err != nil {
			t.Fatalf("apply failed: %s", err)
		}
		if nPair != 4 {
			t.Error("apply returned", nPair, "expected 4")
		}

		// Wait for the output file(s) to be available
		time.Sleep(100 * time.Millisecond)

		goldenPaths, err := filepath.Glob(filepath.Join("testdata", test[5]))
		if err != nil {
			t.Fatal(err)
		}

		for _, goldenPath := range goldenPaths {
			g, err := os.ReadFile(goldenPath)
			if err != nil {
				t.Fatalf("failed reading .golden: %s", err)
			}

			_, filename := filepath.Split(goldenPath)
			o, err := os.ReadFile(filepath.Join(tmp, filename[:len(filename)-len(".golden")]))
			if err != nil {
				t.Fatalf("failed reading output: %s", err)
			}

			if !bytes.Equal(g, o) {
				t.Errorf("written json does not match .golden file X")
			}
		}
	}
}
