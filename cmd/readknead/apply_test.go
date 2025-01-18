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

type TestSet struct {
	fastqsR1     string
	fastqsR2     string
	fqFnameOutR1 string
	fqFnameOutR2 string
	opsR1Path    string
	opsR2Path    string
	goldenPath   string
}

func TestApplyOperations(t *testing.T) {
	tmp := t.TempDir()

	tests := []TestSet{
		{
			fastqsR1:     "sample1_R1.fastq",
			fastqsR2:     "",
			fqFnameOutR1: "sample1_R1.fastq",
			fqFnameOutR2: "",
			opsR1Path:    "clip_trim.json",
			goldenPath:   "sample1_R1.fastq.golden",
		},
		{
			fastqsR1:     "sample2_R1.fastq",
			fastqsR2:     "sample2_R2.fastq",
			fqFnameOutR1: "sample2_R1.fastq",
			fqFnameOutR2: "sample2_R2.fastq",
			opsR1Path:    "paired_end_trim.json",
			goldenPath:   "sample2_R*.fastq.golden",
		},
		{
			fastqsR1:     "sample2_R1.fastq",
			fastqsR2:     "sample2_R2.fastq",
			fqFnameOutR1: "sample2_demultiplex_[DPX]_R1.fastq",
			fqFnameOutR2: "sample2_demultiplex_[DPX]_R2.fastq",
			opsR1Path:    "demultiplex.json",
			goldenPath:   "sample2_demultiplex_*_R*.fastq.golden",
		},
		{
			fastqsR1:     "sample2_R1.fastq",
			fastqsR2:     "sample2_R2.fastq",
			fqFnameOutR1: "sample2_filter_quality_R1.fastq",
			fqFnameOutR2: "sample2_filter_quality_R2.fastq",
			opsR1Path:    "filter_quality.json",
			opsR2Path:    "filter_quality.json",
			goldenPath:   "sample2_filter_quality_R*.fastq.golden",
		},
		{
			fastqsR1:     "sample3_R1.fastq",
			fastqsR2:     "",
			fqFnameOutR1: "sample3_R1.fastq",
			fqFnameOutR2: "",
			opsR1Path:    "trimqual.json",
			goldenPath:   "sample3_R1.fastq.golden",
		},
		{
			fastqsR1:     "sample3_R1.fastq",
			fastqsR2:     "",
			fqFnameOutR1: "sample3_filter_quality_R1.fastq",
			fqFnameOutR2: "",
			opsR1Path:    "filter_quality.json",
			goldenPath:   "sample3_filter_quality_R1.fastq.golden",
		},
		{
			fastqsR1:     "sample4_R1.fastq",
			fastqsR2:     "",
			fqFnameOutR1: "sample4_R1.fastq",
			fqFnameOutR2: "",
			opsR1Path:    "trim_length.json",
			goldenPath:   "sample4_R1.fastq.golden",
		},
	}

	for _, test := range tests {
		var fastqsR1, fastqsR2 []string
		fastqsR1 = []string{filepath.Join("testdata", test.fastqsR1)}
		if len(test.fastqsR2) > 0 {
			fastqsR2 = []string{filepath.Join("testdata", test.fastqsR2)}
		}
		fqPathOut := tmp
		fqCmdIn := []string{}
		fqCmdOut := []string{}
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
		opsR1, err = operations.ReadOps(readAll(filepath.Join("testdata", test.opsR1Path)), param)
		if err != nil {
			t.Fatalf("failed reading json: %s", err)
		}
		if test.opsR2Path != "" {
			opsR2, err = operations.ReadOps(readAll(filepath.Join("testdata", test.opsR2Path)), param)
			if err != nil {
				t.Fatalf("failed reading json: %s", err)
			}
		}

		// Run
		nPair, err := ApplyOperations(fastqsR1, fastqsR2, fqPathOut, test.fqFnameOutR1, test.fqFnameOutR2, fqCmdIn, fqCmdOut, opsR1, opsR2, param, statsInPath, statsOutPath, maxReadLength, reportPath, label, bufSize, nWorker, verboseLevel)
		if err != nil {
			t.Fatalf("apply failed: %s", err)
		}
		if nPair != 4 {
			t.Error("apply returned", nPair, "expected 4")
		}

		// Wait for the output file(s) to be available
		time.Sleep(100 * time.Millisecond)

		goldenPaths, err := filepath.Glob(filepath.Join("testdata", test.goldenPath))
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
