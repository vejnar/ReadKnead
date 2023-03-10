//
// Copyright © 2017 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package main

import (
	"flag"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"strings"
	"time"

	"git.sr.ht/~vejnar/ReadKnead/lib/operations"
	"git.sr.ht/~vejnar/ReadKnead/lib/param"
)

var version = "DEV"

func readAll(path string) []byte {
	fos, err := os.Open(path)
	if err != nil {
		panic(err)
	}
	txt, err := ioutil.ReadAll(fos)
	if err != nil {
		panic(err)
	}
	fos.Close()
	return txt
}

func main() {
	// Arguments: General
	var verbose, printVersion bool
	var reportPath, label string
	var bufSize, nWorker, verboseLevel int
	flag.StringVar(&reportPath, "report_path", "", "Write report to path (stdout with -)")
	flag.StringVar(&label, "label", "", "Label")
	flag.IntVar(&bufSize, "buf_size", 41943040, "Buffer IO size")
	flag.IntVar(&nWorker, "num_worker", 1, "Number of worker(s)")
	flag.IntVar(&verboseLevel, "verbose_level", 0, "Verbose level")
	flag.BoolVar(&verbose, "verbose", false, "Verbose")
	flag.BoolVar(&printVersion, "version", false, "Print version and quit")
	// Arguments: FastQ
	var fqFnamesR1, fqFnamesR2, fqPathOut, fqFnameOutR1, fqFnameOutR2, fqCmdInRaw, fqCmdOutRaw string
	flag.StringVar(&fqFnamesR1, "fq_fnames_r1", "", "Path to read 1 FASTQ files (comma separated)")
	flag.StringVar(&fqFnamesR2, "fq_fnames_r2", "", "Path to read 2 FASTQ files (comma separated)")
	flag.StringVar(&fqPathOut, "fq_path_out", "", "Path to output FASTQ files")
	flag.StringVar(&fqFnameOutR1, "fq_fname_out_r1", "", "Output read 1 FASTQ file")
	flag.StringVar(&fqFnameOutR2, "fq_fname_out_r2", "", "Output read 2 FASTQ file")
	flag.StringVar(&fqCmdInRaw, "fq_command_in", "", "Command line to execute for opening each input file (comma separated)")
	flag.StringVar(&fqCmdOutRaw, "fq_command_out", "", "Command line to execute for opening each output file (comma separated)")
	// Arguments: Stats
	var statsInPath, statsOutPath string
	var maxReadLength, maxQual, asciiMin int
	flag.StringVar(&statsInPath, "stats_in_path", "", "Path to statistics of input FASTQ files")
	flag.StringVar(&statsOutPath, "stats_out_path", "", "Path to statistics of output FASTQ files")
	flag.IntVar(&maxReadLength, "max_read_length", 1000, "Maximum read length. No precision required ")
	flag.IntVar(&maxQual, "max_quality", 43, "Highest nucleotide quality")
	flag.IntVar(&asciiMin, "ascii_min", 33, "ASCII coding of the minimum quality score (Ex: 33 for Phred+33)")
	// Arguments
	var opsR1Raw, opsR2Raw, opsR1Path, opsR2Path string
	flag.StringVar(&opsR1Raw, "ops_r1", "", "Operation(s) for read1")
	flag.StringVar(&opsR2Raw, "ops_r2", "", "Operation(s) for read2")
	flag.StringVar(&opsR1Path, "ops_r1_path", "", "Path to operation(s) for read1")
	flag.StringVar(&opsR2Path, "ops_r2_path", "", "Path to operation(s) for read2")
	// Arguments: Parsing
	flag.Parse()

	// Version
	if printVersion {
		fmt.Println(version)
		os.Exit(0)
	}

	// Verbose
	if verbose && verboseLevel == 0 {
		verboseLevel = 1
	}

	// Time start
	var timeStart time.Time
	if verboseLevel > 0 {
		timeStart = time.Now()
	}

	// FASTQs
	var fastqsR1, fastqsR2 []string
	paired := false
	if fqFnamesR1 != "" {
		fastqsR1 = strings.Split(fqFnamesR1, ",")
	}
	if fqFnamesR2 != "" {
		fastqsR2 = strings.Split(fqFnamesR2, ",")
		paired = true
	}

	// Is there work to do?
	if len(fastqsR1) == 0 && len(fastqsR2) == 0 {
		log.Fatal("No input file.")
	}

	// Does output directory exist?
	if fqPathOut != "" {
		if _, err := os.Stat(fqPathOut); os.IsNotExist(err) {
			log.Fatalln(fqPathOut, "not found")
		}
	}

	// Shared parameters
	param := param.Parameters{AsciiMin: asciiMin, MaxQual: maxQual, Paired: paired}

	// Commands
	var fqCmdIn, fqCmdOut []string
	if len(fqCmdInRaw) > 0 {
		fqCmdIn = strings.Split(fqCmdInRaw, ",")
	}
	if len(fqCmdOutRaw) > 0 {
		fqCmdOut = strings.Split(fqCmdOutRaw, ",")
	}

	// Operations
	var opsR1, opsR2 []operations.Operation
	var err error
	if opsR1Raw != "" {
		opsR1, err = operations.ReadOps([]byte(opsR1Raw), param)
		if err != nil {
			log.Fatal(err)
		}
	}
	if opsR2Raw != "" {
		opsR2, err = operations.ReadOps([]byte(opsR2Raw), param)
		if err != nil {
			log.Fatal(err)
		}
	}
	if opsR1Path != "" {
		opsR1, err = operations.ReadOps(readAll(opsR1Path), param)
		if err != nil {
			log.Fatal(err)
		}
	}
	if opsR2Path != "" {
		opsR2, err = operations.ReadOps(readAll(opsR2Path), param)
		if err != nil {
			log.Fatal(err)
		}
	}

	// Apply
	var nPair uint64
	nPair, err = ApplyOperations(fastqsR1, fastqsR2, fqPathOut, fqFnameOutR1, fqFnameOutR2, fqCmdIn, fqCmdOut, opsR1, opsR2, param, statsInPath, statsOutPath, maxReadLength, reportPath, label, bufSize, nWorker, verboseLevel)
	if err != nil {
		log.Fatal(err)
	}

	// Verbose
	if verboseLevel > 0 {
		timeEnd := time.Now()
		fmt.Printf("Done - %.2fmin - %d reads (or pairs)\n", timeEnd.Sub(timeStart).Minutes(), nPair)
	}
}
