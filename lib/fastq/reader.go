//
// Copyright (C) 2017-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package fastq

import (
	"bufio"
	"io"
	"os"
	"os/exec"
)

type FqReader struct {
	Fos    *os.File
	Pipe   io.ReadCloser
	Reader *bufio.Reader
	Done   bool
}

func Ropen(fpath string, cmd []string, bufSize int) (*FqReader, error) {
	fq := new(FqReader)
	var err error

	// Check input file exists
	if _, err := os.Stat(fpath); os.IsNotExist(err) {
		return fq, err
	}

	if len(cmd) == 0 {
		if fq.Fos, err = os.Open(fpath); err != nil {
			return fq, err
		}
		fq.Reader = bufio.NewReaderSize(fq.Fos, bufSize)
	} else {
		cmd = append(cmd, fpath)
		p := exec.Command(cmd[0], cmd[1:]...)
		if fq.Pipe, err = p.StdoutPipe(); err != nil {
			return fq, err
		}
		if err = p.Start(); err != nil {
			return fq, err
		}
		fq.Reader = bufio.NewReaderSize(fq.Pipe, bufSize)
	}

	return fq, nil
}

// Close closes Ropen
func (fq *FqReader) Close() {
	if fq.Fos != nil {
		fq.Fos.Close()
	}
	if fq.Pipe != nil {
		fq.Pipe.Close()
	}
}

func (fq *FqReader) Iter() (Record, error) {
	var n, s, q []byte
	var err error
	if n, err = fq.Reader.ReadBytes('\n'); err != nil {
		if err == io.EOF {
			fq.Done = true
			return Record{}, nil
		} else {
			return Record{}, err
		}
	}
	if s, err = fq.Reader.ReadBytes('\n'); err != nil {
		return Record{}, err
	}
	if _, err = fq.Reader.ReadBytes('\n'); err != nil {
		return Record{}, err
	}
	if q, err = fq.Reader.ReadBytes('\n'); err != nil {
		return Record{}, err
	}
	return Record{Name: n[1 : len(n)-1], Seq: s[:len(s)-1], Qual: q[:len(q)-1]}, nil
}
