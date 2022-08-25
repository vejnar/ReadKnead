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

type FqWriter struct {
	Fos    *os.File
	Pipe   io.WriteCloser
	Writer *bufio.Writer
}

func Wopen(fpath string, cmd []string, bufSize int) (*FqWriter, error) {
	fq := new(FqWriter)
	var err error

	if len(cmd) == 0 {
		if fq.Fos, err = os.Create(fpath); err != nil {
			return fq, err
		}
		fq.Writer = bufio.NewWriterSize(fq.Fos, bufSize)
	} else {
		cmd = append(cmd, fpath)
		p := exec.Command(cmd[0], cmd[1:]...)
		if fq.Pipe, err = p.StdinPipe(); err != nil {
			return fq, err
		}
		if err = p.Start(); err != nil {
			return fq, err
		}
		fq.Writer = bufio.NewWriterSize(fq.Pipe, bufSize)
	}

	return fq, nil
}

func (fq *FqWriter) WriteRecord(r Record) error {
	var err error
	if err = fq.Writer.WriteByte(byte('@')); err != nil {
		return err
	}
	if _, err = fq.Writer.Write(r.Name); err != nil {
		return err
	}
	if err = fq.Writer.WriteByte(byte('\n')); err != nil {
		return err
	}
	if _, err = fq.Writer.Write(r.Seq); err != nil {
		return err
	}
	if _, err = fq.Writer.Write([]byte{'\n', '+', '\n'}); err != nil {
		return err
	}
	if _, err = fq.Writer.Write(r.Qual); err != nil {
		return err
	}
	if err = fq.Writer.WriteByte(byte('\n')); err != nil {
		return err
	}
	return nil
}

// Close closes Ropen
func (fq *FqWriter) Close() error {
	err := fq.Writer.Flush()
	if err != nil {
		return err
	}
	if fq.Pipe != nil {
		fq.Pipe.Close()
	}
	if fq.Fos != nil {
		fq.Fos.Close()
	}
	return nil
}
