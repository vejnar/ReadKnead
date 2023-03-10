//
// Copyright © 2017 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package operations

import (
	"errors"
	"fmt"

	"git.sr.ht/~vejnar/ReadKnead/lib/fastq"

	"github.com/buger/jsonparser"
)

type Clip struct {
	name         string
	end          int
	length       int
	addClipped   bool
	addSeparator bool
}

func NewClip(data []byte) (*Clip, error) {
	c := Clip{name: "clip"}
	end, err := jsonparser.GetInt(data, "end")
	if err != nil {
		if errors.Is(err, jsonparser.KeyPathNotFoundError) {
			return &c, fmt.Errorf("%w: %v", err, "end")
		}
		return &c, err
	} else {
		c.end = int(end)
	}
	length, err := jsonparser.GetInt(data, "length")
	if err != nil {
		if errors.Is(err, jsonparser.KeyPathNotFoundError) {
			return &c, fmt.Errorf("%w: %v", err, "length")
		}
		return &c, err
	} else {
		c.length = int(length)
	}
	addClipped, err := jsonparser.GetBoolean(data, "add_clipped")
	if err == jsonparser.KeyPathNotFoundError {
		c.addClipped = false
	} else if err != nil {
		return &c, err
	} else {
		c.addClipped = addClipped
	}
	addSeparator, err := jsonparser.GetBoolean(data, "add_separator")
	if err == jsonparser.KeyPathNotFoundError {
		c.addSeparator = true
	} else if err != nil {
		return &c, err
	} else {
		c.addSeparator = addSeparator
	}
	return &c, nil
}

func (op *Clip) Name() string {
	return op.name
}

func (op *Clip) IsThreadSafe() bool {
	return true
}

func (op *Clip) GetDpx(idx int) ([][]byte, int) {
	return [][]byte{}, idx
}

func (op *Clip) Transform(p *fastq.ExtPair, r int, ot *OpStat, verboseLevel int) int {
	if op.end == 5 {
		if r == 1 {
			if verboseLevel > 2 {
				fmt.Printf("%s %s r%d\n%s\n", op.name, p.R1.Name, r, p.R1.Seq)
			}
			if len(p.R1.Seq) < op.length {
				ot.OpsR1[op.name]["too_short"]++
				return 1
			} else {
				if op.addClipped {
					if op.addSeparator {
						p.R1.Name = append(p.R1.Name, byte('#'))
						p.R2.Name = append(p.R2.Name, byte('#'))
					}
					p.R1.Name = append(p.R1.Name, p.R1.Seq[:op.length]...)
					p.R2.Name = append(p.R2.Name, p.R1.Seq[:op.length]...)
				}
				p.R1.Seq = p.R1.Seq[op.length:]
				p.R1.Qual = p.R1.Qual[op.length:]
				if verboseLevel > 2 {
					fmt.Printf("> %s\n", p.R1.Seq)
				}
				return 0
			}
		} else {
			if verboseLevel > 2 {
				fmt.Printf("%s %s r%d\n%s\n", op.name, p.R2.Name, r, p.R2.Seq)
			}
			if len(p.R2.Seq) < op.length {
				ot.OpsR2[op.name]["too_short"]++
				return 1
			} else {
				if op.addClipped {
					if op.addSeparator {
						p.R1.Name = append(p.R1.Name, byte('#'))
						p.R2.Name = append(p.R2.Name, byte('#'))
					}
					p.R1.Name = append(p.R1.Name, p.R2.Seq[:op.length]...)
					p.R2.Name = append(p.R2.Name, p.R2.Seq[:op.length]...)
				}
				p.R2.Seq = p.R2.Seq[op.length:]
				p.R2.Qual = p.R2.Qual[op.length:]
				if verboseLevel > 2 {
					fmt.Printf("> %s\n", p.R2.Seq)
				}
				return 0
			}
		}
	} else if op.end == 3 {
		if r == 1 {
			if verboseLevel > 2 {
				fmt.Printf("%s %s r%d\n%s\n", op.name, p.R1.Name, r, p.R1.Seq)
			}
			clipIndex := len(p.R1.Seq) - op.length
			if clipIndex < 0 {
				ot.OpsR1[op.name]["too_short"]++
				return 1
			} else {
				if op.addClipped {
					if op.addSeparator {
						p.R1.Name = append(p.R1.Name, byte('#'))
						p.R2.Name = append(p.R2.Name, byte('#'))
					}
					p.R1.Name = append(p.R1.Name, p.R1.Seq[clipIndex:]...)
					p.R2.Name = append(p.R2.Name, p.R1.Seq[clipIndex:]...)
				}
				p.R1.Seq = p.R1.Seq[:clipIndex]
				p.R1.Qual = p.R1.Qual[:clipIndex]
				if verboseLevel > 2 {
					fmt.Printf("> %s\n", p.R1.Seq)
				}
				return 0
			}
		} else {
			if verboseLevel > 2 {
				fmt.Printf("%s %s r%d\n%s\n", op.name, p.R2.Name, r, p.R2.Seq)
			}
			clipIndex := len(p.R2.Seq) - op.length
			if clipIndex < 0 {
				ot.OpsR2[op.name]["too_short"]++
				return 1
			} else {
				if op.addClipped {
					if op.addSeparator {
						p.R1.Name = append(p.R1.Name, byte('#'))
						p.R2.Name = append(p.R2.Name, byte('#'))
					}
					p.R1.Name = append(p.R1.Name, p.R2.Seq[clipIndex:]...)
					p.R2.Name = append(p.R2.Name, p.R2.Seq[clipIndex:]...)
				}
				p.R2.Seq = p.R2.Seq[:clipIndex]
				p.R2.Qual = p.R2.Qual[:clipIndex]
				if verboseLevel > 2 {
					fmt.Printf("> %s\n", p.R2.Seq)
				}
				return 0
			}
		}
	}
	return 0
}
