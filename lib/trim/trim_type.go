//
// Copyright (C) 2017-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//

package trim

type TrimType int

const (
	NoTrimType TrimType = iota
	TrimExactType
	TrimAlignType
	TrimTooShortType
)

func (t TrimType) String() string {
	return []string{"no_trim", "trim_exact", "trim_align", "trim_too_short"}[t]
}

var TrimTypes = []TrimType{NoTrimType, TrimExactType, TrimAlignType, TrimTooShortType}
