//
// Copyright (C) 2017-2022 Charles E. Vejnar
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://www.mozilla.org/MPL/2.0/.
//

package bio

var NucleotidesDNA = []byte{'A', 'C', 'G', 'T', 'N', 'M', 'R', 'S', 'V', 'W', 'Y', 'H', 'K', 'D', 'B', 'X', 'a', 'c', 'm', 'g', 'r', 's', 'v', 't', 'w', 'y', 'h', 'k', 'd', 'b', 'n', 'x'}

func IsDNA(nt byte) bool {
	for _, b := range NucleotidesDNA {
		if b == nt {
			return true
		}
	}
	return false
}
