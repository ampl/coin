// Copyright (C) 2012 Stefan Vigerske and others
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Stefan Vigerske

/* add Soplex2 suffix to LPI functions */
#define LPISW_RENAME_SUFFIX Soplex2
#include "lpiswitch_rename.h"

// lpi_spx and lpi_spx2 define a class SPxSCIP, so that the linker may not know which one to call
// thus, rename the class for the Soplex2 interface
#define SPxSCIP SPx2SCIP

// lpi_spx2 does not work when SOPLEX_LEGACY is set
// to keep things simple, just fall back to use lpi_spx in this case, even though it makes the Soplex2 lpiswitch option equivalent to Soplex
#ifdef SOPLEX_LEGACY
#include "lpi/lpi_spx.cpp"
#else
#include "lpi/lpi_spx2.cpp"
#endif
