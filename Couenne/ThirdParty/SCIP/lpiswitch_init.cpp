// Copyright (C) 2012 Stefan Vigerske and others
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Stefan Vigerske
//
// the only reason for not including this into lpiswitch.c is that static initializations by functions are not allowed in C

#include "config_scip.h"
#undef SCIP_VERSION

#include "lpiswitch.h"

static SCIP_RETCODE lpiswinitreturn = SCIPlpiSwitchSetSolver(SCIP_LPISW_LPDEFAULT);
