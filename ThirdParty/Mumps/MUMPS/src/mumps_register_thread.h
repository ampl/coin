/*
 *
 *  This file is part of MUMPS 5.9.1, released
 *  on Mon Jul 20 09:00:43 UTC 2026
 *
 *
 *  Copyright 1991-2026 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
#ifndef MUMPS_REGISTER_THREAD_H
#define MUMPS_REGISTER_THREAD_H
#if defined(USE_GPU)
/*
 * Request pool initialization.
 */
void mumps_init_request_pool(int size_pool);
/*
 * Request pool finalization.
 */
void mumps_finalize_request_pool();
/*
 * Asynchronous host register.
 */
uint64_t mumps_host_register_async(void *ptr, size_t size);
/*
 * Asynchronous host unregister.
 */
uint64_t mumps_host_unregister_async(void *ptr, size_t size);
/*
 * Wait for a given pending request.
 */
int mumps_request_wait(uint64_t index);
/*
 * Wait for all pending requests.
 */
int mumps_request_wait_all();
#else
void mumps_register_thread_return();
{
  /* Registering tools will be available in the future. */
}
#endif
#endif /* MUMPS_REGISTER_THREAD_H */
