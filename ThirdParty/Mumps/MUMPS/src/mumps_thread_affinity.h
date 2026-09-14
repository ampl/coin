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
#ifndef MUMPS_THREAD_AFFINITY_H
#define MUMPS_THREAD_AFFINITY_H
#include "mumps_common.h"
#if defined(MPI_TO_K_OMP)
#define MUMPS_CHECK_OMP_C \
  F_SYMBOL(check_omp_c,CHECK_OMP_C)
void MUMPS_CALL
MUMPS_CHECK_OMP_C(MUMPS_INT *NOMP);
/* This function returns 0 if OpenMP is not used, omp_get_max_threads() otherwise */
#if defined(USE_LIBHWLOC)
#include <hwloc.h>
/**
 * The following global structure is used to store the topology of the
 * machine.
 * In the situation of use of the MPI_TO_K_OMP functionality, it is not
 * expected to need thread safety.
 * If needed someday, requiring a mutex before writing the topology and
 * checking it has not been taken each time we need to read it should be
 * enough.
 */
hwloc_topology_t machine_topology;
#define MUMPS_HWLOC_INIT \
  F_SYMBOL(hwloc_init,HWLOC_INIT)
void MUMPS_CALL
MUMPS_HWLOC_INIT(MUMPS_INT *error);
#define MUMPS_HWLOC_TERMINATE \
  F_SYMBOL(hwloc_terminate,HWLOC_TERMINATE)
void MUMPS_CALL
MUMPS_HWLOC_TERMINATE();
#endif
#if !defined(MTKO_NO_MASK_FUSION)
/**
 * MUMPS_GET_PID returns the process id of the current process.
 *
 * Output:
 * pid            the process id
 */
#define MUMPS_GET_PID \
  F_SYMBOL(get_pid,GET_PID)
void MUMPS_CALL
MUMPS_GET_PID(MUMPS_INT *pid);
#endif
/**
 * MUMPS_GET_CPU_ID returns the current (last) physical cpu id of the running
 * job or thread. (last because at the time we will use or display it,
 * it might have changed already !)
 *
 * Output:
 * cpu_id         the last physical cpu id
 */
#define MUMPS_GET_CPU_ID \
  F_SYMBOL(get_cpu_id,GET_CPU_ID)
void MUMPS_CALL
MUMPS_GET_CPU_ID(MUMPS_INT *cpu_id);
/**
 * MUMPS_PRINT_CPU_MASK displays in the standard output the current
 * allowed cpu set.
 *
 * Input:
 * id             an id used by the caller program to identify the
 *                calling process or thread or whatever.
 *                The cpu set string displayed is prefixed by [id].
 */
#define MUMPS_PRINT_CPU_MASK \
  F_SYMBOL(print_cpu_mask,PRINT_CPU_MASK)
void MUMPS_CALL
MUMPS_PRINT_CPU_MASK(MUMPS_INT *id);
/**
 * MUMPS_PRINT_THREAD displays in the standard output the current
 * allowed cpu set.
 *
 * Input:
 * myid           the mpi id of the calling process
 * threadnum      the openmp thread number of the calling process
 * setmask        to indicate whether it is called after mask restoring (=2), modification (=1) 
                  or original mask (not equal to 1 or 2)
 */
#define MUMPS_PRINT_THREAD_MASK \
  F_SYMBOL(print_thread_mask,PRINT_THREAD_MASK)
void MUMPS_CALL
MUMPS_PRINT_THREAD_MASK(MUMPS_INT *myid, MUMPS_INT *threadnum, MUMPS_INT *setmask);
/**
 * MUMPS_GET_CPU_COUNT returns the number of elements in the current
 * allowed cpu set.
 *
 * Output:
 * cpu_count      the number of elements in the current cpu set
 */
#define MUMPS_GET_CPU_COUNT \
  F_SYMBOL(get_cpu_count,GET_CPU_COUNT)
void MUMPS_CALL
MUMPS_GET_CPU_COUNT(MUMPS_INT *cpu_count);
/**
 * MUMPS_GET_CPU_SET fills in an array of size cpu_count with the
 * cpu_count physical cpu id in the current allowed cpu set.
 *
 * Input:
 * cpu_count      the number of elements expected in the set (size of
 *                the previously allocated cpu_set array)
 *
 * Output:
 * cpu_set        array of size cpu_count containing the indices of
 *                the physical cores in the current cpu set
 * error          if 0 everything went well
 *                if < 0 an error has occured.
 */
#define MUMPS_GET_CPU_SET \
  F_SYMBOL(get_cpu_set,GET_CPU_SET)
void MUMPS_CALL
MUMPS_GET_CPU_SET(MUMPS_INT *cpu_count, MUMPS_INT *cpu_set, MUMPS_INT *error);
/**
 * MUMPS_CHECK_CPU_SET_INCLUSION checks if a cpu set is included into
 * the current cpu set.
 *
 * Input:
 * cpu_count      the number of elements in the set
 * cpu_set        the list of physical cores to use
 *
 * Output:
 * issub_result   the result
 *                if < 0 an error has occured
 *                if 0 cpu_set is not included in the current cpu set
 *                if 1 cpu_set is included in the current cpu set
 */
#define MUMPS_CHECK_CPU_SET_INCLUSION \
  F_SYMBOL(check_cpu_set_inclusion,CHECK_CPU_SET_INCLUSION)
void MUMPS_CALL
MUMPS_CHECK_CPU_SET_INCLUSION(MUMPS_INT *cpu_count, MUMPS_INT *cpu_set, MUMPS_INT *issub_result);
/**
 * MUMPS_ASSIGN_THREADS_FROM_LIST applies a cpu set to the current
 * process.
 *
 * cpu_count      the number of elements in the set
 * cpu_set        the list of physical cores to use
 * is_strict      if set to 1, strict mode is applied
 * error          output parameter to warn if an error has occured
 *                (if != 0)
 *
 * the number of threads is set to the number of elements in the set
 * (cpu_count), but in mode non strict with fusion (default) it is set
 * to the number of elements in the set multiplied by the number of
 * threads
 *
 * 2 modes are available:
 *  - mode strict: (is_strict is equal to 1)
 *                 the calling process and all threads apply the mask
 *  - mode non strict: (is_strict is equal to 0)
 *                 in fusion (default), a mask is built using the
 *                 affinity mask of all processes identified by their
 *                 process ids given in the cpu_set. The first pid is
 *                 the current one.
 *                 This mask is applied to the process then to all the
 *                 threads.
 *
 *                 In non fusion, the initial mask of the calling
 *                 process is augmented with the physical cores of the
 *                 cpu set transmitted in cpu_set.
 *                 Then all threads are assigned a mask containing a
 *                 value of the cpu set.
 *
 */
#define MUMPS_ASSIGN_THREADS_FROM_LIST \
  F_SYMBOL(assign_threads_from_list,ASSIGN_THREADS_FROM_LIST)
void MUMPS_CALL
MUMPS_ASSIGN_THREADS_FROM_LIST(MUMPS_INT *cpu_count, MUMPS_INT *threads_nb,
    MUMPS_INT *cpu_set, MUMPS_INT *is_strict, MUMPS_INT *error);
/** MUMPS_MTKO_DEACTIVATE restore the number of threads that a MPI
 * process used to run before the use of the MPI TO K OMP feature.
 *
 * threads_nb     the number of threads to set back
 *
 *
 * MUMPS_MTKO_DEACTIVATE and MUMPS_MTKO_REACTIVATE may be called to use
 * the MPI TO K OMP feature only on MT computational intensive parts of
 * the code. Call _REACTIVATE before the MT computational intensive part
 * and call _DEACTIVE after.
 */
#define MUMPS_MTKO_DEACTIVATE \
  F_SYMBOL(mtko_deactivate,MTKO_DEACTIVATE)
void MUMPS_CALL
MUMPS_MTKO_DEACTIVATE(MUMPS_INT *threads_nb);
/** MUMPS_MTKO_REACTIVATE set back a number of threads and applies the
 * affinity mask of the calling process to all the threads.
 *
 * threads_nb     the new number of threads
 *
 */
#define MUMPS_MTKO_REACTIVATE \
  F_SYMBOL(mtko_reactivate,MTKO_REACTIVATE)
void MUMPS_CALL
MUMPS_MTKO_REACTIVATE(MUMPS_INT *threads_nb, MUMPS_INT *error);
#if defined(MTKO_ADVANCED)
/**
 * MUMPS_GET_NUMA_NODE_OF_CPU returns the index of the numa node of a
 * given cpu index.
 *
 * Input:
 * cpu_id        index of a physical core
 *
 * Output:
 * numa_node_id  index of the corresponding numa node
 *               values -1 in case of error
 */
#define MUMPS_GET_NUMA_NODE_OF_CPU \
  F_SYMBOL(get_numa_node_of_cpu,GET_NUMA_NODE_OF_CPU)
void MUMPS_CALL
MUMPS_GET_NUMA_NODE_OF_CPU(MUMPS_INT *cpu_id, MUMPS_INT *numa_node_id);
#endif
#endif
#define MUMPS_THREAD_AFFINITY_RETURN \
  F_SYMBOL(thread_affinity_return,THREAD_AFFINITY_RETURN)
void MUMPS_CALL
MUMPS_THREAD_AFFINITY_RETURN();
#endif /* MUMPS_THREAD_AFFINITY_H */
