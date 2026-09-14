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
#if defined(MPI_TO_K_OMP)
/* See function descriptions in mumps_thread_affinity.h. */
/* { */
#if defined(USE_SCHEDAFFINITY)
/*
 * WARNING: "Some of preprocessor macros, such as _GNU_SOURCE, would
 * have no effect if defined after a system header file has already been
 * included. (Read by <features.h>, and mumps_thread_affinity.h includes
 * "mumps_common.h" which includes "mumps_types_c.h" which include
 * <features.h>). Anyway, keep this define before any header file
 * inclusion.
 */
#define _GNU_SOURCE
#endif
/* } */
#endif
#include "mumps_thread_affinity.h"
#if defined(MPI_TO_K_OMP)
/* { */
#if defined(_OPENMP)
#include <omp.h>
#endif
void MUMPS_CALL
MUMPS_CHECK_OMP_C(MUMPS_INT *NOMP)
{
/* Return 0 if OpenMP is not used, omp_get_max_threads() otherwise */
*NOMP=0;
#if defined(_OPENMP)
*NOMP=omp_get_max_threads();
#endif
}
#if defined(USE_SCHEDAFFINITY)
#include <sched.h>
#endif
#if defined(USE_LIBHWLOC)
#include <hwloc.h>
#if defined(MUMPS_WIN32)
#include <processthreadsapi.h> /* to get a handle to a process */
#include <handleapi.h>         /* to close the handle */
#endif
#endif
#include <stdio.h>
#if !defined(MTKO_NO_MASK_FUSION)
#if defined(USE_LIBHWLOC)||defined(USE_SCHEDAFFINITY)||defined(USE_LIBNUMA)
#include <errno.h>
#if defined(MUMPS_WIN32)
#include <process.h>
#else
#include <unistd.h>
#include <sys/types.h>
#endif
void MUMPS_CALL
MUMPS_GET_PID(MUMPS_INT *pid)
{
  *pid = -1;
  *pid = getpid();
}/* mumps_get_pid */
#endif
#endif
#if defined(USE_LIBHWLOC)
/* { */
/* use of hwloc */
void mumps_hwloc_prepare_topo(hwloc_topology_t *machine_topology, int *error)
{
  int err;
  *error = 0;
  /* create a topology */
  err = hwloc_topology_init(machine_topology);
  if (err < 0) {
#if defined(AFFINITY_VERBOSE)
    printf("failed to initialize the topology\n");
#endif
    *error = -1;
  }
  else
  {
    /* load the topology */
    err = hwloc_topology_load(*machine_topology);
    if (err < 0) {
#if defined(AFFINITY_VERBOSE)
      printf("failed to load the topology\n");
#endif
      *error = -2;
    }
  }/* topology init ok */
}/* mumps_hwloc_prepare_topo */
/* void mumps_hwloc_prepare_set(hwloc_topology_t *topology, hwloc_bitmap_t *set, int *error) */
void mumps_hwloc_prepare_set(hwloc_bitmap_t *set, int *error)
{
  int err;
  *error = 0;
    *set = hwloc_bitmap_alloc();
    if (!*set) {
#if defined(AFFINITY_VERBOSE)
      printf("failed to allocate the bitmap set\n");
#endif
      *error = -3;
    }
    else
    {
      err = hwloc_get_cpubind(machine_topology, *set, HWLOC_CPUBIND_PROCESS);
      if (err < 0) {
#if defined(AFFINITY_VERBOSE)
        printf("failed to get cpu binding\n");
#endif
        *error = -4;
      }
    } /* bitmap alloc ok */
}/* mumps_hwloc_prepare_set */
void mumps_hwloc_clean(hwloc_bitmap_t *set, int status)
{
  if ((status == 0 || status < -3) && *set)
  {
    hwloc_bitmap_free(*set);
  }
  if (status == 0 || status < -1)
  {
    /* hwloc_topology_destroy(*machine_topology); */
  }
}/* mumps_hwloc_clean */
void mumps_hwloc_apply_mask_to_process(hwloc_bitmap_t *set, MUMPS_INT *error)
{
  int err;
  *error=0;
  /* apply the whole mask on the 'master' */
  /* we use the default parameter 0 to rely on the more portable call */
  errno = 0;
  err = hwloc_set_cpubind(machine_topology, *set, 0);
  if (err == -1)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Cannot set cpuset: %i (%s) on current process\n",
        errno, (char*)strerror(errno));
#endif
    *error = -1;
  }
}/* mumps_hwloc_apply_mask_to_process */
void mumps_hwloc_apply_mask_to_threads(hwloc_bitmap_t *set, MUMPS_INT *error)
{
  int err;
#if defined(AFFINITY_VERBOSE)
  int thread_num;
#endif
  *error = 0;
#if defined(AFFINITY_VERBOSE)
#if defined(_OPENMP)
#pragma omp parallel private(err, thread_num)
#endif
#else
#if defined(_OPENMP)
#pragma omp parallel private(err)
#endif
#endif
  {
    errno = 0;
    err = hwloc_set_cpubind(machine_topology, *set, HWLOC_CPUBIND_THREAD);
    if (err == -1)
    {
#if defined(AFFINITY_VERBOSE)
      thread_num = 0;
#if defined(_OPENMP)
      thread_num = omp_get_thread_num();
#endif
      printf("Cannot set cpuset, errno %i (%s) for thread %d\n",
          errno, (char*)strerror(errno), thread_num);
#endif/* AFFINITY_VERBOSE */
      *error = -1;
    }
  }
}/* mumps_hwloc_apply_mask_to_threads */
void mumps_hwloc_apply_mask_to_process_and_threads(hwloc_bitmap_t *whole_set, MUMPS_INT *error)
{
  mumps_hwloc_apply_mask_to_process(whole_set, error);
  if (!(*error))
  {
    /* assign each thread to the whole set - security in case
     * the global mask applied to the master is not enough */
    mumps_hwloc_apply_mask_to_threads(whole_set, error);
  }
}/* mumps_hwloc_apply_mask_to_process_and_threads */
void mumps_apply_process_mask_to_threads(MUMPS_INT *error)
{
  hwloc_bitmap_t proc_set;
  int err;
  *error=0;
  proc_set = hwloc_bitmap_alloc();
  if (!proc_set)
  {
#if defined(AFFINITY_VERBOSE)
    printf("failed to allocate the bitmap set\n");
#endif
    *error = -3;
  }
  else
  {
    err = hwloc_get_cpubind(machine_topology, proc_set, HWLOC_CPUBIND_PROCESS);
    if (err)
    {
#if defined(AFFINITY_VERBOSE)
      printf("Cannot get current cpuset\n");
#endif
      *error = -1;
    }
    else
    {
      mumps_hwloc_apply_mask_to_threads(&proc_set, error);
    }
  }/* proc_set bitmap alloc ok */
}/* mumps_apply_process_mask_to_threads */
void MUMPS_CALL
MUMPS_HWLOC_INIT(MUMPS_INT *error)
{
  mumps_hwloc_prepare_topo(&machine_topology, error);
  if (*error)
  {
    printf("hwloc - error in topology init or load %d\n", *error);
  }
}/* mumps_hwloc_init */
void MUMPS_CALL
MUMPS_HWLOC_TERMINATE()
{
  if (machine_topology)
  {
    hwloc_topology_destroy(machine_topology);
    machine_topology=NULL;
  }
} /* mumps_hwloc_terminate */
void MUMPS_CALL
MUMPS_GET_CPU_ID(MUMPS_INT *cpu_id)
{
  hwloc_bitmap_t set;
  int err;
#if defined(AFFINITY_VERBOSE)
  char str_mask[1024];
#endif
  *cpu_id = -1;
  set = hwloc_bitmap_alloc();
  if (set)
  {
    errno = 0;
    err = hwloc_get_last_cpu_location(machine_topology, set, HWLOC_CPUBIND_THREAD);
    if (err < 0)
    {
#if defined(AFFINITY_VERBOSE)
      printf("failed to get last cpu location, errno: %d\n", errno);
#endif
      *cpu_id = -3;
    }
    else
    {
      *cpu_id = hwloc_bitmap_first(set);
#if defined(AFFINITY_VERBOSE)
      hwloc_bitmap_list_snprintf(str_mask, 1024, set);
      printf("hwloc - GET_CPU_ID - THREAD %d: %s\n", *cpu_id, str_mask);
#endif
    } /* get last cpu location ok */
  }
  else
  {
#if defined(AFFINITY_VERBOSE)
    printf("failed to allocate cpu set\n");
#endif
    *cpu_id = -4;
  }
  hwloc_bitmap_free(set);
}/* get_cpu_id */
void MUMPS_CALL
MUMPS_PRINT_CPU_MASK(MUMPS_INT *id)
{
  hwloc_bitmap_t set;
  char str_mask[1024];
  int err;
  mumps_hwloc_prepare_set(&set, &err);
  if (!err)
  {
    /* everything went well */
    hwloc_bitmap_list_snprintf(str_mask, 1024, set);
    printf("[%d] mask is %s\n", *id, str_mask);
  } /* get cpubind ok */
  mumps_hwloc_clean(&set, err);
}/* print_cpu_mask */
void MUMPS_CALL
MUMPS_PRINT_THREAD_MASK(MUMPS_INT *myid, MUMPS_INT *threadnum, MUMPS_INT *setmask)
{
  hwloc_bitmap_t set;
  char str_mask[1024];
  int err;
  mumps_hwloc_prepare_set(&set, &err);
  if (!err)
  {
    /* everything went well */
    hwloc_bitmap_list_snprintf(str_mask, 1024, set);
    /* printf("[%d] mask is %s\n", *id, str_mask);*/
    if (*setmask == 2)
    {
    /* print after resetting mask */
    printf(" MPI TO K OMP: after restoring mask,  [MPI %d, thread %d] mask is %s\n", (int)*myid, (int)*threadnum, str_mask);
    }
    else if (*setmask == 1)
    {
    /* print after resetting mask */
    printf(" MPI TO K OMP: after resetting mask,  [MPI %d, thread %d] mask is %s\n", (int)*myid, (int)*threadnum, str_mask);
    }
    else
    {
    printf(" MPI TO K OMP: original mask,         [MPI %d, thread %d] mask is %s\n", (int)*myid, (int)*threadnum, str_mask);
    }
  } /* get cpubind ok */
  mumps_hwloc_clean(&set, err);
}/* print_thread_mask */
void MUMPS_CALL
MUMPS_GET_CPU_COUNT(MUMPS_INT *cpu_count)
{
  hwloc_bitmap_t set;
  int err;
  *cpu_count = 0;
  mumps_hwloc_prepare_set(&set, &err);
  if (!err)
  {
    /* everything went well */
    *cpu_count = hwloc_bitmap_weight(set);
  } /* get cpubind ok */
  mumps_hwloc_clean(&set, err);
}/* get_cpu_count */
void MUMPS_CALL
MUMPS_GET_CPU_SET(MUMPS_INT *cpu_count, MUMPS_INT *cpu_set, MUMPS_INT *error)
{
  hwloc_bitmap_t set;
  int err;
  int count, max_num, min_num;
  int i;
  *error=0;
  mumps_hwloc_prepare_set(&set, &err);
  if (!err)
  {
    /* everything went well */
    /* retrieve number of elements in the set */
    count = hwloc_bitmap_weight(set);
    if (count > *cpu_count)
    {
#if defined(AFFINITY_VERBOSE)
      printf("warning cpu_count not enough\n");
#endif/* AFFINITY_VERBOSE */
      *error = -2;
    }
    else
    {
      if (count < *cpu_count)
      {
#if defined(AFFINITY_VERBOSE)
        printf("warning the cpu set contains less cpus than expected\n");
#endif/* AFFINITY_VERBOSE */
        *error = -3;
      }
      else
      {
        max_num = hwloc_bitmap_last(set);
        min_num = hwloc_bitmap_first(set);
        count = 0;
        i = min_num;
        while (i <= max_num && count < *cpu_count)
        {
          if (hwloc_bitmap_isset(set, i))
          {
            cpu_set[count] = i;
            count++;
          }
          i++;
        }
        if (count < *cpu_count)
        {
#if defined(AFFINITY_VERBOSE)
          printf("warning the cpu set contains less cpus than expected\n");
#endif/* AFFINITY_VERBOSE */
          *error = -3;
        }
      }
    }
  } /* get cpubind ok */
  else
  {
    *error = err;
  }
  mumps_hwloc_clean(&set, err);
}/* get_cpu_set */
void MUMPS_CALL
MUMPS_CHECK_CPU_SET_INCLUSION(MUMPS_INT *cpu_count, MUMPS_INT *cpu_set, MUMPS_INT *issub_result)
{
  hwloc_bitmap_t set;
  hwloc_bitmap_t minimal_set;
  int err;
  int i;
  *issub_result = 0;
  mumps_hwloc_prepare_set(&set, &err);
  if (!err)
  {
    /* everything went well */
    /* prepare an other bitmap set to check for inclusion */
    minimal_set = hwloc_bitmap_alloc();
    if (minimal_set)
    {
      for(i = 0; i < *cpu_count; i++)
      {
        hwloc_bitmap_set(minimal_set, i);
      }
      *issub_result = hwloc_bitmap_isincluded(minimal_set, set);
      hwloc_bitmap_free(minimal_set);
    }
    else
    {
#if defined(AFFINITY_VERBOSE)
      printf("failed to allocate the bitmap set\n");
#endif
      *issub_result = -3;
    }
  }
  else
  {
    *issub_result = err;
  }
  mumps_hwloc_clean(&set, err);
}/* check_cpu_set_inclusion */
void MUMPS_CALL
MUMPS_ASSIGN_THREADS_FROM_LIST(MUMPS_INT *cpu_count, MUMPS_INT *threads_nb,
    MUMPS_INT *cpu_set, MUMPS_INT *is_strict, MUMPS_INT *error)
{
  /* hwloc_topology_t topology; */
  hwloc_bitmap_t whole_set;
  hwloc_bitmap_t set;
  int i;
  int err;
  int cores_nb;
#if defined(MTKO_NO_MASK_FUSION)
  int thread_num;
#endif
#if !defined(MTKO_NO_MASK_FUSION)
  int j, jstart, jlast;
  hwloc_pid_t hwlocPid;
#endif
  *error = 0;
  whole_set = hwloc_bitmap_alloc();
  if (!whole_set)
  {
#if defined(AFFINITY_VERBOSE)
    printf("failed to allocate the bitmap set\n");
#endif
    *error = -3;
  }
  else
  {
    if (*is_strict == 0)
    {
      /* retrieve the current cpu set */
      err = hwloc_get_cpubind(machine_topology, whole_set, HWLOC_CPUBIND_PROCESS);
      if (err)
      {
        /* in this case of error, we pursue and we let decide to the caller
         * if it stops or not */
#if defined(AFFINITY_VERBOSE)
        printf("Cannot get current cpuset\n");
#endif
        *error = -1;
      }
      /* set the number of threads */
      cores_nb = *cpu_count;
#if !defined(MTKO_NO_MASK_FUSION)
      cores_nb = cores_nb*(*threads_nb);
#endif
#if defined(_OPENMP)
#if defined(AFFINITY_VERBOSE)
      printf("call to omp_set_num_threads of %d\n", cores_nb);
#endif
      omp_set_num_threads(cores_nb);
#endif
#if !defined(MTKO_NO_MASK_FUSION)
      if (*threads_nb > 0)
      {
        /* prepare the whole set getting the mask of the other processes */
        set = hwloc_bitmap_alloc();
        if (!set)
        {
#if defined(AFFINITY_VERBOSE)
          printf("failed to allocate the bitmap set\n");
#endif
          *error = -3;
        }
        else
        {
          /* we can start at i = 1 as the first pid should match the
             current one and its mask is already in whole_set */
          for(i = 1; i < *cpu_count; i++)
          {
            hwloc_bitmap_zero(set);
            errno = 0;
#if defined(MUMPS_WIN32)
            hwlocPid = OpenProcess(PROCESS_QUERY_INFORMATION, FALSE, cpu_set[i]);
#else
            hwlocPid = cpu_set[i];
#endif
            err = hwloc_get_proc_cpubind(machine_topology, hwlocPid, set, 0);
#if defined(MUMPS_WIN32)
            CloseHandle(hwlocPid);
#endif
            if (err == -1)
            {
              /* in this case of error, we pursue and we let decide to the caller
               * if it stops or not */
#if defined(AFFINITY_VERBOSE)
              printf("Cannot get cpuset for pid %d (err %i, %s)\n",
                  cpu_set[i], errno, (char*)strerror(errno));
#endif
              *error = -1;
            }
            /* browse the set to add the cpu ids to the whole set */
            jstart = hwloc_bitmap_first(set);
            jlast = hwloc_bitmap_last(set);
            j = jstart;
            while (j <= jlast)
            {
              hwloc_bitmap_set(whole_set, j);
              if (j != jlast)
              {
                j = hwloc_bitmap_next(set, j);
              }
              else
              {
                j = jlast + 1;
              }
            }
          }/* loop on the process ids of the other processes */
          /* apply whole set to the process and threads*/
          mumps_hwloc_apply_mask_to_process_and_threads(&whole_set, error);
          hwloc_bitmap_free(set);
        }/* bitmap set alloc ok */
      }/* *threads_nb > 0 */
#else
      /* prepare the whole mask */
      for(i = 0; i < *cpu_count; i++)
      {
        hwloc_bitmap_set(whole_set, cpu_set[i]);
      }
      /* apply whole set to the process */
      mumps_hwloc_apply_mask_to_process(&whole_set, error);
      /* assign each thread to a specific cpu id - security in case
       * the global mask applied to the master is not enough */
#if defined(_OPENMP)
#pragma omp parallel private(set, thread_num, err)
#endif
      {
        set = hwloc_bitmap_alloc();
        if (set)
        {
          thread_num = 0;
#if defined(_OPENMP)
          thread_num = omp_get_thread_num();
#endif
          hwloc_bitmap_zero(set);
          hwloc_bitmap_set(set, cpu_set[thread_num]);
          errno = 0;
          err = hwloc_set_cpubind(machine_topology, set, HWLOC_CPUBIND_THREAD);
          if (err == -1)
          {
#if defined(AFFINITY_VERBOSE)
            printf("Cannot set cpuset %d, errno %i (%s) for thread %d\n",
                cpu_set[thread_num], errno, (char*)strerror(errno), thread_num);
#endif/* AFFINITY_VERBOSE */
            *error = -1;
          }
        }
        else
        {
#if defined(AFFINITY_VERBOSE)
          printf("failed to allocate the bitmap set on thread %d\n",
              thread_num);
#endif
          *error = -3;
        }
        hwloc_bitmap_free(set);
      }
#endif
    }/* not strict mode */
    else
    {
      /* mode strict: the mask is applied to the current process */
      /* prepare the whole mask */
      hwloc_bitmap_zero(whole_set);
      for(i = 0; i < *cpu_count; i++)
      {
        hwloc_bitmap_set(whole_set, cpu_set[i]);
      }
      /* set the number of threads */
#if defined(_OPENMP)
      omp_set_num_threads(*threads_nb);
#endif
      /* apply whole set to the process */
      mumps_hwloc_apply_mask_to_process_and_threads(&whole_set, error);
    }/* strict mode */
    hwloc_bitmap_free(whole_set);
    whole_set = NULL;
  }/* bitmap alloc ok */
}/* assign_threads_from_list (set_mask_from_list) */
#if defined(MTKO_ADVANCED)
void MUMPS_CALL
MUMPS_GET_NUMA_NODE_OF_CPU(MUMPS_INT *cpu_id, MUMPS_INT *numa_node_id)
{
  hwloc_obj_t obj;
  int err;
  *numa_node_id = -1;
  obj = hwloc_get_pu_obj_by_os_index(machine_topology, *cpu_id); 
  /* obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_PU, obj->logical_index); */
  while (obj && obj->type != HWLOC_OBJ_NUMANODE)
  {
    obj = obj->parent;
  }
  if (obj)
  {
    *numa_node_id = obj->os_index;
#if defined(AFFINITY_VERBOSE)
    printf("os index is %d for obj of os index %d\n", *numa_node_id, *cpu_id);
#endif
  }
}/* get_numa_node_of_cpu */
#endif
/* } */
#else
/* { */
#if defined(USE_SCHEDAFFINITY)
/* { */
/* use of linux system calls */
void mumps_sprint_cpu_set(cpu_set_t *mask, char *str_mask, int *len)
{
  int i;
  for(i = 0; i < sizeof(cpu_set_t); i++)
  {
    if (CPU_ISSET(i, mask))
    {
      *len += sprintf(str_mask + (*len), " %d", i);
    }
  }
}/* mumps_print_cpu_set */
void mumps_apply_mask_to_process(cpu_set_t *mask, MUMPS_INT *error)
{
  int err;
   *error=0;
  /* apply the whole mask on the 'master' */
  err = sched_setaffinity(getpid(), sizeof(cpu_set_t), mask);
  if (err)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Cannot set cpuset\n");
#endif/* AFFINITY_VERBOSE */
    *error = -1;
  }
}/* mumps_apply_mask_to_process */
void mumps_apply_mask_to_threads(cpu_set_t *mask, MUMPS_INT *error)
{
  int err;
#if defined(AFFINITY_VERBOSE)
  int thread_num;
#endif
  *error=0;
#if defined(AFFINITY_VERBOSE)
#if defined(_OPENMP)
#pragma omp parallel private(err, thread_num)
#endif
#else
#if defined(_OPENMP)
#pragma omp parallel private(err)
#endif
#endif
  {
    errno = 0;
    err = sched_setaffinity(0, sizeof(cpu_set_t), mask);
    if (err != 0) {
#if defined(AFFINITY_VERBOSE)
      thread_num = 0;
#if defined(_OPENMP)
      thread_num = omp_get_thread_num();
#endif
      printf("Cannot set cpuset, errno %i (%s) for thread %d\n",
          errno, (char*)strerror(errno), thread_num);
#endif/* AFFINITY_VERBOSE */
      *error = -1;
    }
    sched_yield();
  }
}/* mumps_apply_mask_to_process */
void mumps_apply_mask_to_process_and_threads(cpu_set_t *mask, MUMPS_INT *error)
{
  mumps_apply_mask_to_process(mask, error);
  if (!(*error))
  {
    /* assign to each thread mask the whole mask - security in case the
     * global mask applied to the master is not enough */
    mumps_apply_mask_to_threads(mask, error);
  }
}/* mumps_apply_mask_to_process_and_threads */
void mumps_apply_process_mask_to_threads(MUMPS_INT *error)
{
  int err;
  cpu_set_t mask;
  *error=0;
  CPU_ZERO(&mask);
  /* retrieve the current cpu set */
  err = sched_getaffinity(0, sizeof(cpu_set_t), &mask);
  if (err)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Cannot get cpuset, only applying new one.\n");
#endif
    *error = -1;
  }
  else
  {
    mumps_apply_mask_to_threads(&mask, error);
  }
}/* mumps_apply_mask_to_process_and_threads */
void MUMPS_CALL
MUMPS_GET_CPU_ID(MUMPS_INT *cpu_id)
{
  *cpu_id = sched_getcpu();
}/* get_cpu_id */
void MUMPS_CALL
MUMPS_PRINT_CPU_MASK(MUMPS_INT *id)
{
  cpu_set_t mask;
  char str_mask[1024];
  int err;
  int len;
  CPU_ZERO(&mask);
  err = sched_getaffinity(0, sizeof(cpu_set_t), &mask);
  if (err)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Cannot get cpuset\n");
#endif
  }
  else
  {
    len = sprintf(str_mask, "[%d] mask is ", (int)*id);
    mumps_sprint_cpu_set(&mask, str_mask, &len);
    printf("%s\n", str_mask);
  }
}/* print_cpu_mask */
void MUMPS_CALL
MUMPS_PRINT_THREAD_MASK(MUMPS_INT *myid, MUMPS_INT *threadnum, MUMPS_INT *setmask)
{
  cpu_set_t mask;
  char str_mask[1024];
  int err;
  int len;
  CPU_ZERO(&mask);
  err = sched_getaffinity(0, sizeof(cpu_set_t), &mask);
  if (err)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Cannot get cpuset\n");
#endif
  }
  else
  {
     if (*setmask == 2)
    {
    /* print after restoring orginal mask */
    len = sprintf(str_mask, " MPI TO K OMP: after restoring mask,  [MPI %d, thread %d] mask is ", (int)*myid, (int)*threadnum);
    }
     else if (*setmask == 1)
    {
    /* print after resetting mask */
    len = sprintf(str_mask, " MPI TO K OMP: after resetting mask,  [MPI %d, thread %d] mask is ", (int)*myid, (int)*threadnum);
    }
     else
    {
    /* print original mask */
    len = sprintf(str_mask, " MPI TO K OMP: original mask,         [MPI %d, thread %d] mask is ", (int)*myid, (int)*threadnum);
    }
    mumps_sprint_cpu_set(&mask, str_mask, &len);
    printf("%s\n", str_mask);
  }
}/* print_thread_mask */
void MUMPS_CALL
MUMPS_GET_CPU_COUNT(MUMPS_INT *cpu_count)
{
  cpu_set_t mask;
  int err;
  *cpu_count = 0;
  CPU_ZERO(&mask);
  err = sched_getaffinity(getpid(), sizeof(cpu_set_t), &mask);
  if (err)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Cannot get cpuset\n");
#endif/* AFFINITY_VERBOSE */
  }
  else
  {
    *cpu_count = CPU_COUNT(&mask);
  }
}/* get_cpu_count */
void MUMPS_CALL
MUMPS_GET_CPU_SET(MUMPS_INT *cpu_count, MUMPS_INT *cpu_set, MUMPS_INT *error)
{
  cpu_set_t mask;
  int err, i;
  int count;
  *error=0;
  CPU_ZERO(&mask);
  err = sched_getaffinity(getpid(), sizeof(cpu_set_t), &mask);
  if (err)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Cannot get cpuset\n");
#endif/* AFFINITY_VERBOSE */
    *error = -1;
  }
  else
  {
    count = 0;
    for(i = 0; i < sizeof(cpu_set_t); i++)
    {
      if (CPU_ISSET(i, &mask))
      {
        if (count < *cpu_count)
        {
          cpu_set[count] = i;
          count++;
        }
        else
        {
#if defined(AFFINITY_VERBOSE)
          printf("warning cpu_count not enough\n");
#endif/* AFFINITY_VERBOSE */
          *error = -2;
        }
      }
    }
    if (count == *cpu_count)
    {
      *error = 0;
    }
    else
    {
      /* not all the expected cpus have been found in the masks */
      *error = -3;
    }
  }
}/* get_cpu_set */
void MUMPS_CALL
MUMPS_CHECK_CPU_SET_INCLUSION(MUMPS_INT *cpu_count, MUMPS_INT *cpu_set, MUMPS_INT *issub_result)
{
  cpu_set_t curr_mask;
  int err, i;
  int included;
  *issub_result = 0;
  CPU_ZERO(&curr_mask);
  err = sched_getaffinity(getpid(), sizeof(cpu_set_t), &curr_mask);
  if (err)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Cannot get cpuset\n");
#endif/* AFFINITY_VERBOSE */
    *issub_result = -1;
  }
  else
  {
    included = 1;
    i = 0;
    while(included && i < *cpu_count)
    {
      if (!CPU_ISSET(cpu_set[i], &curr_mask))
      {
        included = 0;
      }
      else
      {
        i++;
      }
    }
    *issub_result = included;
  }
}/* check_cpu_set_inclusion */
void MUMPS_CALL
MUMPS_ASSIGN_THREADS_FROM_LIST(MUMPS_INT *cpu_count, MUMPS_INT *threads_nb,
    MUMPS_INT *cpu_set, MUMPS_INT *is_strict, MUMPS_INT *error)
{
  cpu_set_t mask;
  cpu_set_t whole_mask;
  int err, i;
  int cores_nb;
#if !defined(MTKO_NO_MASK_FUSION)
  int j;
#endif
#if defined(MTKO_NO_MASK_FUSION)
  int thread_num;
#endif
  *error = 0;
  /* initializes the whole mask */
  CPU_ZERO(&whole_mask);
  if (*is_strict == 0)
  {
    /* retrieve the current cpu set */
    err = sched_getaffinity(0, sizeof(cpu_set_t), &whole_mask);
    if (err)
    {
      /* in this case of error, we pursue and we let decide to the
       * caller if it stops or not */
#if defined(AFFINITY_VERBOSE)
      printf("Cannot get cpuset, only applying new one.\n");
#endif
      *error = -1;
    }
    /* set the number of threads */
    cores_nb = *cpu_count;
#if !defined(MTKO_NO_MASK_FUSION)
    cores_nb = cores_nb * (*threads_nb);
#else
#endif
#if defined(_OPENMP)
#if defined(AFFINITY_VERBOSE)
    printf("call to omp_set_num_threads of %d\n", cores_nb);
#endif
    omp_set_num_threads(cores_nb);
#endif
#if !defined(MTKO_NO_MASK_FUSION)
    /* prepare the whole set getting the mask of the other processes */
    /* we can start at i = 1 as the first pid should match the current
       one and its mask is already in whole_set */
    for(i = 1; i < *cpu_count; i++)
    {
      CPU_ZERO(&mask);
      err = sched_getaffinity(cpu_set[i], sizeof(cpu_set_t), &mask);
      if (err)
      {
#if defined(AFFINITY_VERBOSE)
        printf("Cannot get cpuset for pid %d, only applying others.\n",
            cpu_set[i]);
#endif
        *error = -1;
      }
      else
      {
        /* adding current cpu set to the whole one */
        for(j = 0; j < sizeof(cpu_set_t); j++)
        {
          if (CPU_ISSET(j, &mask))
          {
            CPU_SET(j, &whole_mask);
          }
        }
      }
    }/* loop on the process ids of the other processes */
    /* apply whole set to the process and threads */
    mumps_apply_mask_to_process_and_threads(&whole_mask, error);
#else
    /* prepare the whole mask */
    for(i = 0; i < *cpu_count; i++)
    {
      CPU_SET(cpu_set[i], &whole_mask);
    }
    /* apply whole set to the process */
    mumps_apply_mask_to_process(&whole_mask, error);
    /* assign each thread to a specific cpu id - security in case the
     * global mask applied to the master is not enough */
#if defined(_OPENMP)
#pragma omp parallel private(mask, thread_num, err)
#endif
    {
      CPU_ZERO(&mask);
      thread_num = 0;
#if defined(_OPENMP)
      thread_num = omp_get_thread_num();
#endif
      CPU_SET(cpu_set[thread_num], &mask);
      errno = 0;
      err = sched_setaffinity(0, sizeof(cpu_set_t), &mask);
      if (err !=0) {
#if defined(AFFINITY_VERBOSE)
        printf("Cannot set cpuset %d, errno %i (%s) for thread %d\n",
            cpu_set[thread_num], errno, (char*)strerror(errno), thread_num);
#endif/* AFFINITY_VERBOSE */
        *error = -1;
      }
      sched_yield();
    }
#endif
  }/* *is_strict == 0 */
  else
  {
    /* prepare the whole mask */
    for(i = 0; i < *cpu_count; i++)
    {
      CPU_SET(cpu_set[i], &whole_mask);
    }
    /* set the number of threads */
#if defined(_OPENMP)
    omp_set_num_threads(*threads_nb);
#endif
    /* apply whole set to the process and threads */
    mumps_apply_mask_to_process_and_threads(&whole_mask, error);
  } /* *is_strict == 1 */
}/* set_mask_from_list */
#if defined(MTKO_ADVANCED)
void MUMPS_CALL
MUMPS_GET_NUMA_NODE_OF_CPU(MUMPS_INT *cpu_id, MUMPS_INT *numa_node_id)
{
  char cmd[64];
  char node_str[32];
  FILE * f;
  *numa_node_id = -1;
  snprintf(cmd, 63, "lscpu -p | grep \"^%d,\" | cut -d, -f4\n", *cpu_id);
  f = popen(cmd, "r");
  if (f)
  {
    fgets(node_str, 32, f);
    if (node_str)
    {
#if defined(AFFINITY_VERBOSE)
      printf("Numa node string get for cpu %d is %s\n", *cpu_id, node_str);
#endif
      errno = 0;
      if (sscanf(node_str, "%d", numa_node_id) != EOF)
      {
#if defined(AFFINITY_VERBOSE)
        printf("numa_node_id of cpu %d is %d\n", *cpu_id, *numa_node_id);
#endif
      }
      else
      {
#if defined(AFFINITY_VERBOSE)
        printf("Unable to read an integer from %s\n", node_str);
#endif
      }
    }
    else
    {
#if defined(AFFINITY_VERBOSE)
      printf("The results has not been retrieved properly\n");
#endif
    }
  }
  else
  {
#if defined(AFFINITY_VERBOSE)
    printf("Cannot execute the command %s\n", cmd);
#endif
  }
}/* get_numa_node_of_cpu */
/* } USE_SCHEDAFFINITY */
#endif
#else
/* { */
#if defined(USE_LIBNUMA)
/* { */
/* use of libnuma */
void MUMPS_CALL
MUMPS_GET_CPU_ID(MUMPS_INT *cpu_id)
{
  *cpu_id = sched_getcpu();
}/* get_cpu_id */
void MUMPS_CALL
MUMPS_PRINT_CPU_MASK(MUMPS_INT *id)
{
  struct bitmask *current;
  int err;
  int i;
  int nb_cpus;
  if (numa_available() < 0)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Your system does not support NUMA API\n");
#endif
  }
  else
  {
    nb_cpus = numa_num_configured_cpus();
    current = numa_bitmask_alloc(nb_cpus);
    if (current)
    {
      numa_sched_getaffinity(getpid(), current);
      for(i = 0; i < nb_cpus; i++)
      {
        if (numa_bitmask_isbitset(current, i)) {
          printf("%d ", i);
        }
      }
      printf("\n");
      numa_bitmask_free(current);
    }
    else
    {
#if defined(AFFINITY_VERBOSE)
      printf("error during numa_bitmask allocation\n");
#endif
    }
  }
}/* print_cpu_mask */
void MUMPS_CALL
MUMPS_PRINT_THREAD_MASK(MUMPS_INT *myid, MUMPS_INT *threadnum, MUMPS_INT *setmask)
{
  struct bitmask *current;
  int err;
  int i;
  int nb_cpus;
  if (numa_available() < 0)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Your system does not support NUMA API\n");
#endif
  }
  else
  {
    nb_cpus = numa_num_configured_cpus();
    current = numa_bitmask_alloc(nb_cpus);
    if (current)
    {
      numa_sched_getaffinity(getpid(), current);
      for(i = 0; i < nb_cpus; i++)
      {
        if (numa_bitmask_isbitset(current, i)) {
          printf("%d ", i);
        }
      }
      printf("\n");
      numa_bitmask_free(current);
    }
    else
    {
#if defined(AFFINITY_VERBOSE)
      printf("error during numa_bitmask allocation\n");
#endif
    }
  }
}/* print_thread_mask */
void MUMPS_CALL
MUMPS_GET_CPU_COUNT(MUMPS_INT *cpu_count)
{
  struct bitmask *current;
  int err;
  int nb_cpus, i;
  if (numa_available() < 0)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Your system does not support NUMA API\n");
#endif
  }
  else
  {
    nb_cpus = numa_num_configured_cpus();
    current = numa_bitmask_alloc(nb_cpus);
    if (current)
    {
      numa_sched_getaffinity(getpid(), current);
      *cpu_count = 0;
      for(i = 0; i < nb_cpus; i++)
      {
        if (numa_bitmask_isbitset(current, i)) {
          *cpu_count = *cpu_count + 1;
        }
      }
      numa_bitmask_free(current);
    }
    else
    {
#if defined(AFFINITY_VERBOSE)
      printf("error during numa_bitmask allocation\n");
#endif
    }
  }
}/* get_cpu_count */
void MUMPS_CALL
MUMPS_GET_CPU_SET(MUMPS_INT *cpu_count, MUMPS_INT *cpu_set, MUMPS_INT *error)
{
  struct bitmask *mask;
  int err, i;
  int count, nb_cpus;
  *error=0;
  if (numa_available() < 0)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Your system does not support NUMA API\n");
#endif
    *error = -1;
  }
  else
  {
    nb_cpus = numa_num_configured_cpus();
    mask = numa_bitmask_alloc(nb_cpus);
    if (mask)
    {
      numa_sched_getaffinity(getpid(), mask);
      count = 0;
      for(i = 0; i < nb_cpus; i++)
      {
        if (numa_bitmask_isbitset(mask,i))
        {
          if (count < *cpu_count)
          {
            cpu_set[count] = i;
            count++;
          }
          else
          {
#if defined(AFFINITY_VERBOSE)
            printf("warning cpu_count not enough\n");
#endif/* AFFINITY_VERBOSE */
            *error = -2;
          }
        }
      }
      if (count == *cpu_count)
      {
        *error = 0;
      }
      else
      {
        /* not all the expected cpus have been found in the maks */
        *error = -3;
      }
      numa_bitmask_free(mask);
    }
    else
    {
#if defined(AFFINITY_VERBOSE)
      printf("error during numa_bitmask allocation\n");
#endif
      *error = -4;
    }
  }
}/* get_cpu_set */
/* WARNING: the following function is probably not up to date */
void MUMPS_CALL
MUMPS_ASSIGN_THREADS_FROM_LIST(MUMPS_INT *cpu_count, MUMPS_INT *threads_nb,
    MUMPS_INT *cpu_set, MUMPS_INT *error)
{
  struct bitmask *mask;
  int err, i;
  int count, nb_cpus;
  *error = 0;
  if (numa_available() < 0)
  {
#if defined(AFFINITY_VERBOSE)
    printf("Your system does not support NUMA API\n");
#endif
    *error = -1;
  }
  else
  {
    nb_cpus = numa_num_configured_cpus();
    mask = numa_bitmask_alloc(nb_cpus);
    if (mask)
    {
      /* numa_sched_getaffinity(getpid(), mask); */
#if defined(_OPENMP)
      omp_set_num_threads(*cpu_count);
#pragma omp parallel for private(mask)
#endif
      for(i=0; i < *cpu_count; i++)
      {
        numa_bitmask_setbit(mask, cpu_set[i]);
        err = numa_sched_setaffinity(0, mask);
/*         if (err != 0)
 *         {
 * #if defined(AFFINITY_VERBOSE)
 *           printf("Cannot set cpuset: %i\n",errno);
 * #endif 
 *         }
 */
        sched_yield();
      }
      for(i=0; i < *cpu_count; i++)
      {
        numa_bitmask_setbit(mask, cpu_set[i]);
      }
      err = numa_sched_setaffinity(0, mask);
      if (err)
      {
#if defined(AFFINITY_VERBOSE)
        printf("Cannot set cpuset\n");
#endif/* AFFINITY_VERBOSE */
        *error = -1;
      }
      numa_bitmask_free(mask);
    }
    else
    {
#if defined(AFFINITY_VERBOSE)
      printf("error during numa_bitmask allocation\n");
#endif
      *error = -4;
    }
  }
}/* assign_threads_from_list (set_mask_from_list) */
void mumps_apply_process_mask_to_threads(MUMPS_INT *error)
{
#if defined(AFFINITY_VERBOSE)
  printf("not yet implemented\n");
#endif
  *error = -1;
}/* mumps_apply_process_mask_to_threads */
/* use libnuma */
/* } */
#else
/* { */
/* use no lib for sched affinity, assume everything goes well with no
 * need to set mask */
#if !defined(MTKO_NO_MASK_FUSION)
void MUMPS_CALL
MUMPS_GET_PID(MUMPS_INT *pid)
{
  *pid = 0;
}/* get_pid */
#endif
void MUMPS_CALL
MUMPS_GET_CPU_ID(MUMPS_INT *cpu_id)
{
  *cpu_id = 0;
}/* get_cpu_id */
void MUMPS_CALL
MUMPS_PRINT_CPU_MASK(MUMPS_INT *id)
{
}/* print_cpu_mask */
void MUMPS_CALL
MUMPS_PRINT_THREAD_MASK(MUMPS_INT *myid, MUMPS_INT *threadnum, MUMPS_INT *setmask)
{
}/* print_thread_mask */
void MUMPS_CALL
MUMPS_GET_CPU_COUNT(MUMPS_INT *cpu_count)
{
  *cpu_count = 0;
}/* get_cpu_count */
void MUMPS_CALL
MUMPS_GET_CPU_SET(MUMPS_INT *cpu_count, MUMPS_INT *cpu_set, MUMPS_INT *error)
{
  *error = 0;
}/* get_cpu_set */
void MUMPS_CALL
MUMPS_CHECK_CPU_SET_INCLUSION(MUMPS_INT *cpu_count, MUMPS_INT *cpu_set, MUMPS_INT *issub_result)
{
  *issub_result = 1;
}/* check_cpu_set_inclusion */
void MUMPS_CALL
MUMPS_ASSIGN_THREADS_FROM_LIST(MUMPS_INT *cpu_count, MUMPS_INT *threads_nb,
    MUMPS_INT *cpu_set, MUMPS_INT *is_strict, MUMPS_INT *error)
{
  *error=0;
#if defined(AFFINITY_VERBOSE)
  printf("setting %d omp threads\n", *cpu_count);
#endif
  if (*is_strict == 0)
  {
#if defined(MTKO_NO_MASK_FUSION)
#if defined(_OPENMP)
    omp_set_num_threads(*cpu_count);
#endif
#else
#if defined(_OPENMP)
    omp_set_num_threads(*cpu_count*(*threads_nb));
#endif
#endif
  }
  else
  {
#if defined(_OPENMP)
    omp_set_num_threads(*threads_nb);
#endif
  }
}/* assign_threads_from_list */
void mumps_apply_process_mask_to_threads(MUMPS_INT *error)
{
  /* nothing to do with this configuration */
  *error = 0;
}/* mumps_apply_process_mask_to_threads */
#if defined(MTKO_ADVANCED)
void MUMPS_GET_NUMA_NODE_OF_CPU(MUMPS_INT *cpu_id, MUMPS_INT *numa_node_id)
{
  *numa_node_id = 0;
}/* mumps_get_numa_node_of_cpu */
#endif
/* } */
#endif
/* } */
#endif
/* } */
#endif
void MUMPS_CALL
MUMPS_MTKO_DEACTIVATE(MUMPS_INT *threads_nb)
{
#if defined(_OPENMP)
#if defined(AFFINITY_VERBOSE)
  printf("[MTKO_DEACTIVATE] => call to omp_set_num_threads of %d\n",
      *threads_nb);
#endif
  omp_set_num_threads(*threads_nb);
#endif
}/* mumps_mtko_deactivate */
void MUMPS_CALL
MUMPS_MTKO_REACTIVATE(MUMPS_INT *threads_nb, MUMPS_INT *error)
{
#if defined(_OPENMP)
#if defined(AFFINITY_VERBOSE)
  printf("[MTKO_REACTIVATE] => call to omp_set_num_threads of %d\n",
      *threads_nb);
#endif
  omp_set_num_threads(*threads_nb);
#endif
#if defined(AFFINITY_VERBOSE)
  printf("[MTKO_REACTIVATE] => set back the affinity mask to the threads");
#endif
  *error = 0;
  mumps_apply_process_mask_to_threads(error);
}/* mumps_mtko_reactivate */
/* } */
#else /* MPI_TO_K_OMP */
void MUMPS_CALL
MUMPS_THREAD_AFFINITY_RETURN()
{
  /*
   * Thread affinity tools will be available in the future.
   */
}
#endif /* MPI_TO_K_OMP */
