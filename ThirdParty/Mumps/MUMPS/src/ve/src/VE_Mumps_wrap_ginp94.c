/*
 *  -----------------------------------------------------------------------
 *  NEC Vector Host MUMPS Library Interface
 *  -----------------------------------------------------------------------
 */

#include <VE_Mumps.h>

/*
 *  ---------------------------------------------------------------------------
 *  VE_Mumps_wrap_ginp94
 *  ---------------------------------------------------------------------------
 */

void VE_Mumps_wrap_ginp94( int N               , int64_t * IPE , int * IW             , int64_t   LIW8  , int * PERM    ,
                           int * SizeOfBlocks  , int KEEP60    , int * LISTVAR_SCHUR  , int SIZE_SCHUR  , int KEEP378   ,
                           int * COLCOUNT      , int * PARENT  , int * PORDER         , int * IWTMP1    , int * IWTMP2  ,
                           int * IWTMP3        , int * IWTMP4  , int * IWTMP5         , int * INFO      )
{
  int ret=1,iarg;
  int64_t symid;
  vhcall_handle handle;
  vhcall_args *ca;
  uint64_t retval;
  
  char * env_metis_lib_name;
  char metis_lib_name[ NECVH_MAX_LIBNAME ];

  /*
   *  -----------------------------------------------------------------------
   *  Set the library name ( check environment )
   *  -----------------------------------------------------------------------
   */
  
  if ( env_metis_lib_name = getenv( "VH_MUMPS_LIBRARY" ) ) {
    size_t lenv = strlen( env_metis_lib_name );
    if ( strlen( env_metis_lib_name ) > (size_t)NECVH_MAX_LIBNAME ) {
      fprintf( stdout , "%s @ %d failed : metis library name too long (%s)\n" , basename(__FILE__) , __LINE__ , env_metis_lib_name );
      fflush( stdout );
      (void) mumps_abort();
    }
    else {
      memcpy( ( void * ) metis_lib_name , ( void * ) env_metis_lib_name , lenv * sizeof( char) );
      metis_lib_name[ lenv ] = '\0';
      fprintf( stdout , "%s @ %d input library %s\n" , basename(__FILE__) , __LINE__ , metis_lib_name );
      fflush( stdout );
    }
  }  
  else {
    strcpy( metis_lib_name , "libvh.so" );
  }
  
  /*
   *  -----------------------------------------------------------------------
   *  Load VH C library
   *  -----------------------------------------------------------------------
   */
  
  handle = vhcall_install( metis_lib_name );
  if (handle == (vhcall_handle)-1) {
    perror("vhcall_install");
    fprintf( stdout , "%s @ %d failed : can't install library %s\n" , basename(__FILE__) , __LINE__ , metis_lib_name );
    fprintf( stdout , "Please check your LD_LIBRARY_PATH variable\n");
    fflush( stdout );
    (void) mumps_abort();
  }
  
  /*
   *  -----------------------------------------------------------------------
   *  Find VH C library function
   *  -----------------------------------------------------------------------
   */
  
  symid = vhcall_find(handle, "mumps_wrap_ginp94_");
  if ( symid == -1 ) {
    fprintf( stdout , "%s @ %d failed : can't find symbol MUMPS_WRAP_GINP94\n" , basename(__FILE__) , __LINE__ );
    perror("vhcall_find");
    fflush( stdout );
    (void) mumps_abort();
  }
  
  /*
   *  -----------------------------------------------------------------------
   *  Allocate arguments object for VH library function
   *  -----------------------------------------------------------------------
   */
  
  ca = vhcall_args_alloc();
  if ( ! ca ) {
    perror("vhcall_args_alloc");
    fprintf( stdout , "%s @ %d failed : unable to allocate function's argument list\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  /*
   *  -----------------------------------------------------------------------
   *  Set arguments
   *  -----------------------------------------------------------------------
   */
   
  iarg = 0;
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , iarg++ , &N , sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , iarg++ , IPE , ( N + 1 ) * sizeof( int64_t ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , iarg++ , IW , LIW8 * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , iarg++ , &LIW8 , sizeof( int64_t ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_INOUT , iarg++ , PERM , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , iarg++ , SizeOfBlocks , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , iarg++ , &KEEP60 , sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , iarg++ , LISTVAR_SCHUR , SIZE_SCHUR * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , iarg++ , &SIZE_SCHUR , sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , iarg++ , &KEEP378 , sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_OUT , iarg++ , COLCOUNT , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_OUT , iarg++ , PARENT , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_OUT , iarg++ , PORDER , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_OUT , iarg++ , IWTMP1 , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_OUT , iarg++ , IWTMP2 , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }

  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_OUT , iarg++ , IWTMP3 , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }

  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_OUT , iarg++ , IWTMP4 , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }

  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_OUT , iarg++ , IWTMP5 , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }

  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_INOUT , iarg++ , INFO , 2 * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  /*
   *  -----------------------------------------------------------------------
   *  Invoke VH C library function
   *  -----------------------------------------------------------------------
   */
  
  ret = vhcall_invoke_with_args( symid , ca , &retval );
  if ( ret ) {
    perror("vhcall_invoke_with_args");
    fprintf( stdout , "%s @ %d failed : unable to invoke the symbol from library\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }

  /*
   *  -----------------------------------------------------------------------
   *  Free args
   *  -----------------------------------------------------------------------
   */
  
  vhcall_args_free( ca );
  
  /*
   *  -----------------------------------------------------------------------
   *  Uninstall the library
   *  -----------------------------------------------------------------------
   */

  if ( vhcall_uninstall( handle ) ) {
    perror("vhcall_uninstall");
    fprintf( stdout , "cvhmetis_nodend failed : unable to uninstall the library\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
}
