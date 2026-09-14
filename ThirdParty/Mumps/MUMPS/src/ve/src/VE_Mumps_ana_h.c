/*
 *  -----------------------------------------------------------------------
 *  NEC Vector Host MUMPS Library Interface
 *  -----------------------------------------------------------------------
 */

#include <VE_Mumps.h>

/*
 *  ---------------------------------------------------------------------------
 *  VE_Mumps_ana_h
 *  ---------------------------------------------------------------------------
 */

void VE_Mumps_ana_h( int TOTEL       , bool COMPUTE_PERM , int N         , int64_t IWLEN , int64_t * PE ,
                     int64_t * PFREE , int * LEN         , int * IW      , int * NV      , int * ELEN   ,
                     int * LAST      , int * NCMPA       , int * DEGREE  , int * HEAD    , int * NEXT   ,
                     int * W         , int * PARENT      )
{
  int ret=1;
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
  
  symid = vhcall_find(handle, "mumps_ana_h_");
  if ( symid == -1 ) {
    fprintf( stdout , "%s @ %d failed : can't find symbol MUMPS_ANA_H\n" , basename(__FILE__) , __LINE__ );
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
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , 0 , &TOTEL , sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , 1 , &COMPUTE_PERM , sizeof( bool ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , 2 , &N , sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_IN , 3 , &IWLEN , sizeof( int64_t ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_INOUT , 4 , PE , N * sizeof( int64_t ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_INOUT , 5 , PFREE , sizeof( int64_t ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_INOUT , 6 , LEN , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_INOUT , 7 , IW , IWLEN * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_INOUT , 8 , NV , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_OUT , 9 , ELEN , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_OUT ,10 , LAST , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_OUT ,11 , NCMPA , sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_INOUT ,12 , DEGREE , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_INOUT ,13 , HEAD , TOTEL * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_INOUT ,14 , NEXT , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }
  
  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_INOUT ,15 , W , N * sizeof( int ) );
  if ( ret ) {
    perror("vhcall_args_set_pointer");
    fprintf( stdout , "%s @ %d failed : unable to set argument.\n" , basename(__FILE__) , __LINE__ );
    fflush( stdout );
    (void) mumps_abort();
  }

  ret = vhcall_args_set_pointer( ca , VHCALL_INTENT_OUT ,16 , PARENT , N * sizeof( int ) );
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
