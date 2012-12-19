/* =====================================================================================
// 
//  This is a small C and Python library for reading Plink genotype files,
//  written by Mattias Franberg, version 0.2.2 
//  
//  https://bitbucket.org/mattias_franberg/libplinkio
//
//  This software is not licensed or copyrighted. The varianttools developers
//  have been contacting its author and will include the license information when we
//  hear from the author, or replace it with alternative implementation if the author
//  requests for a removal.
// 
 ===================================================================================== */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "file.h"

/**
 * Total size for a buffer containing the copy command
 * that will be issued to the shell.
 */
#define FILE_COPY_BUFFER_SIZE 4096

file_status_t
file_copy(const char *from_path, const char *to_path)
{
    char *copy_command = malloc( sizeof( char ) * ( strlen( from_path ) + strlen( to_path ) + 5 ) );

    sprintf( copy_command, "cp %s %s", from_path, to_path );
    int status = system( copy_command );
    free( copy_command );

    if( status != -1 )
    {
        return FILE_OK;
    }
    else
    {
        return FILE_ERROR;
    }

}

file_status_t
file_remove(const char *path)
{
    char *rm_command = malloc( sizeof( char ) * ( strlen( path ) + 4 ) ); 
    
    sprintf( rm_command, "rm %s", path );
    int status = system( rm_command );
    free( rm_command );

    if( status != -1 )
    {
        return FILE_OK;
    }
    else
    {
        return FILE_ERROR;
    }
}
