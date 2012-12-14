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



#ifndef __FILE_H__
#define __FILE_H__

/**
 * Defines return values from file operations.
 */
enum file_status {
    /**
     * Means that an operation was successfully performed.
     */
    FILE_OK,

    /**
     * Means that an error occurred.
     */
    FILE_ERROR
};

typedef enum file_status file_status_t;

/**
 * Copies a file to another.
 *
 * @param from_path The path to copy from.
 * @param to_path   The destination path.
 *
 * @return FILE_OK if the file was copied, FILE_ERROR otherwise.
 */
file_status_t file_copy(const char *from_path, const char *to_path);

/**
 * Removes the file of the given path
 *
 * @param path Path to the file to remove.
 *
 * @return FILE_OK if the file was removed, FILE_ERROR otherwise.
 */
file_status_t file_remove(const char *path);

#endif /* __FILE_H__ */
