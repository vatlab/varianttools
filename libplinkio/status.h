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



#ifndef __STATUS_H__
#define __STATUS_H__

enum pio_status_e
{
    /**
     * Function successful.
     */
    PIO_OK,

    /**
     * File reached EOF.
     */
    PIO_END,

    /**
     * Generic error.
     */
    PIO_ERROR
};

typedef enum pio_status_e pio_status_t;

#endif /* End of __STATUS_H__ */
