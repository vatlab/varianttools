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



#ifndef __BIM_PARSE_H__
#define __BIM_PARSE_H__

#include <status.h>

/**
 * Parses the loci and points the given locus array to a
 * the memory that contains them, and writes back the number
 * of loci.
 *
 * @param bim_file Bim file.
 *
 * @return PIO_OK if the loci could be parsed, PIO_ERROR otherwise.
 */
pio_status_t parse_loci(FILE *bim_fp, UT_array *locus);

#endif /* End of __BIM_PARSE_H__ */
