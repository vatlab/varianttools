


#include "hdf5.h"
#include "dirent.h"
#include "regex.h"
#include "string.h"
//#include "../hdf5-blosc/blosc_filter.h"
#include "blosc_filter.h"
#include "genotypes.h"
#include "unistd.h"
#include "stdbool.h"


int findIndex(int *array, size_t size, int target)
{
    int i=0;
    while((i<size) && (array[i] != target)) i++;

    return (i<size) ? (i) : (-1);
}




int get_Genotype_from_hdf5(char* filePath, char* chrName, int variant_id, int *samples,int index, char* genotypeData ,char* conditionalNodeValue, int* sample_IDs, int numberOfSamples)
{
            hid_t       file;                        /* handles */
            hid_t       genotype;
            hid_t       mask;
            hid_t       rowname;
            hid_t       colname;
            hid_t       maskSpace;
            hid_t       genoSpace;
            hid_t       rowSpace;
            hid_t       colSpace;
            hid_t       mask_memspace;
            hid_t       geno_memspace;
            hid_t       row_memspace;
            hid_t       col_memspace;
            hid_t       cparms;
            hsize_t     dims[2];                     /* dataset and chunk dimensions*/

            hsize_t     row_dims[1];
            hsize_t     col_dims[1];
            hsize_t     count[2];
            hsize_t     offset[2];
            herr_t      status, status_n;


            int         rank, row_rank,col_rank,rank_chunk;
            hsize_t i, j;

            char *version, *date;
            // int r;
            // r = register_blosc(&version, &date);
            // printf("Blosc version info: %s (%s) %d\n", version, date,r);
            file = H5Fopen(filePath, H5F_ACC_RDONLY, H5P_DEFAULT);


            char row_name[20]="";
            strcpy(row_name, chrName);
            strcat(row_name, "/rownames");

            char col_name[20]="";
            strcpy(col_name, chrName);
            strcat(col_name, "/colnames");

            char maskData[20]="";
            strcpy(maskData, chrName);
            strcat(maskData, "/Mask");


            genotype = H5Dopen(file, genotypeData, H5P_DEFAULT);
            genoSpace = H5Dget_space(genotype);    /* Get filespace handle first. */
            rank      = H5Sget_simple_extent_ndims(genoSpace);
            status_n  = H5Sget_simple_extent_dims(genoSpace, dims, NULL);
            // printf("genotype rank %d, dimensions %lu x %lu\n",
            //    rank, (unsigned long)(dims[0]), (unsigned long)(dims[1]));


            mask = H5Dopen(file, maskData, H5P_DEFAULT);
            maskSpace = H5Dget_space(mask);



            rowname = H5Dopen(file, row_name, H5P_DEFAULT);
            rowSpace = H5Dget_space(rowname);    /* Get filespace handle first. */
            row_rank      = H5Sget_simple_extent_ndims(rowSpace);
            status_n  = H5Sget_simple_extent_dims(rowSpace, row_dims, NULL);
            // printf("rownames rank %d, dimensions %lu\n",
            //    row_rank, (unsigned long)(row_dims[0]));


            colname = H5Dopen(file, col_name, H5P_DEFAULT);
            colSpace = H5Dget_space(colname);    /* Get filespace handle first. */
            col_rank      = H5Sget_simple_extent_ndims(colSpace);
            status_n  = H5Sget_simple_extent_dims(colSpace, col_dims, NULL);
            // printf("colnames rank %d, dimensions %lu\n",
            //    col_rank, (unsigned long)(col_dims[0]));



            row_memspace = H5Screate_simple(row_rank,row_dims,NULL);
            int         rownames[row_dims[0]];  /* buffer for dataset to be read */
            status = H5Dread(rowname, H5T_NATIVE_INT, row_memspace, rowSpace,
                     H5P_DEFAULT, rownames);

            col_memspace = H5Screate_simple(col_rank,col_dims,NULL);
            int         colnames[col_dims[0]];  /* buffer for dataset to be read */
            status = H5Dread(colname, H5T_NATIVE_INT, col_memspace, colSpace,
                     H5P_DEFAULT, colnames);


            int pos=findIndex(rownames,row_dims[0],variant_id);

            // printf("row position %d.\n", pos);
            if (pos!=-1){

                int numberOfColumns=dims[1];
                offset[0] = pos;
                offset[1] = 0;
                count[0]  = 1;
                count[1]  = numberOfColumns;
                geno_memspace = H5Screate_simple(1,col_dims,NULL);
                int         genotypes[numberOfColumns];
                status = H5Sselect_hyperslab(genoSpace, H5S_SELECT_SET, offset, NULL,
                             count, NULL);
                status = H5Dread(genotype, H5T_NATIVE_INT, geno_memspace, genoSpace,
                         H5P_DEFAULT, genotypes);


                mask_memspace = H5Screate_simple(1,col_dims,NULL);
                int         masks[numberOfColumns];
                status = H5Sselect_hyperslab(maskSpace, H5S_SELECT_SET, offset, NULL,
                             count, NULL);
                status = H5Dread(mask, H5T_NATIVE_INT, mask_memspace, maskSpace,
                         H5P_DEFAULT, masks);
                // printf("Dataset: \n");
                // for (i = 0; i < numberOfColumns; i++) printf("%d ", genotypes[i]);


                // printf("\n");
                // printf("rownames: \n");
                // for (i = 0; i<row_dims[0]; i++) printf("%d ",rownames[i]);

                // printf("\n");
                // printf("colnames: \n");

                for (i = 0; i<numberOfColumns; i++) {
                    if (!conditionalNodeValue){
                        if (sample_IDs[0]==-1){
                            if (genotypes[i]>=0 && masks[i]>0){
                                samples[index]=colnames[i];
                                // printf("%d ",colnames[i]);
                                index+=1;
                            }
                        }else{
                            if (findIndex(sample_IDs,numberOfSamples,colnames[i])!=-1){
                                if (masks[i]>0){
                                    samples[index]=genotypes[i];
                                    index+=1;
                                }
                            }

                        }
                    }else{
                        if (genotypes[i]==atoi(conditionalNodeValue) && masks[i]>0){
                            samples[index]=colnames[i];
                            // printf("%d ",colnames[i]);
                            index+=1;
                        }
                    }
                }
            }



            H5Pclose(cparms);
            H5Dclose(genotype);
            H5Dclose(mask);
            H5Dclose(rowname);
            H5Dclose(colname);
            H5Sclose(genoSpace);
            H5Sclose(maskSpace);
            H5Sclose(rowSpace);
            H5Sclose(colSpace);
            H5Sclose(mask_memspace);
            H5Sclose(geno_memspace);
            H5Sclose(row_memspace);
            H5Sclose(col_memspace);


            H5Fclose(file);
            return index;
}



void get_Genotypes(char* chr,int variant_id,int* samples,int numberOfSamples, char* genoFilter, int* sample_IDs)
{

    regex_t     regex;
    DIR *dir;
    // char *dirName="/Users/jma7/Development/VAT_ref/ismb-2018/data/";
    char dirName[200];
    if (getcwd(dirName, sizeof(dirName)) != NULL) {
       // printf("Current working dir: %s\n", dirName);
   } else {
       perror("getcwd() error");
   }

    struct dirent *ent;
    int reti = regcomp(&regex, "^tmp.*genotypes.h5$", REG_EXTENDED);
    // static int samples[300];
    // int samples[numberOfSamples];
    // samples=malloc(numberOfSamples * sizeof(int));
    memset( samples, -1, numberOfSamples*sizeof(int) );
    int index=0;
    char chrName[20]="";
    strcpy(chrName, "/chr");
    strcat(chrName,chr);

    char genotypeData[20]="";
    strcpy(genotypeData, chrName);
    char *conditionalNodeValue=NULL;

    if (genoFilter != NULL){
        if (sample_IDs[0]==-1){
            char delim[]="=";
            int stringPos=0;
            char *ptr = strtok(genoFilter, delim);
            char *array[2];
            while(ptr != NULL)
            {
                array[stringPos++]=ptr;
                ptr = strtok(NULL, delim);
            }
            char nodeName[20]="";
            strcpy(nodeName, "/");
            strcat(nodeName,array[0]);
            strcat(genotypeData, nodeName);
            // printf("%s\n",genotypeData);
            // printf("%s\n",array[1]);
            conditionalNodeValue=array[1];
        }else{
            char nodeName[20]="";
            strcpy(nodeName, "/");
            strcat(nodeName,genoFilter);
            strcat(genotypeData, nodeName);
        }
    }else{
        strcat(genotypeData, "/GT");
    }

    if ((dir = opendir (dirName)) != NULL) {
      while ((ent = readdir (dir)) != NULL) {
        reti = regexec(&regex, ent->d_name, 0, NULL, 0);
        if (!reti){
            char * filePath = (char *) malloc(2 + strlen(dirName)+ strlen(ent->d_name) );
            strcpy(filePath, dirName);
            strcat(filePath,"/");
            strcat(filePath, ent->d_name);

            // printf("%s,%s,%d\n",filePath,chr,variant_id);
            index=get_Genotype_from_hdf5(filePath,chrName,variant_id,samples,index, genotypeData, conditionalNodeValue, sample_IDs,numberOfSamples);
            // printf("done with file %s ,index %d \n",filePath,index);
        }
      }
      closedir (dir);

    } else {
      perror ("");
      return EXIT_FAILURE;
    }

}


// int main(void){
//     int *samples;
//     samples=get_Genotypes("8",1473,20);
//     int i;
//     printf("%s\n", "hah");
//     for (i = 0; i<10; i++) printf("%d ",samples[i]);
//     return 0;
// }
