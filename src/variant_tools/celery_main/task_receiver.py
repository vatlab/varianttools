from variant_tools.celery_main.start_celery import app
import time
import random
import numpy as np
from variant_tools.accessor import *
from celery import Celery






class HDF5GenotypeImportWorker:
    def __init__(self, chunk,variantIndex,start_sample,end_sample,sample_ids, 
        proc_index,dbLocation,genotype_info,build):


        self.chunk=chunk
        self.variantIndex = variantIndex
        # self.variant_count = variant_count
        self.proc_index = proc_index
        self.start_sample=start_sample
        self.end_sample=end_sample
        
        self.indptr=[]
        self.indices=[]
        self.data=[]
        self.rownames=[]
        self.sample_ids=sample_ids
        self.firstID=0
        if sample_ids[0]!=1 and start_sample!=0:
            self.firstID=sample_ids[0]
            self.start_sample=self.start_sample+1
            self.end_sample=self.end_sample+1
        self.colnames=[self.sample_ids[i-self.firstID] for i in range(self.start_sample,self.end_sample)]
        self.genoCount=0
        self.dbLocation=dbLocation
        self.build=build
        self.info={}
        self.rowData=[]
        self.info["GT"]=[]
        self.info["Mask"]=[]
        self.namedict={}
        self.geno_info=[]
        if "GT" in genotype_info:
            genotype_info.remove("GT")
        if len(genotype_info)>0:
            for info in genotype_info:
                #indptr,indices,data,shape,rownames
                if not isinstance(info,str):
                    self.geno_info.append(info.name.replace("_geno",""))
                else:
                    self.geno_info.append(info.replace("_geno",""))
        for info in self.geno_info:
            self.info[info]=[]
            if "calldata/"+info in self.chunk and np.nansum(self.chunk["calldata/"+info][:10])>0:
                    self.namedict[info]="calldata/"+info
            elif "variants/"+info in self.chunk and np.nansum(self.chunk["variants/"+info][:10])>0:
                    self.namedict[info]="variants/"+info

  
    # check io_vcf_read.pyx function vcf_genotype_parse to see the meaning of coding
    def get_geno(self,variant_id,pos,altIndex):
        self.rownames.append(variant_id)
        # print(self.dbLocation,self.start_sample,self.end_sample,self.firstID)

        if "calldata/GT" in self.chunk:

            GT=self.chunk["calldata/GT"][pos,self.start_sample-self.firstID:self.end_sample-self.firstID]
            GT=GT.astype(float)
            if altIndex==0:
                GT[np.logical_or(GT==3, GT==4)]=np.nan          
            elif altIndex==1:
                # GT_geno[GT_geno==3]=1
                # GT_geno[GT_geno==4]=2
                GT[(GT!=3)&(GT!=4)&(GT!=-1)]=np.nan
                # GT_geno[np.logical_and(GT_geno!=3, GT_geno!=4)]=np.nan
                GT[GT==3]=1
                GT[GT==4]=2
            GT[GT==-10]=np.nan
            self.info["GT"].append(GT)
            self.info["Mask"].append([1.0]*len(GT))
        else:
            # GT_geno=[np.nan]
            GT=[-1]
            self.info["GT"].append(GT)
            self.info["Mask"].append([1.0]*len(GT))
      
        if len(self.geno_info)>0:
            # self.rowData.extend([[variant_id,idx,self.chunk["calldata/DP"][i][idx],self.chunk["calldata/GQ"][i][idx]] for idx in range(self.start_sample,self.end_sample)])
            # self.rowData.extend([[variant_id,idx]+[self.chunk[field][i][idx] for field in self.fields] for idx in range(self.start_sample,self.end_sample)])
            # self.getInfoTable(variant_id,infoDict,altIndex)
            for info in self.geno_info:
                if "variants" in self.namedict[info]:
                    self.info[info].append(np.array([self.chunk[self.namedict[info]][pos]]))
                else:
                    self.info[info].append(self.chunk[self.namedict[info]][pos,self.start_sample-self.firstID:self.end_sample-self.firstID])
                # print(self.namedict[info],info,self.start_sample,self.end_sample,pos,(self.chunk[self.namedict[info]][pos,self.start_sample-self.firstID:self.end_sample-self.firstID]))
        



    def writeIntoHDF(self,chr):
        storageEngine=Engine_Storage.choose_storage_engine(self.dbLocation)
        shape=np.array([len(self.rownames),len(self.colnames)])
        storageEngine.store(np.array(self.info["GT"]),chr,"GT")
        storageEngine.store(np.array(self.info["Mask"]),chr,"Mask")
        storageEngine.store(np.array(self.rownames),chr,"rownames")
        rowmask=np.zeros(len(self.rownames),dtype=np.bool)
        storageEngine.store(np.array(rowmask),chr,"rowmask")
        
        if not storageEngine.checkGroup(chr,"colnames"):
            storageEngine.store(np.array(self.colnames),chr,"colnames")
            colmask=np.zeros(len(self.colnames),dtype=np.bool)
            storageEngine.store(np.array(colmask),chr,"samplemask")
        
        storageEngine.store(shape,chr,"shape")

        self.info["GT"]=[]
        self.info["Mask"]=[]

        if len(self.geno_info)>0:
            for info in self.geno_info:
                storageEngine.store_genoInfo(np.array(self.info[info]),chr,info)
                self.info[info]=[]
 
        storageEngine.close()       
        self.rownames=[]
 

   
    def run(self):
        
        prev_chr=self.chunk["variants/CHROM"][0].replace("chr","")
        prev_variant_id=-1     
        for i in range(len(self.chunk["variants/ID"])):
            infoDict={}
            chr=self.chunk["variants/CHROM"][i].replace("chr","")
            if chr!=prev_chr:
                self.writeIntoHDF(prev_chr)
                prev_chr=chr             
            ref=self.chunk["variants/REF"][i]
            pos=self.chunk["variants/POS"][i]
            for altIndex in range(len(self.chunk["variants/ALT"][i])):
                alt=self.chunk["variants/ALT"][i][altIndex]
                if alt!="":
                    if tuple((chr, ref, alt)) in self.variantIndex:
                        variant_id  = self.variantIndex[tuple((chr, ref, alt))][pos][0]
                        
                        if variant_id!=prev_variant_id:
                            self.get_geno(variant_id,i,altIndex)
                            prev_variant_id=variant_id
                        
                    else:
                        rec=[str(chr),str(pos),ref,alt]  
                        msg=normalize_variant(RefGenome(self.build).crr, rec, 0, 1, 2, 3)
                        if tuple((rec[0], rec[2], rec[3])) in self.variantIndex:
                            variant_id  = self.variantIndex[tuple((rec[0], rec[2], rec[3]))][rec[1]][0]
                            if variant_id!=prev_variant_id:
                                self.get_geno(variant_id,i,altIndex)
                                prev_variant_id=variant_id
        self.writeIntoHDF(chr)


@app.task(bind=True,default_retry_delay=10,serializer="pickle")

def do_work(self, chunk,variantIndex,start_sample,end_sample,sample_ids, 
        proc_index,dbLocation,genotype_info,build):
        worker=HDF5GenotypeImportWorker(chunk,variantIndex,start_sample,end_sample,sample_ids,proc_index,dbLocation,genotype_info,build)
        worker.run()

