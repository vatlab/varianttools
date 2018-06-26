from variant_tools.celery_main.start_celery import app
import time
import random
import re
import glob
import numpy as np
from variant_tools.accessor import *
from variant_tools.tester import *
from celery import Celery
# from variant_tools.association_hdf5 import generateHDFbyGroup,getGenotype_HDF5,generateHDFbyGroup_update
from variant_tools.assoTests import AssoData
from variant_tools.utils import DatabaseEngine,executeUntilSucceed





class HDF5GenotypeImportWorker:
    def __init__(self, chunk,variantIndex,start_sample,end_sample,sample_ids, 
        proc_index,dbLocation,genotype_info,build):


        self.chunk=chunk
        self.variantIndex = variantIndex
        # self.variant_count = variant_count
        self.proc_index = proc_index
        self.start_sample=start_sample
        self.end_sample=end_sample
        self
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



class AssoTestsWorker:
    '''Association test calculator'''
    def __init__(self, param, grp, result_fields,methods,path):

        self.param = param
        self.proj = param.proj
        self.table = param.table
        self.sample_IDs = param.sample_IDs
        self.phenotypes = param.phenotypes
        self.covariates = param.covariates
        self.phenotype_names = param.phenotype_names
        self.covariate_names = param.covariate_names
        self.var_info = param.var_info
        self.geno_info = param.geno_info
        self.tests = self.getAssoTests(methods,len(self.covariates),[])
        self.group_names = param.group_names
        self.missing_ind_ge = param.missing_ind_ge
        self.missing_var_ge = param.missing_var_ge
        self.sample_names = param.sample_names
        

        self.num_extern_tests = param.num_extern_tests
        self.grp = grp
        self.result_fields = result_fields
        self.path=path
   
        # self.shelves = {}
        #

        self.db = DatabaseEngine()
        self.db.connect(self.path+"/"+param.proj.name+'.proj',readonly=True)

        # self.db.attach(param.proj.name + '.proj', '__fromVariant', lock=self.shelf_lock) 
        
        #
        self.g_na = float('NaN')
        if env.treat_missing_as_wildtype:
            self.g_na = 0.0
        # self.args=args


    # def __del__(self):
    #     # self.db.close()
    #     for val in list(self.shelves.values()):
    #         val.close()

    
    def filterGenotype(self, genotype, geno_info, var_info, gname):
        '''
        Filter genotypes for missing calls or lack of minor alleles. Not very efficient because 
        it copies genotype, var_info and geno_info skipping the loci to be removed.
            - genotype is a Individual_list * Variants_list matrix of genotype values
            - var_info is a dictionary with each key being information corresponding Variant_list
            - geno_info is a dictionary with each key having a matrix of the same structure as genotype matrix
        '''
        # Step 1: filter individuals by genotype missingness at a locus

       
        missing_ratios = [sum(list(map(math.isnan, x))) / float(len(x)) for x in genotype]
        which = [x < self.missing_ind_ge for x in missing_ratios]
        # check for non-triviality of phenotype data
        if sum(which) < 5:
            raise ValueError("Sample size too small ({0}) to be analyzed for {1}.".format(sum(which), repr(gname)))
        if len(which) - sum(which) > 0:
            env.logger.debug('In {}, {} out of {} samples will be removed due to '
                              'having more than {}% missing genotypes'.\
                              format(repr(gname), len(which) - sum(which), len(which),
                                     self.missing_ind_ge * 100))
        # Step 2: filter variants by genotype missingness at a locus
        keep_loci = []
        for i in range(len(genotype[0])):
            # tag individuals missing variant calls
            missingness_vi = list(map(math.isnan, [x[i] for x, y in zip(genotype, which) if y]))
            # unique genotype codings on the locus
            gt_codings = list(set([x[i] for x, y in zip(genotype, which) if y and not math.isnan(x[i])]))
            keep_loci.append((float(sum(missingness_vi)) / float(len(missingness_vi))) < self.missing_var_ge and len(gt_codings) > 1)
        if len(keep_loci) - sum(keep_loci) > 0:
            for idx in range(len(genotype)):
                # filter genotype and geno_info
                genotype[idx] = [i for i, j in zip(genotype[idx], keep_loci) if j]
                for k in list(geno_info.keys()):
                    geno_info[k][idx] = [i for i, j in zip(geno_info[k][idx], keep_loci) if j]
            # filter var_info
            for k in list(var_info.keys()):
                var_info[k] = [i for i, j in zip(var_info[k], keep_loci) if j]
            #
            env.logger.debug('In {}, {} out of {} loci will be removed due to '
                              'having no minor allele or having more than {}% missing genotypes'.\
                              format(repr(gname), len(keep_loci) - sum(keep_loci),
                                     len(keep_loci), self.missing_ind_ge * 100))
        # check for non-triviality of genotype matrix
        if len(genotype[0]) == 0:
            raise ValueError("No variant found in genotype data for {}.".format(repr(gname)))
            # raise ValueError("No variant found in genotype data for {}.".format(repr(gname)))
        return genotype, which, var_info, geno_info

    def getAssoTests(self, methods, ncovariates, common_args):
        '''Get a list of methods from parameter methods, passing method specific and common
        args to its constructor. This function sets self.tests as a list of statistical tests'''
        if not methods:
            raise ValueError('Please specify at least one statistical test. '
                             'Please use command "vtools show tests" for a list of tests')
        tests = []
        for m in methods:
            name = m.split()[0]
            args = m.split()[1:] + common_args
            try:
                if '.' in name:
                    # if the method is defined elsewhere
                    m_module, m_name = name.split('.', 1)
                    # also search current working directory
                    my_dir = os.getcwd()
                    env.logger.info('Loading {} from {}'.format(m_module, my_dir))
                    if my_dir not in sys.path:
                        sys.path.append(my_dir)
                        # use the default level, which is -1 for python 2 and 0 for python 3
                        _temp = __import__(m_module, globals(), locals(), [m_name])
                        sys.path.pop()
                    else:
                        _temp = __import__(m_module, globals(), locals(), [m_name])
                    env.logger.info('Loading {}'.format(m_module))
                    method = getattr(_temp, m_name)(ncovariates, args)
                else:
                    method = eval(name)(ncovariates, args)
                # check if method is valid
                if not hasattr(method, 'fields'):
                    raise ValueError('Invalid association test method {}: '
                                     'missing attribute fields'.format(name))
                if not method.fields:
                    env.logger.warning('Association test {} has invalid or empty fields. '
                                        'No result will be generated.'.format(name))
                tests.append(method)
            except Exception as e:
                raise ValueError('Failed to load association test {0}: {1}.'
                                 'Please use command "vtools show tests" to list usable tests'.format(name, e))
        return tests



    def setGenotype(self, which, data, info, grpname):
        geno = [x for idx, x in enumerate(data) if which[idx]]
        self.data.setGenotype(geno)
        self.data.setVar("gname", str(grpname))
        for field in list(info.keys()):
            self.data.setVar('__geno_' + field, [x for idx, x in enumerate(info[field]) if which[idx]])

    def setPhenotype(self, which):
        '''Set phenotype data'''
        if len(self.phenotypes) > 1:
            raise ValueError('Only a single phenotype is allowed at this point')
        # print(len(self.phenotypes[0]))
        # for idx, x in enumerate(self.phenotypes[0]):
        #     if which[idx]:
        #         print(idx,x)
        phen = [x for idx, x in enumerate(self.phenotypes[0]) if which[idx]]
        if self.covariates:
          covt = [[x for idx, x in enumerate(y) if which[idx]] for y in self.covariates]
        if self.covariates:
          self.data.setPhenotype(phen, covt)
        else:
          self.data.setPhenotype(phen)

    def setVarInfo(self, data):
        for field in list(data.keys()):
            if field not in ['chr', 'pos']:
                self.data.setVar('__var_' + field, data[field])



    def getVarInfo(self, group, where_clause):
        var_info = {x:[] for x in self.var_info}
        query = 'SELECT variant_id {0} FROM __asso_tmp WHERE ({1})'.format(
            ','+','.join([x.replace('.', '_') for x in self.var_info]) if self.var_info else '', where_clause)
        #env.logger.debug('Running query: {}'.format(query))
        cur = self.db.cursor()
        # SELECT can fail when the disk is slow which causes database lock problem.
        msg = 'Load variant info for group {} using association worker'.format(group)
        executeUntilSucceed(cur, query, 5, msg, group)
        #
        if not self.var_info:
            data = {x[0]:[] for x in cur.fetchall()}
        else:
            data = {x[0]:x[1:] for x in cur.fetchall()}
        variant_id = sorted(list(data.keys()), key=int)
        for idx, key in enumerate(self.var_info):
            if key not in ['variant.chr', 'variant.pos']:
                var_info[key] = [data[x][idx] if (type(data[x][idx]) in [int, float]) else float('NaN') for x in variant_id]
            else:
                var_info[key] = [data[x][idx] for x in variant_id] 
        return var_info, variant_id

    def setPyData(self, which, geno, var_info, geno_info,
                  missing_code, grpname, recode_missing = True):
        '''set all data to a python dictionary'''
        def fstr(x):
            try:
                float(x)
            except:
                x = str(x)
            return x
        #
        if len(self.phenotypes) > 1:
            raise ValueError('Only a single phenotype is allowed at this point')
        #
        self.pydata['name'] = grpname
        #
        if recode_missing:
            self.pydata['genotype'] = [[missing_code if math.isnan(e) else e for e in x] for idx, x in enumerate(geno) if which[idx]]
        else:
            self.pydata['genotype'] = [x for idx, x in enumerate(geno) if which[idx]]
        #
        try:
            self.pydata['coordinate'] = [(str(x), str(y)) for x, y in zip(var_info['variant.chr'], var_info['variant.pos'])]
        except:
            self.pydata['coordinate'] = []
        # var_info
        self.pydata['var_info'] = []
        self.pydata['var_info_header'] = []
        for k, item in list(var_info.items()):
             if k != 'variant.chr' and k != 'variant.pos':
                 self.pydata['var_info_header'].append(k)
                 self.pydata['var_info'].append(list(map(fstr, item)))
        self.pydata['var_info'] = list(zip(*self.pydata['var_info']))
        # geno_info
        self.pydata['geno_info'] = []
        self.pydata['geno_info_header'] = []
        for k, item in list(geno_info.items()):
            self.pydata['geno_info_header'].append(k)
            if recode_missing:
                self.pydata['geno_info'].append([[missing_code if math.isnan(e) else e for e in x] for idx, x in enumerate(item) if which[idx]])
            else:
                self.pydata['geno_info'].append([x for idx, x in enumerate(item) if which[idx]])
        # convert geno_info to 3 dimensions:
        # D1: samples
        # D2: variants
        # D3: geno_info 
        self.pydata['geno_info'] = list(zip(*self.pydata['geno_info']))
        self.pydata['geno_info'] = [list(zip(*item)) for item in self.pydata['geno_info']]
        unique_names = self.sample_names
        if len(self.sample_names) != len(set(self.sample_names)):
            env.logger.warning("Duplicated sample names found. Using 'sample_ID.sample_name' as sample names") 
            unique_names = ["{0}.{1}".format(i,s) for i,s in zip(self.sample_IDs, self.sample_names)]
        self.pydata['sample_name'] = [str(x) for idx, x in enumerate(unique_names) if which[idx]]
        self.pydata['phenotype_name'] = self.phenotype_names
        self.pydata['phenotype'] = [x for idx, x in enumerate(self.phenotypes[0]) if which[idx]]
        if self.covariates:
            self.pydata['covariate_name'] = self.covariate_names
            # skip the first covariate, a vector of '1''s
            self.pydata['covariates'] = [[x for idx, x in enumerate(y) if which[idx]] for y in self.covariates[1:]]
        #
        if len(self.pydata['genotype']) == 0 or len(self.pydata['phenotype']) == 0 or len(self.pydata['genotype'][0]) == 0:
            raise ValueError("No input data")
        if len(self.pydata['geno_info']) > 0 and len(self.pydata['genotype']) != len(self.pydata['geno_info']):
            raise ValueError("Genotype and genotype information do not match")

    def getChr(self,variantID,cur):
        find_chr="SELECT chr from variant where variant_id={0}".format(variantID)
        chr= [rec[0] for rec in cur.execute(find_chr)]
        return chr[0]

    def getChrs(self,variantIDs,cur):
        idString="("+",".join([str(variantID) for variantID in variantIDs])+")"
        find_chr="SELECT chr,variant_id from variant where variant_id in "+idString
        varDict={}
        for rec in cur.execute(find_chr):
            if rec[0] not in varDict:
                varDict[rec[0]]=[]
            varDict[rec[0]].append(rec[1])
        return varDict

    def transformGeneName(self,geneSymbol):
        if ("-" in geneSymbol):
            geneSymbol=geneSymbol.replace("-","_")
        pattern=re.compile(r'\.')
        if pattern.findall(geneSymbol):
            geneSymbol=geneSymbol.replace(".","_")
        return geneSymbol



    def getGenotype_HDF5(self,group):
        """This function gets genotype of variants in specified group.
        """
        where_clause = ' AND '.join(['{0}={1}'.format(x, self.db.PH) for x in self.group_names])
        cur = self.db.cursor()
        # variant info
        var_info, variant_ids = self.getVarInfo(group, where_clause)
        chr=self.getChr(variant_ids[0],cur)
        chrEnd=self.getChr(variant_ids[-1],cur)

        if chr!=chrEnd:
            varDict=self.getChrs(variant_ids,cur)
            chrs=list(varDict.keys())
        else:
            varDict={chr:variant_ids}
            chrs=[chr]


        # get genotypes / genotype info
        genotype = []
        geno_info = {x:[] for x in self.geno_info}
        # getting samples locally from my own connection

        geneSymbol=self.transformGeneName(group[0])
        files=glob.glob(self.path+"/tmp*_genotypes_multi_genes.h5")
        # files=glob.glob("tmp*_genotypes_multi_genes.h5")

        files=sorted(files, key=lambda name: int(name.split("/")[-1].split("_")[1]))
 
        for fileName in files:   
            accessEngine=Engine_Access.choose_access_engine(fileName)

            if len(chrs)==1:
                # colnames=accessEngine.get_colnames(chr,geneSymbol)
                # snpdict=accessEngine.get_geno_info_by_group(geneSymbol,chr)

                # for chr in chrs:
                    
                colnames=accessEngine.get_colnames(chr)
                snpdict=accessEngine.get_geno_by_group(chr,geneSymbol)
                accessEngine.close()
                for ID in colnames:
                    data=snpdict[ID]
                    
                    gtmp = [data.get(x, [self.g_na] + [float('NaN')]*len(self.geno_info)) for x in varDict[chr]]
                    # handle -1 coding (double heterozygotes)     
                    genotype.append([2.0 if x[0] == -1.0 else x[0] for x in gtmp])

                    #
                    # handle genotype_info
                    #
                    for idx, key in enumerate(self.geno_info):
                        geno_info[key].append([x[idx+1] if (type(x[idx+1]) in [int, float]) else float('NaN') for x in gtmp])
            else:    
                colnames=accessEngine.get_colnames(chrs[0])
                    
                for ID in colnames:
                    gtmp=[]
                    for chr in chrs:
                        snpdict=accessEngine.get_geno_by_group(chr,geneSymbol)
                        data=snpdict[ID]
                        
                        gtmp.extend([data.get(x, [self.g_na] + [float('NaN')]*len(self.geno_info)) for x in varDict[chr]])
                        # handle -1 coding (double heterozygotes)     
                    genotype.append([2.0 if x[0] == -1.0 else x[0] for x in gtmp])
                    #
                    # handle genotype_info
                    #
                    for idx, key in enumerate(self.geno_info):
                        geno_info[key].append([x[idx+1] if (type(x[idx+1]) in [int, float]) else float('NaN') for x in gtmp])

                accessEngine.close()
        gname = ':'.join(list(map(str, group)))
        return self.filterGenotype(genotype, geno_info, var_info, gname)


    def run(self):
        
        grp = self.grp
        #
        try:
            grpname = ":".join(list(map(str, grp)))
        except TypeError:
            grpname = "None"
        # if grp is None:
        #     break
        # env.logger.debug('Retrieved association unit {}'.format(repr(grpname)))
        #
        #
        self.data = AssoData()
        self.pydata = {}
        values = list(grp)

        try:
            
            genotype, which, var_info, geno_info = self.getGenotype_HDF5(grp)
      

            # if I throw an exception here, the program completes in 5 minutes, indicating
            # the data collection part takes an insignificant part of the process.
            # 
            # set C++ data object
            if (len(self.tests) - self.num_extern_tests) > 0:
                self.setGenotype(which, genotype, geno_info, grpname)
                self.setPhenotype(which)
                self.setVarInfo(var_info)
            # set Python data object, for external tests
            if self.num_extern_tests:
                self.setPyData(which, genotype, var_info, geno_info, None, grpname)
            # association tests
            for test in self.tests:
                test.setData(self.data, self.pydata)
                result = test.calculate(env.association_timeout)
                # env.logger.debug('Finished association test on {}'.format(repr(grpname)))
                values.extend(result)
        # except KeyboardInterrupt as e:
        #     # die silently if stopped by Ctrl-C
        #     break
        except Exception as e:
            # env.logger.debug('An ERROR has occurred in process {} while processing {}: {}'.\
            #                   format(self.index, repr(grpname), e),exc_info=True)
            # self.data might have been messed up, create a new one
            print(e)
            self.data = AssoData()
            self.pydata = {}
            # return no result for any of the tests if an error message is captured.
            values.extend([float('NaN') for x in range(len(self.result_fields) - len(list(grp)))])
        print(values)

@app.task(bind=True,default_retry_delay=10,serializer="pickle")
def run_grp_association(self, param, grp,  result_fields,methods,path):
        worker=AssoTestsWorker(param, grp,  result_fields,methods,path)
        worker.run()

