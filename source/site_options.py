# 
# Site configuration files for site-wise installation of variant tools.
# 

# ; separated pragmas for sqlite database that can be used to optimize the
# performance of sqlite database operations. Please check the sqlite manual
# for more details.
sqlite_pragma=''

# Number of processes to read from input files during multi-process import.
# A smaller number might provide better performance for file systems with
# slow random I/O access.
import_num_of_readers=2

# Number of processes to read genotype data for association tests. The default
# value is the minimum of value of option --jobs and 8. A smaller number
# might provide better performance for file system with slow random I/O
# access.
associate_num_of_readers=None

# root of tempory directory to store temporary files (default to system default)
# Setting it to a different physical disk than user projects can generally 
# improve the performance of variant tools. It should also be set to a directory
# with at least 500G of free diskspace if the default temporary partition
# is small.
temp_dir=None

# A ;-separated list of URL that host the variant tools repository. This option
# should only be changed if you have created a local mirror of the variant tools
# repository. Adding the URL before the default URL might provide better 
# downloading performance for your users. Removing the default URL is possible
# but not recommended.
search_path='http://bioinformatics.mdanderson.org/Software/VariantTools/repository/'

# Default local_resource directory for users. Users by default will have their
# own copies of downloaded resources but a system admin can set a shared directory
# so that users can share those files. User-specific directory will be used
# if the system wide directory is not writable.
local_resource='~/.variant_tools'


