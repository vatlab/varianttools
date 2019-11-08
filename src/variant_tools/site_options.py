#
# Site configuration file for site-wise installation of variant tools.
#

# ; separated pragmas for sqlite database that can be used to optimize the
# performance of sqlite database operations. Please check the sqlite manual
# for more details.
sqlite_pragma = ''

# Number of processes to read from input files during multi-process import.
# A smaller number might provide better performance for file systems with
# slow random I/O access.
import_num_of_readers = 2

# Number of processes to read genotype data for association tests. The default
# value is the minimum of value of option --jobs and 8. A smaller number
# might provide better performance for file system with slow random I/O
# access.
associate_num_of_readers = None

# root of tempory directory to store temporary files (default to system default)
# Setting it to a different physical disk than user projects can generally
# improve the performance of variant tools. It should also be set to a directory
# with at least 500G of free diskspace if the default temporary partition
# is small.
temp_dir = None

# A ;-separated list of URL that host the variant tools repository. This option
# should only be changed if you have created a local mirror of the variant tools
# repository. Adding the URL before the default URL might provide better
# downloading performance for your users. Removing the default URL is possible
# but not recommended.
search_path = 'http://bioinformatics.mdanderson.org/Software/VariantTools/repository/;http://bioinformatics.mdanderson.org/Software/VariantTools/archive/'

# A directory for shared resource files. It can be configured as
# 1. NO shared resource (default). All users maintain their own resource
#   directory ($local_resource, which is usually ~/.variant_tools).
#
# 2. A read-only directory with a mirror of the variant tools repository, with
#   .DB.gz files decompressed in the annoDB directory. This is important because
#   otherwise each user will have to decompress the files in their local resource
#   directory. The system admin can choose to remove outdated databases to reduce
#   the use of disk space. This option requires regular update of the resources.
#
# 3. A directory that is writable by all users. The resources will be downloaded
#   to this directory by users, and shared by all users. This option is easier
#   to implement and requires less maintenance. The system admin can choose to
#   mirror the variant tools repository and let the users to keep it up to date.
shared_resource = None
