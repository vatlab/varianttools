import argparse

def checking_association(sql_file, hdf5_file):
	with open(str(sql_file)) as f1:
		result_sql = {}
		for line in f1:
			fields = line.split()
			if len(fields) == 7:
				result_sql[fields[0]] = fields[1:7]

	with open(str(hdf5_file)) as f2:
		result_h5 = {}
		for line in f2:
			fields = line.split()
			if len(fields) == 7:
				result_h5[fields[0]] = fields[1:7]

	count=0

	for key in result_h5 :
		if result_h5[key] != result_sql[key]:
			count+=1
			print("Values in %s are different:\n"  % (key),"          %s \n" % (result_sql["refgene_name2"]) ,"  in sql: %s \n  in hdf5: %s" % (result_sql[key], result_h5[key]))

	if count!=0:
		print("there are total %d different results." % count)
	else:
		print("All result are the same!")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="the two files path and name, for sqlite and -for hdf5")
    parser.add_argument("-sql",
                        help="for sqlite and -for sqlite")
    parser.add_argument("-h5",
                        help="for hdf5 and -for hdf5")

    args = parser.parse_args()
    checking_association(args.sql, args.h5)
