# Example of usage:
# pyhton scheduler.py -p /mnt/avoton/fog/data/dataset4test/split16 -d Ecoli_2D -n 8

from multiprocessing import Pool, Queue
import argparse
import sys
import subprocess
import operator
import os
import ntpath
import random
import shutil

# Define usage and takes input parameters
def msg_usage(name=None):
    return ''' python scheduler.py -path INPUT_PATH -dataset DATASET_NAME -n MAX_NUM_CORE'''

def initialize():
    parser = argparse.ArgumentParser(description="Scheduler for bioinfo analys on SoCs",usage=msg_usage())
    parser.add_argument('-p', '--path', type=str, default=os.getcwd(), help="Path where the folders of the splitted DB reside. Default: current directory")
    parser.add_argument('-d', '--dataset', type=str, default='Ecoli_2D',  help="Dataset. Default: Ecoli_2D")
    parser.add_argument('-n', '--max_num_core', type=int, default=1, help="Max number of cores available in the platform. Default: 1")


    args = parser.parse_args()
    path=args.path
    dataset=args.dataset
    max_num_core=args.max_num_core

    return path, dataset, max_num_core

# Process each folder:
#	- launch Deepnano
#	- launch Kraken
#	- append labels to a copy of aws_scheme.js to be uploaded to AWS

def process_folder(queue):
	#print os.getpid(),"working" 
	while not queue.empty():
		folder = queue.get()
	        print os.getpid(), "got", folder

		# Process folder
		short_dir=os.path.basename(os.path.normpath(folder))
		# Launch Deepnano for all files in folder
		# Safer to use THEANO FLAGS
		num_random=str(random.randint(0,32767))
		# Create dir for outputs if it does not exist
		if not os.path.exists('outputs'):
			os.makedirs('outputs')
		print os.getpid(), " performing base-calling with Deepnano"
		cmd_deepnano="OMP_NUM_THREADS=1  THEANO_FLAGS=base_compiledir=/tmp/" + num_random + "/theano.NOBACK python  basecall_no_metrichor.py  --directory "+ folder + " --output outputs/" + short_dir +".fasta 1>deepnano_exectime_" + short_dir + ".txt"
                ip2 = subprocess.check_call(cmd_deepnano, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

		# Upload Minikraken DB to RAM
		#

		# Launch Kraken, with Deepnano output file as input.
		# Kraken prints output as standard error. We save it to file_seq
                file_seq = "outputs/sequences_" + short_dir +".kraken"
		print os.getpid(), " performing bacteria identification with (mini)Kraken"
                cmd_kraken="../kraken-0.10.5-beta//kraken  --db ../minikraken_20141208/  outputs/" + short_dir +".fasta>" + file_seq + " 2>kraken_perf_"+ short_dir +".txt" 
                ip2 = subprocess.check_call(cmd_kraken, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)

		# Append sequences outputted by Kraken to a copy of file aws_scheme.js 
                print os.getpid(), " appending labels to file aws_label.js"
		with open(file_seq, "r") as seq, open("aws_label.js", "a") as aws_file:
                       for line in seq:
                                # read sequence line according to Kraken output format
                                splitted = line.split()
                                kind = splitted[0]
                                fasta_id = splitted[1]
                                kraken_id = splitted[2]
                                bp = splitted[3]
                                LCA_mapping = ""
                                for index in range(4,len(splitted)): # everything that follows
                                        LCA_mapping += splitted[index] + " "
                                LCA_mapping = LCA_mapping.rstrip()
				label = "device.publish('air_filter/nanopore_2', JSON.stringify({\"kind\": \"" +  kind + "\", \"fasta_id\": \"" + fasta_id + "\", \"kraken_id\": " + kraken_id + " , \"bp\": " + bp + ", \"LCA_mapping\": \"" + LCA_mapping + "\" }));"
				aws_file.write(label+ "\n")

	os.getpid(), " done. Bye!"

# do multiprocessing for folders in queue
if __name__ == '__main__':
	
	# Read input parameters
        parameters=[0.]*3
	parameters=initialize()

        path = parameters[0]
        dataset = parameters[1]
        n_proc=int(parameters[2])

	# Labels will be appended by each process to a file which is a copy of aws_scheme.js
        shutil.copyfile('aws_scheme.js','aws_label.js')

	print "Starting job scheduler to run ", n_proc, " jobs in folders ", path, " for the dataset ", dataset 
	queue = Queue()

        # Put folders into queue
        list_dir = [dir_name for dir_name in os.listdir(path) if dir_name.startswith(dataset)]
        for folder in list_dir:
                queue.put(path+"/"+folder)

        pool = Pool(n_proc, process_folder, (queue,))

	pool.close() # signal that no more tasks will be submitted to pool
	pool.join() # wait until all processes are done

	# Complete upload to AWS by closing the device
	with open("aws_label.js", "a") as aws_file:
		aws_file.write("});")
	print "All processes have done. File aws_label.js is ready to fly."
	#remove content of /outputs
	




