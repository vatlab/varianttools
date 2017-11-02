#!/usr/bin/env python3
import os
import psutil
import threading
import time
from datetime import datetime
import stat
import subprocess
import sys



class ProcessMonitor(threading.Thread):
    def __init__(self, task_id, monitor_interval,resource_monitor_interval, command):
        threading.Thread.__init__(self)
        self.task_id = task_id
        self.pid = os.getpid()
        self.monitor_interval = monitor_interval
        self.resource_monitor_interval=resource_monitor_interval
        self.daemon = True
        self.command=command
        # self.pulse_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.pulse')
        #self.pulse_file = task_id + '.pulse'
        # How to know whether we are import hdf5 or sqlite when we are running the import command? The pulse file should be further than the h.5 files
        if 'HDF' in self.command or 'hdf' in self.command:
            self.pulse_file = command.split(' ')[1] + '_j' + command.split(' ')[[i for i, x in enumerate(command.split(' ')) if '-j' == x][0]+1] + '_hdf5' + '.pulse'
        else:
            self.pulse_file = command.split(' ')[1] + '_j' + command.split(' ')[[i for i, x in enumerate(command.split(' ')) if '-j' == x][0]+1] + '_sqlite' + '.pulse'
        # remove previous status file, which could be readonly if the job is killed
        if os.path.isfile(self.pulse_file):
            if not os.access(self.pulse_file, os.W_OK):
                os.chmod(self.pulse_file, stat.S_IREAD | stat.S_IWRITE)
            os.remove(self.pulse_file)
        with open(self.pulse_file, 'w') as pd:
            pd.write('#task: {}\n'.format(task_id))
            pd.write('#started at {}\n#\n'.format(datetime.now().strftime("%A, %d. %B %Y %I:%M%p")))
            pd.write('#time\tproc_cpu\tproc_mem\tchildren\tchildren_cpu\tchildren_mem\tdisk_use\n')

    def _check(self):
        current_process = psutil.Process(self.pid)
        par_cpu = current_process.cpu_percent()
        par_mem = current_process.memory_info()[0]
        ch_cpu = 0
        ch_mem = 0
        children = current_process.children(recursive=True)
        n_children = len(children)
        #disk_use = psutil.disk_usage('/Users/manchongleong/Documents/testVariantTools')
        #disk_use = os.path.getsize('/Users/manchongleong/Documents/testVariantTools/')
        for child in children:
            ch_cpu += child.cpu_percent()
            ch_mem += child.memory_info()[0]
        return par_cpu, par_mem, n_children, ch_cpu, ch_mem

    def size(self):
        disk_use = 0
        for root, dirs, files in os.walk(os.getcwd(), topdown=False):
            disk_use += sum(os.path.getsize(os.path.join(root,file)) for file in files if os.path.splitext(file)[1] in ('.h5', '.DB'))
            dirs = [x for x in dirs if not x.startswith('.')]
        return disk_use

    def run(self):
        counter = 0
        start_time = time.time()

        while True:
            try:
                if not os.access(self.pulse_file, os.W_OK):
                    # the job should be killed
                    p = psutil.Process(self.pid)
                    p.kill()
                # most of the time we only update
                # if counter % self.resource_monitor_interval:
                #     os.utime(self.pulse_file, None)
                # else:
                cpu, mem, nch, ch_cpu, ch_mem = self._check()
                disk_use = self.size()
                with open(self.pulse_file, 'a') as pd:
                    pd.write('{}\t{:.2f}\t{}\t{}\t{}\t{}\t{}\n'.format(time.time()-start_time, cpu, mem, nch, ch_cpu, ch_mem, disk_use))
                time.sleep(self.monitor_interval)
                counter += 1
            except Exception as e:
                # if the process died, exit the thread
                # the warning message is usually:
                # WARNING: psutil.NoSuchProcess no process found with pid XXXXX
                #env.logger.warning(e)
                print('Monitor of {} failed with message {}'.format(self.task_id, e))
                # env.logger.debug('Monitor of {} failed with message {}'.format(self.task_id, e))
                break

def main():
    monitor_interval = 2
    resource_monitor_interval = 60
    task_id=datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    command=sys.argv[1]
    print(command)
    m = ProcessMonitor(task_id, monitor_interval=monitor_interval,
             resource_monitor_interval=resource_monitor_interval, command=sys.argv[1])
    m.start()
    subprocess.call(command,shell=True)

if __name__ == "__main__":
    main()
