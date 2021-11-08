from .tools import *
from math import ceil

def run_llr(command, dir, mode, options, parameters):
    # Create script
    with open(parameters['exec_script'], 'w') as exec_script_desc:
        NODES=int(ceil(options['mpi']/2.))
        options['ppn'] = {"jollyjumper":24, "tornado":36}[options['partition']]
        exec_script_desc.write(
            "#PBS -l nodes="+str(NODES)+":ppn="+str(options['ppn'])+" \n"
            +"#PBS -q default \n"
            +"#PBS -j oe\n"
            +"#PBS -l walltime="+options['max_time']+"\n"
            +"module purge\n"
            +"unset MODULEPATH;\n"
            +"module use /opt/exp_soft/vo.llr.in2p3.fr/modulefiles_el7\n"
            +"module load hdf5/1.10.5-icc-omp4.1.1\n"
            +"module load python/3.7.0\n"
            +"module load h5py/hdf5_1.10.5-icc-omp4.1.1-py3.7.0\n"
            +"module load mpi4py/omp4.1.1-ib-icc_py3.7.0\n"
            +"module load compilers/gcc/9.x.x\n"
            +"module load fftw/3.3.10-omp-4.1.1-icc-19 \n"
            +"export OMP_NUM_THREADS="+str(options['omp'])+" \n"
            +"export OMP_SCHEDULE=DYNAMIC \n"
            +"export KMP_AFFINITY=verbose \n"
            +"export PATH=$PATH:/opt/exp_soft/vo.llr.in2p3.fr/GALOP/beck \n"
            +"export LIBPXR=/home/llr/galop/derouil/applications.ompi216.Py3/picsar/lib \n"
            +"export LD_LIBRARY_PATH=$LIBPXR:$LD_LIBRARY_PATH \n"
            +"export FFTW3_LIB=/opt/exp_soft/vo.llr.in2p3.fr/fftw/3.3.10/opm-4.1.1-intel-19-el7/lib\n"
            +"export FFTW3_INC=/opt/exp_soft/vo.llr.in2p3.fr/fftw/3.3.10/opm-4.1.1-intel-19-el7/include\n"
            +"ulimit -s unlimited \n"
            +"#Specify the number of sockets per node in -mca orte_num_sockets \n"
            +"#Specify the number of cores per sockets in -mca orte_num_cores \n"
            +"cd "+dir+" \n"
            +"module list 2> module.log\n"
            +command+" \n"
            +"echo $? > exit_status_file \n"
        )
    if options['partition'] == "jollyjumper":
        JOB = "PBS_DEFAULT=llrlsi-jj.in2p3.fr qsub  "+parameters['exec_script']
    elif options['partition'] == "tornado":
        JOB = "PBS_DEFAULT=poltrnd.in2p3.fr qsub  "+parameters['exec_script']
    launch_job(command, JOB, dir, options['max_time_seconds'], parameters['output_file'], repeat=2, verbose=options['verbose'])
