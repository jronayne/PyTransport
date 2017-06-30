import subprocess

### only need to do this once for each new potential -- it creates a potential.h file for the potential, and then compliles a python extension model InfEasyPy
#subprocess.call(["python", "potentialSetup.py"]) 

#subprocess.check_output(["ls", "-l"])
subprocess.Popen(['/usr/local/bin/mpiexec', '-n', '20', 'python', 'alpbetMpi.py']) 

