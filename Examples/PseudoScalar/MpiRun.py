import subprocess


#subprocess.Popen(['/usr/local/bin/mpiexec', '-n', '10', 'python', 'MpiEqBi.py']) 
subprocess.Popen(['mpiexec', '-n', '10', 'python', 'MpiAlpBetBi.py']) 


