#This file is part of PyTransport.

#PyTransport is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#PyTransport is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with PyTransport.  If not, see <http://www.gnu.org/licenses/>.


# This file contains python scripts used to setup the complided PyTrans module

import sympy as sym
import numpy as np
import math
import subprocess
import sys
import os
import shutil



def directory():
    dir = os.path.dirname(__file__)
    filename = os.path.join(dir, 'PyTrans/PyTrans.cpp')
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    f = open(filename,"w")
    for line in lines:
        if not  line.endswith("//evolve\n") and not line.endswith("//moments\n") and not line.endswith("//model\n") and not line.endswith("//stepper\n"):
            f.write(line)
        if line.endswith("//evolve\n"):
            fileT = os.path.join(dir, 'CppTrans/evolve.h')
            f.write('#include' + '"'+ fileT +'"' + '//evolve' +'\n')
        if line.endswith("//moments\n"):
            fileT = os.path.join(dir, 'CppTrans/moments.h')
            f.write('#include' + '"'+ fileT +'"' + '//moments' +'\n')
        if line.endswith("//model\n"):
            fileT = os.path.join(dir, 'CppTrans/model.h')
            f.write('#include' + '"'+ fileT +'"' + '//model' +'\n')
        if line.endswith("//stepper\n"):
            fileT = os.path.join(dir, 'CppTrans/stepper/rkf45.hpp')
            f.write('#include' + '"'+ fileT +'"' + '//stepper' +'\n')

    f.close()

def pathSet():
    dir = os.path.dirname(__file__)
    path1 = os.path.join(dir, 'PyTrans/lib/python/')
    path2 = os.path.join(dir, 'PyTransScripts/')
    sys.path.append(path1)
    sys.path.append(path2)
 

def compile():   
    directory()
    name = ""
    dir = os.path.dirname(__file__)
    location = os.path.join(dir, 'PyTrans/')
    filename1 = os.path.join(dir, 'PyTrans/moduleSetup.py')
    f = open(filename1,"r")
    lines = f.readlines()
    f.close()
    f = open(filename1,"w")
    for line in lines:
        if not  line.endswith("#setup"):
            f.write(line)
        if line.endswith("#setup"):
            f.write('setup(name="PyTrans'+name+'", version="1.0", ext_modules=[Extension("PyTrans'+name+'", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs])#setup')
    f.close()    
    
    filename = os.path.join(dir, 'PyTrans/PyTrans.cpp')
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    f = open(filename,"w")
    for line in lines:
        if not  line.endswith("//FuncDef\n") and not line.endswith("//initFunc\n"):
            f.write(line)
        if line.endswith("//FuncDef\n"):
            f.write('static PyMethodDef PyTrans'+name+'_funcs[] = {{"H", (PyCFunction)  _H,    METH_VARARGS, PyTrans_docs},{"nF", (PyCFunction)MT_fieldNumber,        METH_VARARGS, PyTrans_docs},{"nP", (PyCFunction)MT_paramNumber,        METH_VARARGS, PyTrans_docs},{"V", (PyCFunction)MT_V,            METH_VARARGS, PyTrans_docs},{"dV", (PyCFunction)MT_dV,                METH_VARARGS, PyTrans_docs},  {"ddV", (PyCFunction)MT_ddV,                METH_VARARGS, PyTrans_docs},  {"backEvolve", (PyCFunction)MT_backEvolve,        METH_VARARGS, PyTrans_docs},    {"sigEvolve", (PyCFunction)MT_sigEvolve,        METH_VARARGS, PyTrans_docs},    {"alphaEvolve", (PyCFunction)MT_alphaEvolve,        METH_VARARGS, PyTrans_docs},    {NULL}};//FuncDef\n')
        if line.endswith("//initFunc\n"):
            f.write('void initPyTrans'+name+'(void)    {        Py_InitModule3("PyTrans'+name+'", PyTrans'+name+'_funcs,                       "Extension module for inflationary statistics");        import_array();   }//initFunc\n')
    f.close()    
    
    subprocess.call(["python", "/PyTrans/moduleSetup.py", "install", "--home="+location],cwd=location)
    sys.path.append(location+"/lib/python/")
    sys.path.append(location+"../PyTransScripts")
    shutil.rmtree(location+"/build/")



def compileName(name):    
    directory()
    dir = os.path.dirname(__file__)
    location = os.path.join(dir, 'PyTrans/')
    filename1 = os.path.join(dir, 'PyTrans/moduleSetup.py')
    f = open(filename1,"r")
    lines = f.readlines()
    f.close()
    f = open(filename1,"w")
    for line in lines:
        if not  line.endswith("#setup"):
            f.write(line)
        if line.endswith("#setup"):
            f.write('setup(name="PyTrans'+name+'", version="1.0", ext_modules=[Extension("PyTrans'+name+'", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs])#setup')
    f.close()

    filename = os.path.join(dir, 'PyTrans/PyTrans.cpp')
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    f = open(filename,"w")
    for line in lines:
        if not  line.endswith("//FuncDef\n") and not line.endswith("//initFunc\n"):
            f.write(line)
        if line.endswith("//FuncDef\n"):
            f.write('static PyMethodDef PyTrans'+name+'_funcs[] = {{"H", (PyCFunction)MT_H,    METH_VARARGS, PyTrans_docs},{"nF", (PyCFunction)MT_fieldNumber,        METH_VARARGS, PyTrans_docs},{"nP", (PyCFunction)MT_paramNumber,        METH_VARARGS, PyTrans_docs},{"V", (PyCFunction)MT_V,            METH_VARARGS, PyTrans_docs},{"dV", (PyCFunction)MT_dV,                METH_VARARGS, PyTrans_docs},  {"ddV", (PyCFunction)MT_ddV,                METH_VARARGS, PyTrans_docs},  {"backEvolve", (PyCFunction)MT_backEvolve,        METH_VARARGS, PyTrans_docs},    {"sigEvolve", (PyCFunction)MT_sigEvolve,        METH_VARARGS, PyTrans_docs},    {"alphaEvolve", (PyCFunction)MT_alphaEvolve,        METH_VARARGS, PyTrans_docs},    {NULL}};//FuncDef\n')
        if line.endswith("//initFunc\n"):
            f.write('void initPyTrans'+name+'(void)    {        Py_InitModule3("PyTrans'+name+'", PyTrans'+name+'_funcs,                       "Extension module for inflationary statistics");        import_array();   }//initFunc\n')
    f.close()

    subprocess.call(["python", filename1, "install", "--home=" + location],cwd=location)
    sys.path.append(location+"/lib/python/")
    sys.path.append(location+"../PyTransScripts")
    shutil.rmtree(location+"/build/")



def deleteModule(name):
    location = os.path.join(dir, 'PyTrans/')
    [os.remove(os.path.join(location,f)) for f in os.listdir(location) if f.startswith("PyTrans"+name)]

#os.remove(location+"/lib/python/PyTrans"+name+".so")
#os.remove(location+"/lib/python/PyTrans"+name+"-1.0-py2.7.egg-info")


    
def tol(rtol, atol):
    dir = os.path.dirname(__file__)
    filename = os.path.join(dir, 'PyTrans/PyTrans.cpp')
    f = open(filename,"r")  

    lines = f.readlines()
    f.close()
   
    f = open(filename,"w")  
    for line in lines:
        if not  line.endswith("//tols\n"):
            f.write(line)
        if line.endswith("//tols\n"):
            f.write('    double rtol='+str(rtol)+', atol='+str(atol)+';//tols\n')
    f.close()



def potential(V,nF,nP):
    f=sym.symarray('f',nF)
    p=sym.symarray('p',nP)

    vd=sym.symarray('vd',nF)
    vdd=sym.symarray('vdd',nF*nF)
    vddd=sym.symarray('vddd',nF*nF*nF)

    for i in range(nF):
        vd[i] = sym.simplify(V.diff(f[i]))
        for j in range(nF):
            vdd[i+j*nF] = sym.simplify(V.diff(f[i]).diff(f[j]) )
                        #sym.simplify(V.diff(f[i]).diff(f[j]))
            for k in range(nF):
                vddd[i+j*nF+k*nF*nF] = sym.simplify(V.diff(f[i]).diff(f[j]).diff(f[k]))
                    #sym.simplify(V.diff(f[i]).diff(f[j]).diff(f[k]))


    import os
    dir = os.path.dirname(__file__)
    filename1 = os.path.join(dir, 'CppTrans/potentialProto.h')
    filename2 = os.path.join(dir, 'CppTrans/potential.h')
    f = open(filename1, 'r')
    g = open(filename2, 'w')


    for line in f:
        g.write(line)
        if line == "// #FP\n":
            g.write('nF='+str(nF)+';\n'+'nP='+str(nP)+';\n')
        #if line == "// pdef\n":
            #for i in range(nP):
            #    g.write('p['+str(i)+']='+str(pp[i])+';\n')
        if line == "// Pot\n":
            expr=str(sym.ccode(V))
            for l in range(max(nP,nF)): 
                expr=expr.replace("_"+str(l),"["+str(l)+"]")
            g.write('\n sum='+str(expr)+';\n')
        if line == "// dPot\n":
            for i in range(nF):
                expr=str(sym.ccode(vd[i]))
                for l in range(max(nF,nP)): 
                    expr=expr.replace("_"+str(l),"["+str(l)+"]")
                g.write('\n sum['+str(i)+']='+str(expr)+';\n')    

        if line == "// ddPot\n":
              for i in range(nF):
                  for j in range(nF):
                      expr=str(sym.ccode(vdd[i+nF*j]))
                      for l in range(max(nF,nP)): 
                          expr=expr.replace("_"+str(l),"["+str(l)+"]")
                      g.write('\n sum['+str(i)+'+nF*'+str(j)+']='+str(expr)+';\n')    
        if line == "// dddPot\n":
              for i in range(nF):
                  for j in range(nF):
                      for k in range(nF):
                          expr=str(sym.ccode(vddd[i+nF*j+nF*nF*k]))
                          for l in range(max(nF,nP)): 
                              expr=expr.replace("_"+str(l),"["+str(l)+"]")
                          g.write('\n sum['+str(i)+'+nF*'+str(j)+'+nF*nF*'+str(k)+']='+str(expr)+';\n')    
                          
       
    g.close()
    f.close()
    