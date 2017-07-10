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
import sys
import site
import subprocess
import platform
import os
import shutil
import time
from gravipy import *


def directory(NC):
    dir = os.path.dirname(__file__)
    filename = os.path.join(dir, 'PyTrans', 'PyTrans.cpp')
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    f = open(filename,"w")
    if NC==False:
        for line in lines:
            if not  line.endswith("//evolve\n") and not line.endswith("//moments\n") and not line.endswith("//model\n") and not line.endswith("//stepper\n"):
                f.write(line)
            if line.endswith("//evolve\n"):
                fileT = os.path.join(dir, 'CppTrans', 'evolve.h')
                f.write('#include' + '"'+ fileT +'"' + '//evolve' +'\n')
            if line.endswith("//moments\n"):
                fileT = os.path.join(dir, 'CppTrans', 'moments.h')
                f.write('#include' + '"'+ fileT +'"' + '//moments' +'\n')
            if line.endswith("//model\n"):
                fileT = os.path.join(dir, 'CppTrans', 'model.h')
                f.write('#include' + '"'+ fileT +'"' + '//model' +'\n')
            if line.endswith("//stepper\n"):
                fileT = os.path.join(dir, 'CppTrans', 'stepper', 'rkf45.hpp')
                f.write('#include' + '"'+ fileT +'"' + '//stepper' +'\n')
    else:
        for line in lines:
            if not  line.endswith("//evolve\n") and not line.endswith("//moments\n") and not line.endswith("//model\n") and not line.endswith("//stepper\n"):
                f.write(line)
            if line.endswith("//evolve\n"):
                fileT = os.path.join(dir, 'CppTrans', 'NC', 'evolve.h')
                f.write('#include' + '"'+ fileT +'"' + '//evolve' +'\n')
            if line.endswith("//moments\n"):
                fileT = os.path.join(dir, 'CppTrans', 'NC', 'moments.h')
                f.write('#include' + '"'+ fileT +'"' + '//moments' +'\n')
            if line.endswith("//model\n"):
                fileT = os.path.join(dir, 'CppTrans', 'NC', 'model.h')
                f.write('#include' + '"'+ fileT +'"' + '//model' +'\n')
            if line.endswith("//stepper\n"):
                fileT = os.path.join(dir, 'CppTrans', 'stepper', 'rkf45.hpp')
                f.write('#include' + '"'+ fileT +'"' + '//stepper' +'\n')
    f.close()

def pathSet():
    dir = os.path.dirname(__file__)
    site.addsitedir(dir)

    p = platform.system()
    if p is 'Windows':
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python', 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python', 'Lib', 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python' + sys.version[:3].translate(None, '.'), 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python' + sys.version[:3].translate(None, '.'), 'Lib', 'site-packages'))
    else:
        site.addsitedir(os.path.join(dir, 'PyTrans', 'lib', 'python', 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'lib', 'python' + sys.version[:3], 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'lib', 'site-python'))

    site.addsitedir(os.path.join(dir, 'PyTransScripts'))


def compileName(name,NC=False):    
    directory(NC)
    dir = os.path.dirname(__file__)
    location = os.path.join(dir, 'PyTrans')
    filename1 = os.path.join(dir, 'PyTrans', 'moduleSetup.py')
    f = open(filename1,"r")
    lines = f.readlines()
    f.close()
    f = open(filename1,"w")
    for line in lines:
        if not  line.endswith("#setup\n"):
            f.write(line)
        if line.endswith("#setup\n"):
            f.write('setup(name="PyTrans'+name+'", version="1.0", ext_modules=[Extension("PyTrans'+name+'", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs])#setup\n')
    f.close()
    filename = os.path.join(dir, 'PyTrans', 'PyTrans.cpp')
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    f = open(filename,"w")
    for line in lines:
        if not  line.endswith("//FuncDef\n") and not line.endswith("//initFunc\n") and not line.endswith("//modDef\n"):
            f.write(line)
        if line.endswith("//FuncDef\n"):
            f.write('static PyMethodDef PyTrans'+name+'_funcs[] = {{"H", (PyCFunction)MT_H,    METH_VARARGS, PyTrans_docs},{"nF", (PyCFunction)MT_fieldNumber,        METH_VARARGS, PyTrans_docs},{"nP", (PyCFunction)MT_paramNumber,        METH_VARARGS, PyTrans_docs},{"V", (PyCFunction)MT_V,            METH_VARARGS, PyTrans_docs},{"dV", (PyCFunction)MT_dV,                METH_VARARGS, PyTrans_docs},  {"ddV", (PyCFunction)MT_ddV,                METH_VARARGS, PyTrans_docs},  {"backEvolve", (PyCFunction)MT_backEvolve,        METH_VARARGS, PyTrans_docs},    {"sigEvolve", (PyCFunction)MT_sigEvolve,        METH_VARARGS, PyTrans_docs},    {"alphaEvolve", (PyCFunction)MT_alphaEvolve,        METH_VARARGS, PyTrans_docs},    {NULL}};//FuncDef\n')
        if line.endswith("//modDef\n"):
            f.write('     //modDef\n')
        if line.endswith("//initFunc\n"):
            f.write('void initPyTrans'+name+'(void)    {        Py_InitModule3("PyTrans'+name+'", PyTrans'+name+'_funcs,                       "Extension module for inflationary statistics");        import_array();   }//initFunc\n')
        
    f.close()

    my_env = os.environ.copy()
    my_env["PYTHONUSERBASE"] = location
    p = subprocess.Popen(["python", filename1, "install", "--user"], cwd=location, stdin=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)
    stdout, stderr = p.communicate()

    p = platform.system()
    if p is 'Windows':
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python', 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python', 'Lib', 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python' + sys.version[:3].translate(None, '.'), 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python' + sys.version[:3].translate(None, '.'), 'Lib', 'site-packages'))
    else:
        site.addsitedir(os.path.join(dir, 'PyTrans', 'lib', 'python', 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'lib', 'python' + sys.version[:3], 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'lib', 'site-python'))

    site.addsitedir(os.path.join(location, '..', 'PyTransScripts'))

    shutil.rmtree(os.path.join(location, 'build'), ignore_errors=True)

def compileName3(name,NC=False):
    directory(NC)
    dir = os.path.dirname(__file__)
    location = os.path.join(dir, 'PyTrans')
    filename1 = os.path.join(dir, 'PyTrans', 'moduleSetup.py')
    f = open(filename1,"r")
    lines = f.readlines()
    f.close()
    f = open(filename1,"w")
    for line in lines:
        if not  line.endswith("#setup\n"):
            f.write(line)
        if line.endswith("#setup\n"):
            f.write('setup(name="PyTrans'+name+'", version="1.0", ext_modules=[Extension("PyTrans'+name+'", [filename, filename2 ])], include_dirs=[numpy.get_include(), dirs])#setup\n')
    f.close()

    filename = os.path.join(dir, 'PyTrans', 'PyTrans.cpp')
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    f = open(filename,"w")
    for line in lines:
        if not  line.endswith("//FuncDef\n") and not line.endswith("//initFunc\n") and not line.endswith("//modDef\n"):
            f.write(line)
        if line.endswith("//FuncDef\n"):
            f.write('static PyMethodDef PyTrans'+name+'_methods[] = {{"H", (PyCFunction)MT_H,    METH_VARARGS, PyTrans_docs},{"nF", (PyCFunction)MT_fieldNumber,        METH_VARARGS, PyTrans_docs},{"nP", (PyCFunction)MT_paramNumber,        METH_VARARGS, PyTrans_docs},{"V", (PyCFunction)MT_V,            METH_VARARGS, PyTrans_docs},{"dV", (PyCFunction)MT_dV,                METH_VARARGS, PyTrans_docs},  {"ddV", (PyCFunction)MT_ddV,                METH_VARARGS, PyTrans_docs},  {"backEvolve", (PyCFunction)MT_backEvolve,        METH_VARARGS, PyTrans_docs},    {"sigEvolve", (PyCFunction)MT_sigEvolve,        METH_VARARGS, PyTrans_docs},    {"alphaEvolve", (PyCFunction)MT_alphaEvolve,        METH_VARARGS, PyTrans_docs},   {NULL, NULL, 0, NULL}};//FuncDef\n')

        if line.endswith("//modDef\n"):
            f.write('static struct PyModuleDef PyTransModule = {PyModuleDef_HEAD_INIT, "PyTrans'+name+'", PyTrans_docs, -1, PyTrans'+name+'_methods}; //modDef\n')

        if line.endswith("//initFunc\n"):
            f.write('PyMODINIT_FUNC PyInit_PyTrans'+name+'(void)    {    PyObject *m = PyModule_Create(&PyTransModule); import_array(); return m;} //initFunc\n')
    f.close()

    my_env = os.environ.copy()
    my_env["PYTHONUSERBASE"] = location
    p = subprocess.Popen(["python", filename1, "install", "--user"], cwd=location, stdin=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)
    stdout, stderr = p.communicate()

    p = platform.system()
    if p is 'Windows':
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python', 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python', 'Lib', 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python' + sys.version[:3].translate(None, '.'), 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'Python' + sys.version[:3].translate(None, '.'), 'Lib', 'site-packages'))
    else:
        site.addsitedir(os.path.join(dir, 'PyTrans', 'lib', 'python', 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'lib', 'python' + sys.version[:3], 'site-packages'))
        site.addsitedir(os.path.join(dir, 'PyTrans', 'lib', 'site-python'))

    site.addsitedir(os.path.join(location, '..', 'PyTransScripts'))

    shutil.rmtree(os.path.join(location, 'build'), ignore_errors=True)



def deleteModule(name):
    location = os.path.join(dir, 'PyTrans')
    [os.remove(os.path.join(location,f)) for f in os.listdir(location) if f.startswith("PyTrans"+name)]

#os.remove(location+"/lib/python/PyTrans"+name+".so")
#os.remove(location+"/lib/python/PyTrans"+name+"-1.0-py2.7.egg-info")


    
def tol(rtol, atol):
    dir = os.path.dirname(__file__)
    filename = os.path.join(dir, 'PyTrans', 'PyTrans.cpp')
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



def potential(V,nF,nP,G=0,simple=False):
    f=sym.symarray('f',nF)
    p=sym.symarray('p',nP)

    vd=sym.symarray('vd',nF)
    vdd=sym.symarray('vdd',nF*nF)
    vddd=sym.symarray('vddd',nF*nF*nF)


    if G!=0:
        g, Ga, Ri, Rm =fieldmetric(G,nF,nP,simple)
        FMP=0
        for i in range(nF):
            if simple==True:
				vd[i] = sym.simplify(V.diff(f[i]))
            else:
                vd[i] = V.diff(f[i])
        for i in range(nF):
            for j in range(nF):
                for l in range(nF):
                    FMP=FMP+Ga(-(l+1),i+1,j+1) * vd[l]
                if simple==True:	
                    vdd[i+j*nF] = sym.simplify(V.diff(f[i]).diff(f[j])-FMP)
                else:
                    vdd[i+j*nF] = V.diff(f[i]).diff(f[j])-FMP
                FMP=0
        for i in range(nF):
            for j in range(nF):
                for k in range(nF):
                    for l in range(nF):
                       FMP=FMP+Ga(-(l+1),i+1,k+1)*vdd[l+j*nF] + Ga(-(l+1),j+1,k+1)*vdd[i+l*nF] +sym.expand(Ga(-(1+l),1+i,1+j)).diff(f[k])*vd[l]+sym.expand(Ga(-(1+l),1+i,j+1))*vd[l].diff(f[k])#	+sym.expand(Ga(-(l+1),i+1,j+1)).diff(f[k]) * vd[l] +Ga(-(l+1),i+1,j+1)* (sym.expand(vd[l]).diff(f[k])) 				
                    if simple==True:
                        vddd[i+j*nF+k*nF*nF] =sym.simplify(V.diff(f[i]).diff(f[j]).diff(f[k]) -FMP)
                    else:
                        vddd[i+j*nF+k*nF*nF] =V.diff(f[i]).diff(f[j]).diff(f[k]) -FMP
                    FMP=0
    else:
        for i in range(nF):
            if simple==True:
                vd[i] = sym.simplify(V.diff(f[i]))
            else:
                vd[i] = V.diff(f[i])
            for j in range(nF):
                if simple==True:
                      vdd[i+j*nF] = sym.simplify(V.diff(f[i]).diff(f[j]) )
                else:
                      vdd[i+j*nF] = V.diff(f[i]).diff(f[j])
                for k in range(nF):
                    if simple==True:
                        vddd[i+j*nF+k*nF*nF] = sym.simplify(V.diff(f[i]).diff(f[j]).diff(f[k]))
                    else:
                        vddd[i+j*nF+k*nF*nF] = V.diff(f[i]).diff(f[j]).diff(f[k])
            
    import os
    dir = os.path.dirname(__file__)
    filename1 = os.path.join(dir, 'CppTrans', 'potentialProto.h')
    filename2 = os.path.join(dir, 'CppTrans', 'potential.h')
    f = open(filename1, 'r')
    g = open(filename2, 'w')

    for line in f:
        g.write(line)
        if line == "// #Rewrite\n":
            g.write('// Potential file rewriten at' + ' ' + time.strftime("%c") +'\n')
        if line == "// #FP\n":
            g.write('nF='+str(nF)+';\n'+'nP='+str(nP)+';\n')

        if line == "// Pot\n":
            expr=str(sym.ccode(V))
            if (expr!=str(0)):
                if (expr!=str(0.0)):
                    for l in range(max(nP,nF)):
                        l=max(nP,nF)-1-l
                        expr=expr.replace("_"+str(l),"["+str(l)+"]")
                    g.write('\n sum='+str(expr)+';\n')
        if line == "// dPot\n":
            for i in range(nF):
                expr=str(sym.ccode(vd[i]))
                if (expr!=str(0)):
                    if (expr!=str(0.0)):
                        for l in range(max(nF,nP)):
                            l=max(nP,nF)-1-l
                            expr=expr.replace("_"+str(l),"["+str(l)+"]")
                        g.write('\n sum['+str(i)+']='+str(expr)+';\n')
    
        if line == "// ddPot\n":
            for i in range(nF):
                for j in range(nF):
                    expr=str(sym.ccode(vdd[i+nF*j]))
                    if (expr!=str(0)):
                        if (expr!=str(0.0)):
                            for l in range(max(nF,nP)):
                                l=max(nP,nF)-1-l
                                expr=expr.replace("_"+str(l),"["+str(l)+"]")
                            g.write('\n sum['+str(i)+'+nF*'+str(j)+']='+str(expr)+';\n')
        if line == "// dddPot\n":
            for i in range(nF):
                for j in range(nF):
                    for k in range(nF):
                        expr=str(sym.ccode(vddd[i+nF*j+nF*nF*k]))
                        if (expr!=str(0)):
                            if (expr!=str(0.0)):
                                for l in range(max(nF,nP)):
                                    l=max(nP,nF)-1-l
                                    expr=expr.replace("_"+str(l),"["+str(l)+"]")
                                g.write('\n sum['+str(i)+'+nF*'+str(j)+'+nF*nF*'+str(k)+']='+str(expr)+';\n')


    g.close()
    f.close()

def fieldmetric(G,nF,nP,simple=False):
    f=sym.symarray('f',nF)
    p=sym.symarray('p',nP)

    COR = Coordinates('\chi', f)
    g = MetricTensor('g',COR , G)
    Ga = Christoffel('Ga', g)
    Ri = Ricci('Ri', g)
    Rm = Riemann('Rm',g)
    import os
    dir = os.path.dirname(__file__)
    filename1 = os.path.join(dir, 'CppTrans', 'fieldmetricProto.h')
    filename2 = os.path.join(dir, 'CppTrans', 'fieldmetric.h')
    e = open(filename1, 'r')
    h = open(filename2, 'w')


    for line in e:
        h.write(line)
        if line == "// #FP\n":
            #h.write('nF='+str(nF)+';\n')
            h.write('nF='+str(nF)+';\n'+'nP='+str(nP)+';\n')

        if line == "// metric\n":
            for i in  range(2*nF):
                for j in  range(2*nF):
                    if i<nF:
                        ii=-i-1
                    else:
                        ii=i-(nF-1)
                    if j<nF:
                        jj=-j-1
                    else:
                        jj=j-(nF-1)
                    if simple==True:
                        expr=str(ccode(sym.simplify(g(ii,jj)*1.0)))
                    else:
                        expr=str(ccode(g(ii,jj)*1.0))
                    if g(ii,jj)!=0.0 or g(ii,jj)!=0:
                        for m in range(max(nF,nP)): 
                            expr=expr.replace("_"+str(m),"["+str(m)+"]")
                        h.write('\n FM['+str((2*nF)*i+j)+']='+str(expr)+';\n')


        if line == "// Christoffel\n":
            for i in range(2*nF):
                for j in range(2*nF):
                    for k in range(2*nF):
                        if i<nF:
                            ii=-i-1
                        else:
                            ii=i-(nF-1)
                        if j<nF:
							jj=-j-1
                        else:
						    jj=j-(nF-1)
                        if k<nF:
                            kk=-k-1
                        else:
                            kk=k-(nF-1)
                        if kk<0 or jj<0 or ii>0:
                            expr=str(0.0)
                        else:
                            if simple==True:
                                expr=str(ccode(sym.simplify(Ga(ii,jj,kk)*1.0)))
                            else:
                                expr=str(ccode(Ga(ii,jj,kk)*1.0))
                        if (expr!=str(0)):
                            if (expr!=str(0.0)):
                                for m in range(max(nF,nP)): 
                                    expr=expr.replace("_"+str(m),"["+str(m)+"]")
                                h.write('\n CS['+str((2*nF)*(2*nF)*i+(2*nF)*j+k)+']='+str(expr)+';\n')
		
    
        if line == "// Riemann\n":
                for i in range(nF):
                    for j in range(nF):
                        for k in range(nF):
                            for l in range(nF):
                                ii=i+1
                                jj=j+1
                                kk=k+1
                                ll=l+1
                                if simple==True:
                                    expr=str(ccode(sym.simplify(Rm(ii,jj,kk,ll)*1.0)))
                                else:
                                    expr=str(ccode(Rm(ii,jj,kk,ll)*1.0))
                                if (expr!=str(0)):
                                    if (expr!=str(0.0)):
                                        for m in range(max(nF,nP)): 
                                            expr=expr.replace("_"+str(m),"["+str(m)+"]")
                                        h.write('\n RM['+str((nF)*(nF)*(nF)*i+(nF)*(nF)*j+(nF)*k+l)+']='+str(expr)+';\n')
                          
       
        if line == "// Riemanncd\n":
                for i in range(nF):
                    for j in range(nF):
                        for k in range(nF):
                            for l in range(nF):
                                for m in range(nF):    
                                    ii=i+1
                                    jj=j+1
                                    kk=k+1
                                    ll=l+1
                                    mm=m+1
                                    if simple==True:
                                        expr=str(ccode(sym.simplify(Rm.covariantD(ii,jj,kk,ll,mm)*1.0)))
                                    else:
                                        expr=str(ccode(Rm.covariantD(ii,jj,kk,ll,mm)*1.0))
                                    if (expr!=str(0)):
                                        if (expr!=str(0.0)):
                                            for n in range(max(nF,nP)): 
                                                expr=expr.replace("_"+str(n),"["+str(n)+"]")
                                            h.write('\n RMcd['+str((nF)*(nF)*(nF)*(nF)*i+(nF)*(nF)*(nF)*j+(nF)*(nF)*k+(nF)*l+m)+']='+str(expr)+';\n')

    h.close()
    e.close()
    return g, Ga, Ri, Rm