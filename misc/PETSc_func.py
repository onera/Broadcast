
### A few PETSc functions
# Slightly modified from S. Beneddine & M. Lugrin

import scipy as sy
from scipy.sparse.linalg import eigs
from scipy.sparse.linalg import LinearOperator
from mpi4py import MPI
import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
from slepc4py import SLEPc
Print = PETSc.Sys.Print
import numpy as np
import sys
import copy
import timeit

def createMatIDPetsc(im, jm, comm=PETSc.COMM_WORLD):
    A = PETSc.Mat()
    A.create(comm)
    A.setSizes([im,jm])
    A.setUp()
    A.setOption(A.Option.ROW_ORIENTED,False)
    A.setType('mpiaij')
    A.setPreallocationNNZ((33*7,33*7))
    l,d=A.getVecs()
    d.set(1.)
    A.setDiagonal(d,PETSc.InsertMode.INSERT_VALUES)
    A.assemble()
    return A

def createVecFromArray(array):

    return PETSc.Vec().createWithArray(array, comm=PETSc.COMM_WORLD)

def createMatPetsc(Imat,Jmat,Amat,sizMat_i,sizMat_j,nnz, typeData=np.complex64):
  Imat = np.array(Imat, dtype=np.int32)
  Jmat = np.array(Jmat, dtype=np.int32)
  A = PETSc.Mat()
  A.create(PETSc.COMM_WORLD)
  comm = PETSc.COMM_WORLD.tompi4py()
  rank = comm.Get_rank()
  A.setSizes([sizMat_i,sizMat_j])
  A.setUp()
  A.setOption(A.Option.ROW_ORIENTED,False)
  A.setType('mpiaij')
  A.setPreallocationNNZ((33*7,33*7))
  iStart = 0
  iEnd = 0
  Print ("   ** processing mat")
  Print (" ")
  nnz_loc = len(Imat)
# we fill the matrix by chunks using setValues instead of setValue, significant increase of performance
  while iEnd < nnz_loc:
    while Jmat[iEnd] == Jmat[iStart]:
      iEnd += 1
      if iEnd >= nnz_loc:
        break
    if rank == 0:
      progress = "%.1f %% (%i out of %i)"%(100.0 * iEnd/ nnz_loc, iEnd, nnz_loc)
      sys.stdout.write(progress + chr(13))
    A.setValues(Imat[iStart:iEnd],[Jmat[iStart]],typeData(Amat[iStart:iEnd]))
    iStart = iEnd 
  Print ("   ** assembling mat                          ")
  A.assemble()
  Print ("   ** assembling done!                         \n ")
  return A


def createMatPetscCSR(Imat,Jmat,Amat,ni,nj,nnz,comm=PETSc.COMM_WORLD.tompi4py()) :
    # comm = PETSc.COMM_WORLD.tompi4py()
    rank = comm.Get_rank()
    A=PETSc.Mat()
    # A.create(PETSc.COMM_WORLD)
    A.create(comm)
    size=comm.Get_size()
    A.setSizes([ni,nj])
    A.setUp()
    A.setOption(A.Option.ROW_ORIENTED,False)
    A.setOption(A.Option.NEW_NONZERO_ALLOCATION_ERR,False)
    A.setType('mpiaij')
    A.setPreallocationNNZ((nnz,nnz))
    #passe par le format CSR via scipy pour ranger correctement
    M=sy.sparse.csr_matrix((Amat, (Imat, Jmat)), shape=(ni, nj))
    #Remplissage
    (rstart,rend)=A.getOwnershipRange()
    # print(rstart,rend)
    #set directement depuis la CSR les bons data pour chaque proc en fonction du ownership range
    A.setValuesCSR(M.indptr[rstart:rend+1] - M.indptr[rstart],M.indices[M.indptr[rstart]:M.indptr[rend]],M.data[M.indptr[rstart]:M.indptr[rend]])
    A.assemble()
    # toprint=A.getInfo()
    # if rank==0 : 
    # 	print toprint
    return A  


# def gatherVector2ArrayPetsc(vec,root=0,broadcast=False):
def gatherVector2ArrayPetsc(vec,comm=MPI.COMM_WORLD,broadcast=False):
  # comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  nProc = comm.Get_size()
  sol = None
  indi,indf = vec.getOwnershipRange()
  allIndi = comm.gather(indi,root=0)
  allIndf = comm.gather(indf,root=0)
  localVectR = np.array(vec.getArray().real, dtype=np.float64)
  localVectI = np.array(vec.getArray().imag, dtype=np.float64)
  localSize = localVectR.size
  allLocalSize = comm.gather(localSize, root=0)
  #allLocalSol = comm.gather(,root=0)
  #print "    -> gather done!"

  if rank == 0:
    sizeVec = np.amax(allIndf)
    sol = np.zeros((sizeVec),dtype=complex)
    i1, i2 = allIndi[0], allIndf[0]
    sol[i1:i2] = localVectR + 1j * localVectI
    for k in range(1,nProc):
      i1, i2 = allIndi[k], allIndf[k]
      receiveBufferR = np.empty(allLocalSize[k], dtype=np.float64)
      receiveBufferI = np.empty(allLocalSize[k], dtype=np.float64)
      comm.Recv([receiveBufferR, MPI.DOUBLE], source=k, tag=k*10)
      comm.Recv([receiveBufferI, MPI.DOUBLE], source=k, tag=k*10+1)
      sol[i1:i2] = receiveBufferR + 1j * receiveBufferI
  else:
    comm.Send([localVectR, MPI.DOUBLE], dest=0, tag=rank*10)
    comm.Send([localVectI, MPI.DOUBLE], dest=0, tag=rank*10+1)
      
  comm.Barrier()
  if broadcast:
      sol = comm.bcast(sol,root=0)
  
  return sol


def kspLUPetsc(A, comm=PETSc.COMM_WORLD):
  ksp = PETSc.KSP()
  ksp.create(comm)
  # ksp.cancelMonitor()
  ksp.setType('preonly')
  ksp.getPC().setType('lu')
  # ksp.getPC().setType('bjacobi')
  # ksp.getPC().setType('asm')
  #ksp.getPC().setFactorSolverPackage('mumps')
  ksp.getPC().setFactorSolverType('mumps')
  if A is not None:
    ksp.setOperators(A)
  ksp.setFromOptions()
  # ksp.setUp()
  #ksp.getPC().setMumpsIcntl(7, 3)
  return ksp

def kspILUPetsc(A, comm=PETSc.COMM_WORLD):
  ksp = PETSc.KSP()
  ksp.create(comm)
  # ksp.cancelMonitor()
  ksp.setType('preonly')
  ksp.getPC().setType('ilu')
  # ksp.getPC().setFactorLevels(0)
  if A is not None:
    ksp.setOperators(A)
  ksp.setFromOptions()
  return ksp

def kspGMRESplusILUPetsc(A, comm=PETSc.COMM_WORLD):
  ksp = PETSc.KSP()
  ksp.create(comm)
  # ksp.cancelMonitor()
  ksp.setType('fgmres')
  ksp.setGMRESRestart(300)
  # ksp.setTolerances(rtol=1.e-2, max_it=300)
  # ksp.setTolerances(rtol=0.5, max_it=300)
  ksp.setTolerances(rtol=1.e-4, max_it=300)
  # ksp.setTolerances(rtol=1.e-8, max_it=300)
  # ksp.setOptionsPrefix("gmresin_")
  ## 
  ksp.getPC().setType('ilu')
  # ksp.getPC().setFactorSolverType('superlu') # superlu not installed 
  ksp.getPC().setFactorLevels(0)
  ## Or with nested GMRES
  # ksp.getPC().setType('ksp')
  # kspInner = ksp.getPC().getKSP()
  # kspInner.setType('fgmres')
  # ksp.setGMRESRestart(500)
  # kspInner.setTolerances(rtol=1.e-2, max_it=10)
  # kspInner.setOptionsPrefix("gmresin2_")
  # kspInner.getPC().setType('ilu')
  # kspInner.getPC().setFactorLevels(1)
  ## 
  if A is not None:
    ksp.setOperators(A)
  ksp.setFromOptions()
  return ksp

def kspGMRESplusLUPetsc(A, comm=PETSc.COMM_WORLD):
  ksp = PETSc.KSP()
  ksp.create(comm)
  # ksp.cancelMonitor()
  ksp.setType('fgmres')
  ksp.setGMRESRestart(300)
  ksp.setTolerances(rtol=1.e-4, max_it=50)
  # ksp.setOptionsPrefix("gmresin_")
  ## 
  ksp.getPC().setType('lu')
  ksp.getPC().setFactorSolverType('mumps') 
  if A is not None:
    ksp.setOperators(A)
  ksp.setFromOptions()
  return ksp

def kspASMPetsc(A, A2):
  ksp = PETSc.KSP()
  ksp.create(PETSc.COMM_WORLD)
  ksp.setType('lgmres')
  pc = ksp.getPC()
  # ksp.setOperators(A,A2)
  ksp.setOperators(A)
  tol = 1e-3
  ksp.setTolerances(rtol=tol,max_it=1000)
  #ksp.setNormType(1) # norm preconditioned (default)
  pc.setType('asm')
#  pc.setASMType(1)
  pc.setASMOverlap(1)
  #pc.setASMTotalSubdomains(3)
  pc.setASMLocalSubdomains(1)
  pc.setUp()
  for ksp2 in pc.getASMSubKSP():
    ksp2.setType('preonly')
    # ksp2.getPC().setType('ilu')
    #ksp2.getPC().setFactorLevels(1)
    
    ksp2.getPC().setType('lu')
    #ksp2.getPC().setFactorSolverPackage('mumps')
    ksp2.getPC().setFactorSolverType('mumps')
    # ksp2.setOperators(A2)
    
    #ksp2.getPC().setFactorSolverPackage('mumps')
    #ksp2.getPC().setType('lu')
    ksp2.setFromOptions()
  ksp.setFromOptions()
  #ksp.getPC().setMumpsIcntl(7, 3)
  return ksp  



def iterNewton(residual, A, ksp, comm=None):
  """ Compute one newton iteration, i.e. 
  dW = A^{-1} * residual """
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  x, b = A.getVecs()
  rangeVec = b.getOwnershipRange()
  for k in range(rangeVec[0],rangeVec[1]):
    b[k] = residual[k]
  b.assemble()
  
  ksp.solve(b,x)
  sol = gatherVector2ArrayPetsc(x, comm, broadcast=True)
  # sol = x

  return ksp, sol


def OPinv(rhs,ksp,broadcast=True):

  A,_ = ksp.getOperators()
  x, b = A.getVecs()
  rangeVec = b.getOwnershipRange()
  for k in range(rangeVec[0],rangeVec[1]):
    b[k] = rhs[k]
  b.assemble()
  ksp.solve(b,x)
  sol = gatherVector2ArrayPetsc(x,MPI.COMM_WORLD,broadcast=broadcast)
  # sol = x
  return sol  


def printResultsEps(eps,nev=None):
  # Print()
  Print("******************************")
  Print("*** SLEPc Solution Results ***")
  Print("******************************")
  # Print()
  its = eps.getIterationNumber()
  Print("Number of iterations of the method:%d" % its)
  # eps_type = eps.getType()
  # Print("Solution method:%s" % eps_type)
  # nev, ncv, mpd = eps.getDimensions()
  # Print("Number of requested eigenvalues:%d" % nev)
  # tol, maxit = eps.getTolerances()
  # Print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
  nconv = eps.getConverged()
  Print("Number of converged eigenpairs: %d" % nconv)


############# RESOLVANT

class mon_operateur_P(object) :
    def __init__(self,ksp_R,Qq,Qvolinv,P,D,E,F,G,H): #plus vraiment besoin de Qvol maintenant
        self.ksp_R = ksp_R
        self.Qq=Qq
        self.Qvolinv=Qvolinv
        self.P=P
        self.D=D
        self.E=E
        self.F=F
        self.G=G
        self.H=H
    #mult est remplace par la resolution du systeme correspondant
    def mult(self, mat, X, Y):
        self.P.mult(X,self.D)
        self.ksp_R.solve(self.D,self.E)
        self.Qq.mult(self.E,self.F)
        self.F.conjugate()

        # self.ksp_R.solveTranspose(self.F,Y)
        # Y.conjugate()
        self.ksp_R.solveTranspose(self.F,self.G)
        self.G.conjugate()
        # self.Qvolinv.mult(self.G,self.H)
        # self.P.multTranspose(self.H,Y)
        self.P.multTranspose(self.G,Y)

class mon_operateur_P_control(object) :
    def __init__(self,ksp_R,Qq,Qvolinv,P,dRdp,D,E,F,G,H,K): #plus vraiment besoin de Qvol maintenant
        self.ksp_R = ksp_R
        self.Qq=Qq
        self.Qvolinv=Qvolinv
        self.P=P
        self.dRdp=dRdp
        self.D=D
        self.E=E
        self.F=F
        self.G=G
        self.H=H
        self.K=K
    #mult est remplace par la resolution du systeme correspondant
    def mult(self, mat, X, Y):
        self.P.mult(X,self.D)
        self.dRdp.mult(self.D,self.E)
        self.ksp_R.solve(self.E,self.F)
        self.Qq.mult(self.F,self.G)
        self.G.conjugate()
        self.ksp_R.solveTranspose(self.G,self.H)
        self.H.conjugate()
        self.dRdp.multTranspose(self.H,self.K)
        self.P.multTranspose(self.K,Y)

class mon_operateur_P_diagon(object) :
    def __init__(self,ksp_R,Qq,P,Pv,Pv_inv,D,E,F,G,H,K,L,S): #plus vraiment besoin de Qvol maintenant
        self.ksp_R = ksp_R
        self.Qq=Qq
        self.P=P
        self.Pv=Pv
        self.Pv_inv=Pv_inv
        self.D=D
        self.E=E
        self.F=F
        self.G=G
        self.H=H
        self.K=K
        self.L=L
        self.S=S
    #mult est remplace par la resolution du systeme correspondant
    def mult(self, mat, X, Y):
        self.P.mult(X,self.D)
        self.Pv_inv.mult(self.D,self.E)
        self.ksp_R.solve(self.E,self.F)
        self.Pv.mult(self.F,self.G)
        self.Qq.mult(self.G,self.H)
        self.Pv.multHermitian(self.H,self.K)
        self.K.conjugate()

        # self.ksp_R.solveTranspose(self.F,Y)
        # Y.conjugate()
        self.ksp_R.solveTranspose(self.K,self.L)
        self.L.conjugate()
        self.Pv_inv.multHermitian(self.L,self.S)
        self.P.multTranspose(self.S,Y)              

class mon_operateur_svd(object) :
    def __init__(self,ksp_R,Qvolinv,P,D,E,F): #plus vraiment besoin de Qvol maintenant
        self.ksp_R = ksp_R
        self.Qvolinv=Qvolinv
        self.P=P
        self.D=D
        self.E=E
        self.F=F
    ##mult est remplace par la resolution du systeme correspondant
    def mult(self, mat, X, Y):
        self.P.mult(X,self.D)
        self.Qvolinv.mult(self.D,self.E)
        self.ksp_R.solve(self.E,Y)
    def multHermitian(self, mat, X, Y):
        X.conjugate()
        self.ksp_R.solve(X,self.D)
        self.D.conjugate()
        self.Qvolinv.multTranspose(self.D,self.E)
        self.P.multTranspose(self.E,Y)
    ### BELOW is without Qf    
    # def mult(self, mat, X, Y): 
    #     self.P.mult(X,self.D)
    #     self.ksp_R.solve(self.D,Y)
    # def multHermitian(self, mat, X, Y):
    # 	X.conjugate()
    # 	self.ksp_R.solve(X,self.D)
    # 	self.D.conjugate()
    #     self.P.multTranspose(self.D,Y)    
            


def createShellResolvent(MAT,frequency,Qq,Qvol, P): #Qvol est plus utilise avec la nouvelle methode
        #(Iw-J)
    R_inv=-MAT
    # R_inv= MAT
    R_inv.setOption(R_inv.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    # l,d=MAT.getVecs()
    # d.set(2. * np.pi * 1.j * frequency)
    # print 'ok00'
    # R_inv.setDiagonal(d,PETSc.InsertMode.ADD_VALUES)
    # R_inv.setDiagonal(d,PETSc.InsertMode.INSERT_VALUES)
    IDE = createMatIDPetsc(R_inv.getSize()[0], R_inv.getSize()[1])
    R_inv = R_inv + 2. * np.pi * 1.j * frequency * IDE
       #pour solve
    # print 'ksp'
    # print 'ok1'
    ksp_R=kspLUPetsc(R_inv)
       #Definition de la matrice shell :
    D,E=R_inv.createVecs()
    F,G=R_inv.createVecs()
    H,K=R_inv.createVecs()
    # pde=mon_operateur(ksp_R,Qq,Qvol,P,D,E,F,G,H) #D,E,F sont juste la pour stocker des resultats intermediaires
    # pde=mon_operateur_svd(ksp_R,Qvol,P,D,E,F) #D,E,F sont juste la pour stocker des resultats intermediaires
    pde=mon_operateur_P(ksp_R,Qq,Qvol,P,D,E,F,G,H) #D,E,F sont juste la pour stocker des resultats intermediaires
    A = PETSc.Mat().create(comm=PETSc.COMM_WORLD)
    sizeMatj = P.getSize()[1]
    sizeMati = P.getSize()[1]
    # sizeMati = P.getSize()[0]
    A.setSizes([sizeMati, sizeMatj])
    A.setType('python')
    A.setPythonContext(pde)
    A.setUp()

    return A,ksp_R


def eigPetsc2(comm,A,Qvol2,P,nev,typeEps=SLEPc.EPS.Type.ARNOLDI,tol=1e-4,verbose=1,maxits=5): #Qvol =norme L2
  rank = comm.Get_rank()
  size = comm.Get_size()
  eps = SLEPc.EPS()
  eps.create(PETSc.COMM_WORLD)
  eps.setType(typeEps)
  eps.setProblemType(SLEPc.EPS.ProblemType.GHEP) #GHEP pour pouvoir passer Qvol directement
  # eps.setProblemType(SLEPc.EPS.ProblemType.GNHEP)
  # eps.setProblemType(SLEPc.EPS.ProblemType.NHEP)
  # eps.setOperators(A,Qvol2.matMult(P)) 
  # eps.setOperators(A,Qvol2) #Qvol en second argument pour la norme l2
  t00 = timeit.time.time()
  eps.setOperators(A,P.transposeMatMult(Qvol2.matMult(P)))
  t11 = timeit.time.time()
  t_setop = t11 - t00
  # print 't set operator : ', t_setop
  # eps.setOperators(A) 
  eps.setDimensions(nev=nev)
  eps.setWhichEigenpairs(1) # largest magnitude
  # eps.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_MAGNITUDE) 
  eps.setTolerances(tol=tol,max_it=maxits)
  # eps.setTolerances(tol=1.e-5,max_it=10)
  eps.setTrueResidual(True) # Do not base convergence on residual of transform eigenproblem (trough shift-invert)
  # eps.setConvergenceTest(SLEPc.EPS.Conv.ABS)
  eps.setFromOptions()
  Print ("   ** solving the eigenvalue problem...                         ")
  import psutil
  eps.solve()
  if rank==0:
    print(psutil.virtual_memory())
  t22 = timeit.time.time()
  t_solve = t22 - t11
  # print 't solve eig : ', t_solve
  Print ("   ** solved!      \n                   ")
  # if verbose > 0:
  printResultsEps(eps)
  nconv = eps.getConverged()
  xr, tmp = A.getVecs()
  # vec, tmp = A.getVecs()
  vals = []
  vecs = []
  # for i in range(nconv):
  for i in range(nev):	
    val = eps.getEigenpair(i, xr)
    vec = gatherVector2ArrayPetsc(xr,comm,broadcast=True)
    if rank == 0:
    # vec = vecR + 1j * vecI
      print(val)
      # vals.append(val)
      # vecs.append(vec)
    vals.append(val)
    vecs.append(vec)

  comm.Barrier()
  # print 'vecs = ', vecs
  return vals, vecs, eps
  # if rank == 0:
  #   return vals,vecs,eps
  # else:
  #   return None, None, eps


def computeReponse(comm, eigenvalue, eigenvector_forcage, ksp_R, P):
    comm.Barrier()
    rank = comm.Get_rank()
    eigenvector_forcage = comm.bcast(eigenvector_forcage, root=0)
    eigenvalue = comm.bcast(eigenvalue, root=0)
    eigenvector_reponse=copy.deepcopy(eigenvector_forcage)
    eigenvector_forcageP=copy.deepcopy(eigenvector_forcage)
    for idMode in np.arange(0,len(eigenvalue)):
        # f_pet = PETSc.Vec().createWithArray(eigenvector_forcage[idMode], size=np.shape(eigenvector_forcage[idMode])[0], comm=PETSc.COMM_WORLD)
        Pf,Pf2 = P.getVecs()
        rangePf = Pf.getOwnershipRange()
        for k in range(rangePf[0], rangePf[1]):
            Pf[k] = eigenvector_forcage[idMode][k]
        Pf.assemble()
        # print Pf.getSizes()
        P.mult(Pf, Pf2)
        eigenvector_reponse[idMode] = OPinv(Pf2,ksp_R)
        Pf2 = gatherVector2ArrayPetsc(Pf2,comm,broadcast=True)
        eigenvector_forcageP[idMode]= Pf2

    return eigenvector_reponse, eigenvector_forcageP


def computeReponse2(comm, eigenvalue, eigenvector_forcage, ksp_R, P):
    comm.Barrier()
    rank = comm.Get_rank()
    eigenvector_forcage = comm.bcast(eigenvector_forcage[0], root=0)
    Pf,Pf2 = P.getVecs()
    rangePf = Pf.getOwnershipRange()
    for k in range(rangePf[0], rangePf[1]):
        Pf[k] = eigenvector_forcage[k]
    Pf.assemble()
    # rangePf2 = Pf2.getOwnershipRange()
    # for k in range(rangePf2[0], rangePf2[1]):
    #     Pf2[k] = eigenvector_forcage[k]
    # Pf2.assemble()
    # print Pf.getSizes()
    P.mult(Pf, Pf2)
    eigenvector_reponse  = OPinv(Pf2,ksp_R)
    eigenvector_forcageP = gatherVector2ArrayPetsc(Pf2,comm,broadcast=True)
    # eigenvector_reponse = eigenvector_forcageP

    return eigenvector_reponse, eigenvector_forcageP


def computeReponse2_control(comm, eigenvalue, eigenvector_forcage, ksp_R, P, dRdp):
    comm.Barrier()
    rank = comm.Get_rank()
    eigenvector_forcage = comm.bcast(eigenvector_forcage[0], root=0)
    Pf,Pf2 = P.getVecs()
    x,y = dRdp.getVecs()
    rangePf = Pf.getOwnershipRange()
    for k in range(rangePf[0], rangePf[1]):
        Pf[k] = eigenvector_forcage[k]
    Pf.assemble()
    # rangePf2 = Pf2.getOwnershipRange()
    # for k in range(rangePf2[0], rangePf2[1]):
    #     Pf2[k] = eigenvector_forcage[k]
    # Pf2.assemble()
    # print Pf.getSizes()
    P.mult(Pf, Pf2)
    dRdp.mult(Pf2, y)
    eigenvector_reponse  = OPinv(y,ksp_R)
    eigenvector_forcageP = gatherVector2ArrayPetsc(Pf2,comm,broadcast=True)
    # eigenvector_reponse = eigenvector_forcageP

    return eigenvector_reponse, eigenvector_forcageP


def computeResponseSB(ksp_R,forcingMode,broadcast=True):
    AA,_ = ksp_R.getOperators()
    u, _ = AA.getVecs()
    f = PETSc.Vec().createMPI(forcingMode.size)
    rangeVec = f.getOwnershipRange()
    for k in range(rangeVec[0],rangeVec[1]):
      f[k] = forcingMode[k]
    f.assemble()
    ksp_R.solve(f,u)
    sol = gatherVector2ArrayPetsc(u,MPI.COMM_WORLD,broadcast=broadcast)
    return sol


def svdPetsc2(comm,A,P,nev,tol=1e-4,verbose=1,maxits=5): #Qvol =norme L2
  rank = comm.Get_rank()
  size = comm.Get_size()
  svd = SLEPc.SVD()
  svd.create(PETSc.COMM_WORLD)
  # svd.setType(SLEPc.SVD.Type.CYCLIC)
  # svd.setType(SLEPc.SVD.Type.LANCZOS)
  t00 = timeit.time.time()
  svd.setOperator(A)
  t11 = timeit.time.time()
  t_setop = t11 - t00
  svd.setImplicitTranspose(True)
  # print 't set operator : ', t_setop
  # eps.setOperators(A) 
  svd.setDimensions(nsv=nev)
  svd.setFromOptions()
  Print ("   ** solving the eigenvalue problem...                         ")
  svd.solve()
  t22 = timeit.time.time()
  t_solve = t22 - t11
  # print 't solve svd : ', t_solve
  Print ("   ** solved!      \n                   ")
  # printResultsEps(eps)
  nconv = svd.getConverged()
  vecv, vecu = A.getVecs()
  # xi, tmp = A.getVecs()
  vals = []
  vecsu = []
  vecsv = []
  # for i in range(nconv):
  for i in range(1):	
    val = svd.getSingularTriplet(i, vecu, vecv)
    vecu = gatherVector2ArrayPetsc(vecu,comm,broadcast=True)
    vecv = gatherVector2ArrayPetsc(vecv,comm,broadcast=True)
    # print vec
    if rank == 0:
    # vec = vecR + 1j * vecI
      print(val)
    vals.append(val)
    vecsu.append(vecu)
    vecsv.append(vecv)
  comm.Barrier()
  # vecsv = comm.bcast(vecsv, root=0)
  # eigenvector_forcageP=copy.deepcopy(vecsv)
  # for idMode in np.arange(0,len(vals)):
  #      f_pet = PETSc.Vec().createWithArray(vecsv[idMode], size=np.shape(vecsv[idMode])[0], comm=PETSc.COMM_WORLD)
  #      Pf2,Pf=A.createVecs()
  #      P.mult(f_pet, Pf)
  #      Pf = gatherVector2ArrayPetsc(Pf,comm,broadcast=True)
  #      eigenvector_forcageP[idMode]= Pf
  return vals, vecsu, vecs, svd
  # return vals, vecsu, eigenvector_forcageP, svd



#################################################
#################################################
## 3D ##

def createShellResolvent3D(MAT, frequency, wavenumber, Qq, Qvol, P, Dz, Dzz, Pv, Pv_inv): #Qvol est plus utilise avec la nouvelle methode
        #(Iw-J)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ranknum = comm.Get_size() 
    R_inv=-MAT
    # R_inv= MAT
    R_inv.setOption(R_inv.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    l,d=MAT.getVecs()
    d.set(2. * np.pi * 1.j * frequency)
    R_inv.setDiagonal(d,PETSc.InsertMode.ADD_VALUES)
       #pour solve  
    ## 1st method : DIRECT
    # M11 = R_inv - wavenumber/3**0.5 * Dz + wavenumber**2 * Dzz
    # M12 = - 2*wavenumber/3**0.5 * Dz
    # M21 = 2*wavenumber/3**0.5 * Dz
    # M22 = R_inv + wavenumber/3**0.5 * Dz + wavenumber**2 * Dzz
    # R_inv22 = createMatPetscBlockCSR_2x2_diag2mat(M11, M22, R_inv.getSize()[0], R_inv.getSize()[1], 2*5*11**2)
    # R_inv2 = createMatPetscBlockCSR_2x2_offdiag2mat(M12, M21, R_inv.getSize()[0], R_inv.getSize()[1], 2*5*11**2)
    # R_inv2.aypx(1., R_inv22)
    ## 2nd method : DIAGONALISED 2x2
    # M11 = R_inv + wavenumber**2 * Dzz - wavenumber * 1.j * Dz
    # M22 = R_inv + wavenumber**2 * Dzz + wavenumber * 1.j * Dz
    tinit5 = timeit.time.time()
    # R_inv2 = createMatPetscBlockCSR_2x2_diag2mat(M11, M22, R_inv.getSize()[0], R_inv.getSize()[1], 5*11**2)
    ## DIAGONALISED 3x3
    # R_inv2 = createMatPetscBlockCSR_NxN_diagNmat([R_inv, M11, M22], R_inv.getSize()[0], R_inv.getSize()[1], 5*11**2, 3)

    # TSM_mat = TSM_circulant3(R_inv, Dz, Dzz, wavenumber)
    if ranknum %2 == 0:
        TSM_mat = TSM_circulant_Neven(R_inv, Dz, Dzz, wavenumber, ranknum)
    else:
        TSM_mat = TSM_circulant_Nodd(R_inv, Dz, Dzz, wavenumber, ranknum)    
    listmat = eigen_circulant(TSM_mat)
    R_inv2 = createMatPetscBlockCSR_NxN_diagNmat(listmat, R_inv.getSize()[0], R_inv.getSize()[1], 5*11**2, ranknum)

    tinit6 = timeit.time.time()
    print('Time R-1: ', tinit6 - tinit5)

    PETSc.COMM_WORLD.tompi4py().barrier()
    ksp_R=kspLUPetsc(R_inv2)
    # ksp_R=kspASMPetsc(R_inv2, R_inv22) #Only Diag as preconditionning
       #Definition de la matrice shell :
    D,E=R_inv2.createVecs()
    F,G=R_inv2.createVecs()
    H,K=R_inv2.createVecs()
    L,S=R_inv2.createVecs()
    ## 1st method : DIRECT
    # pde=mon_operateur_P(ksp_R,Qq,Qvol,P,D,E,F,G,H) #D,E,F sont juste la pour stocker des resultats intermediaires
    ## 2nd method : DIAGONALISED
    pde=mon_operateur_P_diagon(ksp_R,Qq,P,Pv,Pv_inv,D,E,F,G,H,K,L,S) #D,E,F sont juste la pour stocker des resultats intermediaires

    A = PETSc.Mat().create(comm=PETSc.COMM_WORLD)
    sizeMatj = P.getSize()[1]
    sizeMati = P.getSize()[1]
    # sizeMati = P.getSize()[0]
    A.setSizes([sizeMati, sizeMatj])
    A.setType('python')
    A.setPythonContext(pde)
    A.setUp()

    return A,ksp_R

def TSM_circulant3(R_inv, Dz, Dzz, wavenumber):
    ''' Compute the 1st line of a 3x3 block circulant matrix computed with TSM coeffs '''
    TSM_mat = [R_inv + 2./3 * wavenumber**2 * Dzz , - wavenumber/3**0.5 * Dz - 1./3 * wavenumber**2 * Dzz , wavenumber/3**0.5 * Dz - 1./3 * wavenumber**2 * Dzz]
    return TSM_mat

def TSM_circulant_Neven(R_inv, Dz, Dzz, wavenumber, ranknum):
    ''' Compute the 1st line of a NxN (N even number) block circulant matrix computed with TSM coeffs '''
    TSM_mat = [R_inv + (ranknum -1)*(ranknum - 2)/12. * wavenumber**2 * Dzz]
    for k in range(1, ranknum):
        TSM_mat.append( -0.5 * (-1)**k * 1./np.tan(np.pi / ranknum * k) * wavenumber * Dz - (ranknum/4. * (-1)**k + 0.5 * (-1)**(k+1) * 1./np.sin(np.pi / ranknum * k)**2) * wavenumber**2 * Dzz )
    return TSM_mat  

def TSM_circulant_Nodd(R_inv, Dz, Dzz, wavenumber, ranknum):
    ''' Compute the 1st line of a NxN (N odd number) block circulant matrix computed with TSM coeffs '''
    TSM_mat = [R_inv + (ranknum**2 -1)/12. * wavenumber**2 * Dzz]
    for k in range(1, ranknum):
        TSM_mat.append( -0.5 * (-1)**k * 1./np.sin(np.pi / ranknum * k) * wavenumber * Dz - 0.5 * (-1)**(k+1) * np.cos(np.pi / ranknum * k) / np.sin(np.pi / ranknum * k)**2 * wavenumber**2 * Dzz )
    return TSM_mat       

def eigen_circulant(TSM_mat):
    ''' Compute the eigenvalues of a block circulant matrix '''  
    ranknum = len(TSM_mat)
    z = np.exp(1.j * 2 * np.pi / ranknum)
    listmat = []
    for k in range(ranknum):
        summ = TSM_mat[0]
        for j in range(1, ranknum):
            summ = summ + z**(k*j) * TSM_mat[j]
        listmat.append(summ)
    return listmat

def createShellResolvent3D_1B(MAT, frequency, wavenumber, Qq, Qvol, P, Dz, Dzz): #Qvol est plus utilise avec la nouvelle methode
        #(Iw-J)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    R_inv=-MAT
    # R_inv= MAT
    R_inv.setOption(R_inv.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    l,d=MAT.getVecs()
    d.set(2. * np.pi * 1.j * frequency)
    R_inv.setDiagonal(d,PETSc.InsertMode.ADD_VALUES)
       #pour solve  
    R_inv2 = R_inv + wavenumber**2 * Dzz - wavenumber * 1.j * Dz #good
    # R_inv2 = R_inv - wavenumber**2 * Dzz + wavenumber * 1.j * Dz
    ## Direct
    ksp_R=kspLUPetsc(R_inv2)
    ## Or inner GMRES
    # ksp_R=kspGMRESplusLUPetsc(R_inv2)
    ##
       #Definition de la matrice shell :
    D,E=R_inv2.createVecs()
    F,G=R_inv2.createVecs()
    H,K=R_inv2.createVecs()
    # L,S=R_inv2.createVecs()
    pde=mon_operateur_P(ksp_R,Qq,Qvol,P,D,E,F,G,H) #D,E,F sont juste la pour stocker des resultats intermediaires
    A = PETSc.Mat().create(comm=PETSc.COMM_WORLD)
    sizeMatj = P.getSize()[1]
    sizeMati = P.getSize()[1]
    # sizeMati = P.getSize()[0]
    A.setSizes([sizeMati, sizeMatj])
    A.setType('python')
    A.setPythonContext(pde)
    A.setUp()

    return A,ksp_R    

def createShellResolvent3D_1B_control(MAT, frequency, wavenumber, Qq, Qvol, P, Dz, Dzz, dRdp): #Qvol est plus utilise avec la nouvelle methode
        #(Iw-J)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    R_inv=-MAT
    # R_inv= MAT
    R_inv.setOption(R_inv.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    l,d=MAT.getVecs()
    d.set(2. * np.pi * 1.j * frequency)
    R_inv.setDiagonal(d,PETSc.InsertMode.ADD_VALUES)
       #pour solve  
    R_inv2 = R_inv + wavenumber**2 * Dzz - wavenumber * 1.j * Dz #good
    # R_inv2 = R_inv - wavenumber**2 * Dzz + wavenumber * 1.j * Dz
    ## Direct
    ksp_R=kspLUPetsc(R_inv2)
    ## Or inner GMRES
    # ksp_R=kspGMRESplusLUPetsc(R_inv2)
    ##
       #Definition de la matrice shell :
    D,DD=dRdp.createVecs()
    K,KK=dRdp.createVecs()
    F,G=R_inv2.createVecs()
    H,E=R_inv2.createVecs()
    # L,S=R_inv2.createVecs()
    pde=mon_operateur_P_control(ksp_R,Qq,Qvol,P,dRdp,D,E,F,G,H,K) #D,E,F sont juste la pour stocker des resultats intermediaires
    A = PETSc.Mat().create(comm=PETSc.COMM_WORLD)
    sizeMatj = P.getSize()[1]
    sizeMati = P.getSize()[1]
    # sizeMati = P.getSize()[0]
    A.setSizes([sizeMati, sizeMatj])
    A.setType('python')
    A.setPythonContext(pde)
    A.setUp()

    return A,ksp_R  

def communicateIAsquare(arr, comm, start=0, IA=False):
  ''' gather from rank 1 to rank 0 then send a copy of the complete information stored in rank 0 to rank 1, both ranks get the complete information '''

  rank = comm.Get_rank()

  if rank==1:
    comm.Send(arr[start:], dest=0, tag=1)

  if rank==0:
    arr2 = np.empty_like(arr[start:])
    comm.Recv(arr2, source=1, tag=1)
    if IA:
      arr2 = arr2 + arr[-1]
    arr_new  = np.concatenate((arr, arr2))
    comm.Send(arr_new, dest=1, tag=2)

  if rank==1:
    arr2 = np.empty_like(arr[start:])
    arr_new  = np.concatenate((arr, arr2))
    comm.Recv(arr_new, source=0, tag=2)

  return arr_new

def communicateIA(arr, comm, start=0, IA=False):
  rank = comm.Get_rank()

  if rank==1:
    comm.send(arr[start:], dest=0, tag=1)

  if rank==0:
    arr2 = comm.recv(source=1, tag=1)
    if IA:
      arr2 = arr2 + arr[-1]
    arr_new  = np.concatenate((arr, arr2))
    comm.send(arr_new, dest=1, tag=2)

  if rank==1:
    arr_new = comm.recv(source=0, tag=2)

  return arr_new      

def createMatPetscBlockCSR_2x2_diag(MAT,ni,nj,nnz, square=True) :
    comm = PETSc.COMM_WORLD.tompi4py()
    rank = comm.Get_rank()
    A=PETSc.Mat()
    A.create(PETSc.COMM_WORLD)
    size=comm.Get_size()
    A.setSizes([2*ni,2*nj])
    A.setUp()
    A.setOption(A.Option.ROW_ORIENTED,False)
    # A.setOption(A.Option.NEW_NONZERO_ALLOCATION_ERR,False)
    # A.setType('mpibaij')
    A.setType('mpiaij')
    A.setPreallocationNNZ((nnz,nnz))
    IA, JA, Aij = MAT.getValuesCSR()
    # print np.shape(IA)
    # print np.shape(JA)
    # print np.shape(Aij)
    if square:
      comm.Barrier()
      IA_new = communicateIAsquare(IA, comm, start=1, IA=True)
      comm.Barrier()
      JA_new = communicateIAsquare(JA, comm)
      comm.Barrier()
      Aij_new = communicateIAsquare(Aij, comm)
    else:
      comm.Barrier()
      IA_new = communicateIA(IA, comm, start=1, IA=True)
      comm.Barrier()
      JA_new = communicateIA(JA, comm)
      comm.Barrier()
      Aij_new = communicateIA(Aij, comm)  
    if rank==1:
      JA_new = JA_new + nj

    # print np.shape(IA)
    # print np.shape(JA)
    # print np.shape(Aij)
    # print np.shape(IA_new)
    # print np.shape(JA_new)
    # print np.shape(Aij_new)
    comm.Barrier()
    #Remplissage
    # (rstart,rend)=A.getOwnershipRange()
    rstart=0
    rend=2*(np.shape(IA)[0]-1)
    #set directement depuis la CSR les bons data pour chaque proc en fonction du ownership range
    # A.setValuesBlockedCSR(IA_new[rstart:rend+1], JA_new[IA_new[rstart]:IA_new[rend]], Aij_new[IA_new[rstart]:IA_new[rend]])
    A.setValuesCSR(IA_new[rstart:rend+1], JA_new[IA_new[rstart]:IA_new[rend]], Aij_new[IA_new[rstart]:IA_new[rend]])
    A.assemble()
    # toprint=A.getInfo()
    # if rank==0 : 
    #   print toprint
    return A   

def gatherIAsquare(arr, comm, start=0, IA=False):
  ''' Gather rank 1 information in rank 0 '''
  rank = comm.Get_rank()

  if rank==1:
    comm.Send(arr[start:], dest=0, tag=1)

  if rank==0:
    arr2 = np.empty_like(arr[start:])
    comm.Recv(arr2, source=1, tag=1)
    if IA:
      arr2 = arr2 + arr[-1]
    arr_new  = np.concatenate((arr, arr2))

    return arr_new

def gatherIA(arr, comm, start=0, IA=False):
  ''' Gather rank 1 information in rank 0 '''
  rank = comm.Get_rank()

  if rank==1:
    comm.send(arr[start:], dest=0, tag=1)

  if rank==0:
    arr2 = comm.recv(source=1, tag=1)
    if IA:
      arr2 = arr2 + arr[-1]
    arr_new  = np.concatenate((arr, arr2))

    return arr_new           

def createMatPetscBlockCSR_2x2_diag2mat(MAT11, MAT22,ni,nj,nnz) :
    comm = PETSc.COMM_WORLD.tompi4py()
    rank = comm.Get_rank()
    A=PETSc.Mat()
    A.create(PETSc.COMM_WORLD)
    size=comm.Get_size()
    A.setSizes([2*ni,2*nj])
    A.setUp()
    A.setOption(A.Option.ROW_ORIENTED,False)
    # A.setOption(A.Option.NEW_NONZERO_ALLOCATION_ERR,False)
    # A.setType('mpibaij')
    A.setType('mpiaij')
    A.setPreallocationNNZ((nnz,nnz))

    IA11, JA11, Aij11 = MAT11.getValuesCSR()
    comm.Barrier()
    IA_new1 = gatherIA(IA11, comm, start=1, IA=True)
    comm.Barrier()
    JA_new1 = gatherIA(JA11, comm)
    comm.Barrier()
    Aij_new1 = gatherIA(Aij11, comm)

    IA22, JA22, Aij22 = MAT22.getValuesCSR()
    comm.Barrier()
    IA_new4 = communicateIA(IA22, comm, start=1, IA=True)
    comm.Barrier()
    JA_new4 = communicateIA(JA22, comm)
    comm.Barrier()
    Aij_new4 = communicateIA(Aij22, comm)

    if rank ==0:
      IA_new = IA_new1
      JA_new = JA_new1
      Aij_new = Aij_new1

    if rank ==1:
      IA_new = IA_new4
      JA_new = JA_new4 + nj
      Aij_new = Aij_new4 

    # print np.shape(IA)
    # print np.shape(JA)
    # print np.shape(Aij)
    # print np.shape(IA_new)
    # print np.shape(JA_new)
    # print np.shape(Aij_new)
    comm.Barrier()
    #Remplissage
    # (rstart,rend)=A.getOwnershipRange()
    rstart=0
    rend=np.shape(IA_new)[0]-1
    #set directement depuis la CSR les bons data pour chaque proc en fonction du ownership range
    # A.setValuesBlockedCSR(IA_new[rstart:rend+1], JA_new[IA_new[rstart]:IA_new[rend]], Aij_new[IA_new[rstart]:IA_new[rend]])
    A.setValuesCSR(IA_new[rstart:rend+1], JA_new[IA_new[rstart]:IA_new[rend]], Aij_new[IA_new[rstart]:IA_new[rend]])
    A.assemble()
    # toprint=A.getInfo()
    # if rank==0 : 
    #   print toprint
    return A


def createMatPetscBlockCSR_2x2_offdiag2mat(MAT12, MAT21,ni,nj,nnz) :
    comm = PETSc.COMM_WORLD.tompi4py()
    rank = comm.Get_rank()
    A=PETSc.Mat()
    A.create(PETSc.COMM_WORLD)
    size=comm.Get_size()
    A.setSizes([2*ni,2*nj])
    A.setUp()
    A.setOption(A.Option.ROW_ORIENTED,False)
    # A.setOption(A.Option.NEW_NONZERO_ALLOCATION_ERR,False)
    # A.setType('mpibaij')
    A.setType('mpiaij')
    A.setPreallocationNNZ((nnz,nnz))

    IA12, JA12, Aij12 = MAT12.getValuesCSR()
    comm.Barrier()
    IA_new2 = gatherIA(IA12, comm, start=1, IA=True)
    comm.Barrier()
    JA_new2 = gatherIA(JA12, comm)
    comm.Barrier()
    Aij_new2 = gatherIA(Aij12, comm)

    IA21, JA21, Aij21 = MAT21.getValuesCSR()
    comm.Barrier()
    IA_new3 = communicateIA(IA21, comm, start=1, IA=True)
    comm.Barrier()
    JA_new3 = communicateIA(JA21, comm)
    comm.Barrier()
    Aij_new3 = communicateIA(Aij21, comm)

    if rank ==0:
      IA_new = IA_new2
      JA_new = JA_new2 + nj
      Aij_new = Aij_new2

    if rank ==1:
      IA_new = IA_new3
      JA_new = JA_new3
      Aij_new = Aij_new3

    # print np.shape(IA)
    # print np.shape(JA)
    # print np.shape(Aij)
    # print np.shape(IA_new)
    # print np.shape(JA_new)
    # print np.shape(Aij_new)
    comm.Barrier()
    #Remplissage
    # (rstart,rend)=A.getOwnershipRange()
    rstart=0
    rend=np.shape(IA_new)[0]-1
    #set directement depuis la CSR les bons data pour chaque proc en fonction du ownership range
    # A.setValuesBlockedCSR(IA_new[rstart:rend+1], JA_new[IA_new[rstart]:IA_new[rend]], Aij_new[IA_new[rstart]:IA_new[rend]])
    A.setValuesCSR(IA_new[rstart:rend+1], JA_new[IA_new[rstart]:IA_new[rend]], Aij_new[IA_new[rstart]:IA_new[rend]])
    A.assemble()
    # toprint=A.getInfo()
    # if rank==0 : 
    #   print toprint
    return A    

def createJacEigenVectorsMat(ni, nj):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    I = createMatIDPetsc(ni, nj)

    M11 = I
    M12 = I
    z = -0.5 + 1.j * 0.5 * 3**0.5
    M21 = z * I
    z_conj = -0.5 - 1.j * 0.5 * 3**0.5
    M22 = z_conj * I
    A2  = createMatPetscBlockCSR_2x2_diag2mat(M11, M22, ni, nj, 2)
    A22 = createMatPetscBlockCSR_2x2_offdiag2mat(M12, M21, ni, nj, 2)
    A2.aypx(1., A22)
    # A2.scale(1./2**0.5)
    A2 = A2 /2**0.5

    return A2

def createJacEigenVectorsMatInv(ni, nj): 
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    I = createMatIDPetsc(ni, nj)

    zconj = -0.5 - 1.j * 0.5 * 3**0.5
    M11 = zconj * I
    M12 = - I
    z = -0.5 + 1.j * 0.5 * 3**0.5
    M21 = - z * I
    M22 = I
    A2  = createMatPetscBlockCSR_2x2_diag2mat(M11, M22, ni, nj, 2)
    A22 = createMatPetscBlockCSR_2x2_offdiag2mat(M12, M21, ni, nj, 2)
    A2.aypx(1., A22)
    # A2.scale(1.j * 2**0.5 / 3**0.5)
    A2 = A2 * (1.j * 2**0.5 / 3**0.5)

    return A2    

def sendIA_k(arr, comm, ranknum):
  ''' send the complete information stored in rank 0 to rank ranknum, both ranks get the complete information '''
  rank = comm.Get_rank()

  if rank==0:
    comm.send(arr, dest=ranknum, tag=ranknum)

    return arr

  if rank==ranknum:
    arr_new = comm.recv(source=0, tag=ranknum)

    return arr_new  


def gatherIA_k(arr, comm, ranknum, start=0, IA=False):
  ''' Gather rank number ranknum information in rank 0 '''
  rank = comm.Get_rank()

  if rank==ranknum:
    comm.send(arr[start:], dest=0, tag=ranknum)

  if rank==0:
    arr2 = comm.recv(source=ranknum, tag=ranknum)
    if IA:
      arr2 = arr2 + arr[-1]
    arr_new  = np.concatenate((arr, arr2))

    return arr_new  

def createMatPetscBlockCSR_NxN_diag(MAT,ni,nj,nnz, ranknum):
    ''' Create a diagonal matrix with the same block matrix repeated on the diagonal in a big sparse (ni*ranknum X nj*ranknum) matrix '''
    comm = PETSc.COMM_WORLD.tompi4py()
    rank = comm.Get_rank()
    A=PETSc.Mat()
    A.create(PETSc.COMM_WORLD)
    size=comm.Get_size()
    A.setSizes([ranknum*ni,ranknum*nj])
    A.setUp()
    A.setOption(A.Option.ROW_ORIENTED,False)
    # A.setOption(A.Option.NEW_NONZERO_ALLOCATION_ERR,False)
    # A.setType('mpibaij')
    A.setType('mpiaij')
    A.setPreallocationNNZ((nnz,nnz))
    IA_new, JA_new, Aij_new = MAT.getValuesCSR()
    for k in range(1, ranknum):
        comm.Barrier()
        if rank==k or rank==0:
            IA_new = gatherIA_k(IA_new, comm, k, start=1, IA=True)
            JA_new = gatherIA_k(JA_new, comm, k)
            Aij_new = gatherIA_k(Aij_new, comm, k)    
    for k in range(1, ranknum):
        comm.Barrier()
        if rank==k or rank==0:
            IA_new = sendIA_k(IA_new, comm, k)
            JA_new = sendIA_k(JA_new, comm, k)
            Aij_new = sendIA_k(Aij_new, comm, k) 
            if rank == k:
                JA_new = JA_new + k * nj
    # print np.shape(IA_new)
    # print np.shape(JA_new)
    # print np.shape(Aij_new)
    comm.Barrier()
    #Remplissage
    # (rstart,rend)=A.getOwnershipRange()
    rstart=0
    rend=ni
    #set directement depuis la CSR les bons data pour chaque proc en fonction du ownership range
    # A.setValuesBlockedCSR(IA_new[rstart:rend+1], JA_new[IA_new[rstart]:IA_new[rend]], Aij_new[IA_new[rstart]:IA_new[rend]])
    A.setValuesCSR(IA_new[rstart:rend+1], JA_new[IA_new[rstart]:IA_new[rend]], Aij_new[IA_new[rstart]:IA_new[rend]])
    A.assemble()
    # toprint=A.getInfo()
    # if rank==0 : 
    #   print toprint
    return A

def createMatPetscBlockCSR_NxN_column(MATlist,ni,nj,nnz, ranknum, column):
    ''' Fill in only 1 column (e.g. 1st column if column==0) of a big sparse (ni*ranknum X nj*ranknum) matrix '''
    comm = PETSc.COMM_WORLD.tompi4py()
    rank = comm.Get_rank()
    A=PETSc.Mat()
    A.create(PETSc.COMM_WORLD)
    size=comm.Get_size()
    A.setSizes([ranknum*ni,ranknum*nj])
    A.setUp()
    A.setOption(A.Option.ROW_ORIENTED,False)
    # A.setOption(A.Option.NEW_NONZERO_ALLOCATION_ERR,False)
    # A.setType('mpibaij')
    A.setType('mpiaij')
    A.setPreallocationNNZ((nnz,nnz))
    IA_saved = []
    JA_saved = []
    Aij_saved = []
    for i in range(len(MATlist)):
        comm.Barrier()
        MAT = MATlist[i]
        IA_new, JA_new, Aij_new = MAT.getValuesCSR()
        for k in range(1, ranknum):
            comm.Barrier()
            if rank==k or rank==0:
                IA_new = gatherIA_k(IA_new, comm, k, start=1, IA=True)
                JA_new = gatherIA_k(JA_new, comm, k)
                Aij_new = gatherIA_k(Aij_new, comm, k)
        comm.Barrier()
        IA_saved.append(IA_new)
        JA_saved.append(JA_new)
        Aij_saved.append(Aij_new)
    for i in range(1, len(MATlist)):
        comm.Barrier()
        if rank==i or rank==0:
            IA_new = sendIA_k(IA_saved[i], comm, i)
            JA_new = sendIA_k(JA_saved[i], comm, i)
            Aij_new = sendIA_k(Aij_saved[i], comm, i) 
    if rank == 0:
        IA_new = IA_saved[0]
        JA_new = JA_saved[0]
        Aij_new = Aij_saved[0]

    JA_new = JA_new + column * nj
    # print np.shape(IA_new)
    # print np.shape(JA_new)
    # print np.shape(Aij_new)
    comm.Barrier()
    #Remplissage
    # (rstart,rend)=A.getOwnershipRange()
    rstart=0
    rend=ni
    #set directement depuis la CSR les bons data pour chaque proc en fonction du ownership range
    # A.setValuesBlockedCSR(IA_new[rstart:rend+1], JA_new[IA_new[rstart]:IA_new[rend]], Aij_new[IA_new[rstart]:IA_new[rend]])
    A.setValuesCSR(IA_new[rstart:rend+1], JA_new[IA_new[rstart]:IA_new[rend]], Aij_new[IA_new[rstart]:IA_new[rend]])
    A.assemble()
    # toprint=A.getInfo()
    # if rank==0 : 
    #   print toprint
    return A       

def createMatPetscBlockCSR_NxN_diagNmat(MATlist,ni,nj,nnz, ranknum):
    ''' Fill in only the block matrices on the diagonal of a big sparse (ni*ranknum X nj*ranknum) matrix '''
    comm = PETSc.COMM_WORLD.tompi4py()
    rank = comm.Get_rank()
    A=PETSc.Mat()
    A.create(PETSc.COMM_WORLD)
    size=comm.Get_size()
    A.setSizes([ranknum*ni,ranknum*nj])
    A.setUp()
    A.setOption(A.Option.ROW_ORIENTED,False)
    # A.setOption(A.Option.NEW_NONZERO_ALLOCATION_ERR,False)
    # A.setType('mpibaij')
    A.setType('mpiaij')
    A.setPreallocationNNZ((nnz,nnz))
    IA_saved = []
    JA_saved = []
    Aij_saved = []
    for i in range(len(MATlist)):
        comm.Barrier()
        MAT = MATlist[i]
        IA_new, JA_new, Aij_new = MAT.getValuesCSR()
        for k in range(1, ranknum):
            comm.Barrier()
            if rank==k or rank==0:
                IA_new = gatherIA_k(IA_new, comm, k, start=1, IA=True)
                JA_new = gatherIA_k(JA_new, comm, k)
                Aij_new = gatherIA_k(Aij_new, comm, k)
        comm.Barrier()
        IA_saved.append(IA_new)
        JA_saved.append(JA_new)
        Aij_saved.append(Aij_new)
    for i in range(1, len(MATlist)):
        comm.Barrier()
        if rank==i or rank==0:
            IA_new = sendIA_k(IA_saved[i], comm, i)
            JA_new = sendIA_k(JA_saved[i], comm, i)
            Aij_new = sendIA_k(Aij_saved[i], comm, i) 
    if rank == 0:
        IA_new = IA_saved[0]
        JA_new = JA_saved[0]
        Aij_new = Aij_saved[0]

    JA_new = JA_new + rank * nj
    # print np.shape(IA_new)
    # print np.shape(JA_new)
    # print np.shape(Aij_new)
    comm.Barrier()
    #Remplissage
    # (rstart,rend)=A.getOwnershipRange()
    rstart=0
    rend=ni
    #set directement depuis la CSR les bons data pour chaque proc en fonction du ownership range
    # A.setValuesBlockedCSR(IA_new[rstart:rend+1], JA_new[IA_new[rstart]:IA_new[rend]], Aij_new[IA_new[rstart]:IA_new[rend]])
    A.setValuesCSR(IA_new[rstart:rend+1], JA_new[IA_new[rstart]:IA_new[rend]], Aij_new[IA_new[rstart]:IA_new[rend]])
    A.assemble()
    # toprint=A.getInfo()
    # if rank==0 : 
    #   print toprint
    return A     

def createJacEigenVectorsMat_3x3(ni, nj):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    I = createMatIDPetsc(ni, nj)

    z = -0.5 + 1.j * 0.5 * 3**0.5
    z_conj = -0.5 - 1.j * 0.5 * 3**0.5
    M11 = I
    M12 = I
    M13 = I

    M21 = I
    M22 = z * I
    M23 = z_conj * I

    M31 = I
    M32 = z_conj * I
    M33 = z * I
    A2  = createMatPetscBlockCSR_NxN_column([M11, M21, M31], ni, nj, 3, 3, 0)
    A22 = createMatPetscBlockCSR_NxN_column([M12, M22, M32], ni, nj, 3, 3, 1)
    A222= createMatPetscBlockCSR_NxN_column([M13, M23, M33], ni, nj, 3, 3, 2)
    A2.aypx(1., A22)
    A2.aypx(1., A222)
    # A2.scale(1./2**0.5)
    A2 = A2 /3**0.5

    return A2

def createJacEigenVectorsMatInv_3x3(ni, nj):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    I = createMatIDPetsc(ni, nj)

    z = -0.5 + 1.j * 0.5 * 3**0.5
    z_conj = -0.5 - 1.j * 0.5 * 3**0.5
    M11 = I
    M12 = I
    M13 = I

    M21 = I
    M22 = z_conj * I
    M23 = z * I

    M31 = I
    M32 = z * I
    M33 = z_conj * I
    A2  = createMatPetscBlockCSR_NxN_column([M11, M21, M31], ni, nj, 3, 3, 0)
    A22 = createMatPetscBlockCSR_NxN_column([M12, M22, M32], ni, nj, 3, 3, 1)
    A222= createMatPetscBlockCSR_NxN_column([M13, M23, M33], ni, nj, 3, 3, 2)
    A2.aypx(1., A22)
    A2.aypx(1., A222)
    # A2.scale(1./2**0.5)
    A2 = A2 /3**0.5

    return A2

def createJacEigenVectorsMat_NxN(ni, nj):
    ''' DFT '''
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ranknum = comm.Get_size() 
    I = createMatIDPetsc(ni, nj)

    z = np.exp(1.j * 2 * np.pi / ranknum)
    listmat = [I]*ranknum
    A2  = createMatPetscBlockCSR_NxN_column(listmat, ni, nj, ranknum, ranknum, 0)
    for k in range(1, ranknum):
        listmat=[I]
        for l in range(1, ranknum):
            listmat.append(z**(l*k) * I)
        A22 = createMatPetscBlockCSR_NxN_column(listmat, ni, nj, ranknum, ranknum, k)
        A2.aypx(1., A22)
    # A2.scale(1./2**0.5)
    A2 = A2 /ranknum**0.5

    return A2

def createJacEigenVectorsMatInv_NxN(ni, nj):
    ''' Inverse DFT '''
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ranknum = comm.Get_size() 
    I = createMatIDPetsc(ni, nj)

    zconj = np.exp( - 1.j * 2 * np.pi / ranknum)
    listmat = [I]*ranknum
    A2  = createMatPetscBlockCSR_NxN_column(listmat, ni, nj, ranknum, ranknum, 0)
    for k in range(1, ranknum):
        listmat=[I]
        for l in range(1, ranknum):
            listmat.append(zconj**(l*k) * I)
        A22 = createMatPetscBlockCSR_NxN_column(listmat, ni, nj, ranknum, ranknum, k)
        A2.aypx(1., A22)
    # A2.scale(1./2**0.5)
    A2 = A2 /ranknum**0.5

    return A2    


##########################################################
################ SPECTRUM

def eigPetsc(A,ksp,target,nev, restrixJac,typeEps=SLEPc.EPS.Type.ARNOLDI,tol=1e-5,verbose=1,maxits=15,dump=False):
  comm = MPI.COMM_WORLD
  #SLEPc.EPS.Type.ARNOLDI
  #SLEPc.EPS.Type.KRYLOVSCHUR
  rank = comm.Get_rank()
  
  eps = SLEPc.EPS()
  eps.create(PETSc.COMM_WORLD)
  eps.setType(typeEps)
  eps.setProblemType(SLEPc.EPS.ProblemType.GNHEP)
  # eps.setProblemType(SLEPc.EPS.ProblemType.NHEP)
  
  if target is not None:
    S = SLEPc.ST().create()
    S.setType(SLEPc.ST.Type.SINVERT)
    S.setKSP(ksp)
    S.setShift(target)

  if target is not None:
    eps.setST(S)
  eps.setOperators(A)
  eps.setDimensions(nev=nev)
  
  if target is not None:
    eps.setTarget(target)
  eps.setTolerances(tol=tol,max_it=maxits)
  eps.setTrueResidual(True) #slow down the computation ???
  eps.setConvergenceTest(SLEPc.EPS.Conv.ABS)
  # eps.setConvergenceTest(SLEPc.EPS.Conv.NORM)   
  
  if target is not None:
    # eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
    eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_REAL)
  
  eps.setFromOptions()
  xr, tmp = A.getVecs()
  xi, tmp = A.getVecs()
  Print ("   ** solving the eigenvalue problem...                         ")
  import psutil
  eps.solve()
  if rank==0:
    print(psutil.virtual_memory())
  Print ("   ** solved!      \n                   ")
  if verbose > 0:
    printResultsEps(eps)
    
  nconv = eps.getConverged()
  #print "%i eigval converged"%nconv
  vals = []
  vecR = []
  vecI = []
  for i in range(nconv):
 # for i in range(nconv):
    #print "processing eigval #%i"%i
    vals.append(eps.getEigenpair(i, xr, xi))
    #print "   -> gathering real part to proc 0"
    # vecR.append(gatherVector2ArrayPetsc(xr))
    tmp, xr2 = restrixJac.getVecs()
    tmp, xi2 = restrixJac.getVecs()
    restrixJac.mult(xr, xr2) 
    vecR.append(gatherVector2ArrayPetsc(xr2))
    #print "   -> gathering imaginary part to proc 0"
    # vecI.append(gatherVector2ArrayPetsc(xi))
    restrixJac.mult(xi, xi2) 
    vecI.append(gatherVector2ArrayPetsc(xi2))
  if rank == 0:
    return vals,vecR,vecI,eps
  else:
    return None, None, None, eps

def saveSpectrumTecplot(eps,rank,filename):
  nconv = eps.getConverged()
  if nconv > 0:
    if rank==0:
      file = open(filename,"w")
      # file = open(filename,"a")
      file.write('VARIABLES= "sigma" "omega" "freq" "convergence"\n')
    for i in range(nconv):
      k = eps.getEigenpair(i)
      error = eps.computeError(i)
      if rank==0:
        file.write('%e %e %e %e\n'%(k.real,k.imag,k.imag/(2*np.pi),error))
    if rank==0:
      file.close()


def computeEigenvector(comm, eigenvector, P):
    comm.Barrier()
    rank = comm.Get_rank()
    eigenvector_f = comm.bcast(eigenvector[0], root=0)
    Pf,Pf2 = P.getVecs()
    rangePf = Pf.getOwnershipRange()
    for k in range(rangePf[0], rangePf[1]):
        Pf[k] = eigenvector_f[k]
    Pf.assemble()
    P.mult(Pf, Pf2)
    eigenvector_fP = gatherVector2ArrayPetsc(Pf2,comm,broadcast=True)

    return eigenvector_fP


class mon_operateur_sinvert(object) :
    def __init__(self,ksp_A): #plus vraiment besoin de Qvol maintenant
        self.ksp_A = ksp_A
    ##mult est remplace par la resolution du systeme correspondant
    def mult(self, mat, X, Y):
        self.ksp_A.solve(X,Y) 

class mon_operateur_sinvert_Transpose(object) :
    def __init__(self,ksp_A): #plus vraiment besoin de Qvol maintenant
        self.ksp_A = ksp_A
    ##mult est remplace par la resolution du systeme correspondant
    def mult(self, mat, X, Y):
        self.ksp_A.solveTranspose(X,Y)                   


def createShiftInvert(MAT,target): 
        #(Iw-J)
    MAT.setOption(MAT.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    # l,d=MAT.getVecs()
    # d.set(2. * np.pi * 1.j * frequency)
    # R_inv.setDiagonal(d,PETSc.InsertMode.ADD_VALUES)
    # R_inv.setDiagonal(d,PETSc.InsertMode.INSERT_VALUES)
    IDE = createMatIDPetsc(MAT.getSize()[0], MAT.getSize()[1])
    MAT = MAT - target * IDE
    ksp_A=kspLUPetsc(MAT)
       #Definition de la matrice shell :
    # D,E=R_inv.createVecs()
    # F,G=R_inv.createVecs()
    pde=mon_operateur_sinvert(ksp_A) #D,E,F sont juste la pour stocker des resultats intermediaires
    A = PETSc.Mat().create(comm=PETSc.COMM_WORLD)
    sizeMatj = MAT.getSize()[1]
    sizeMati = MAT.getSize()[1]
    # sizeMati = P.getSize()[0]
    A.setSizes([sizeMati, sizeMatj])
    A.setType('python')
    A.setPythonContext(pde)
    A.setUp()

    return A,ksp_A

def createShiftInvert_Transpose(MAT,target): 
        #(Iw-J)
    MAT.setOption(MAT.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    # l,d=MAT.getVecs()
    # d.set(2. * np.pi * 1.j * frequency)
    # R_inv.setDiagonal(d,PETSc.InsertMode.ADD_VALUES)
    # R_inv.setDiagonal(d,PETSc.InsertMode.INSERT_VALUES)
    IDE = createMatIDPetsc(MAT.getSize()[0], MAT.getSize()[1])
    MAT = MAT - target * IDE
    ksp_A=kspLUPetsc(MAT)
       #Definition de la matrice shell :
    # D,E=R_inv.createVecs()
    # F,G=R_inv.createVecs()
    pde=mon_operateur_sinvert_Transpose(ksp_A) #D,E,F sont juste la pour stocker des resultats intermediaires
    A = PETSc.Mat().create(comm=PETSc.COMM_WORLD)
    sizeMatj = MAT.getSize()[1]
    sizeMati = MAT.getSize()[1]
    # sizeMati = P.getSize()[0]
    A.setSizes([sizeMati, sizeMatj])
    A.setType('python')
    A.setPythonContext(pde)
    A.setUp()

    return A,ksp_A    

def eigPetsc3(comm,A,Qvol2,P,nev,typeEps=SLEPc.EPS.Type.ARNOLDI,tol=1e-4,verbose=1,maxits=5): #Qvol =norme L2
  rank = comm.Get_rank()
  size = comm.Get_size()
  eps = SLEPc.EPS()
  eps.create(PETSc.COMM_WORLD)
  eps.setType(typeEps)
  # eps.setProblemType(SLEPc.EPS.ProblemType.GHEP) #GHEP pour pouvoir passer Qvol directement
  eps.setProblemType(SLEPc.EPS.ProblemType.GNHEP)
  # eps.setProblemType(SLEPc.EPS.ProblemType.NHEP)
  # eps.setOperators(A,Qvol2.matMult(P)) 
  # eps.setOperators(A,Qvol2) #Qvol en second argument pour la norme l2
  t00 = timeit.time.time()
  eps.setOperators(A,P.transposeMatMult(Qvol2.matMult(P)))
  t11 = timeit.time.time()
  t_setop = t11 - t00
  # print 't set operator : ', t_setop
  # eps.setOperators(A) 
  eps.setDimensions(nev=nev)
  eps.setWhichEigenpairs(1) # largest magnitude
  # eps.setWhichEigenpairs(SLEPc.EPS.Which.LARGEST_MAGNITUDE) 
  eps.setTolerances(tol=tol,max_it=maxits)
  # eps.setTolerances(tol=1.e-5,max_it=10)
  eps.setTrueResidual(True) # Do not base convergence on residual of transform eigenproblem (trough shift-invert)
  # eps.setConvergenceTest(SLEPc.EPS.Conv.ABS)
  eps.setFromOptions()
  Print ("   ** solving the eigenvalue problem...                         ")
  import psutil
  eps.solve()
  if rank==0:
    print(psutil.virtual_memory())
  t22 = timeit.time.time()
  t_solve = t22 - t11
  # print 't solve eig : ', t_solve
  Print ("   ** solved!      \n                   ")
  # if verbose > 0:
  printResultsEps(eps)
  nconv = eps.getConverged()
  xr, tmp = A.getVecs()
  # vec, tmp = A.getVecs()
  vals = []
  vecs = []
  # for i in range(nconv):
  for i in range(nev):  
    val = eps.getEigenpair(i, xr)
    vec = gatherVector2ArrayPetsc(xr,comm,broadcast=True)
    if rank == 0:
    # vec = vecR + 1j * vecI
      print(val)
      # vals.append(val)
      # vecs.append(vec)
    vals.append(val)
    vecs.append(vec)

  comm.Barrier()
  # print 'vecs = ', vecs
  return vals, vecs, eps
  # if rank == 0:
  #   return vals,vecs,eps
  # else:
  #   return None, None, eps
