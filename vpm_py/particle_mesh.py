import numpy as np

class ParticleMesh:
    def __init__(
        self,
        Nx: int= 10,
        Ny: int= 10,
        Nz: int= 10,
    ) -> None:
        self.Ux = np.zeros((Nx, Ny, Nz))
        self.Uy = np.zeros((Nx, Ny, Nz))
        self.Uz = np.zeros((Nx, Ny, Nz))
        self.RHS = np.zeros((3, Nx, Ny, Nz))
    
    def update_U(self, Ux, Uy, Uz):
        self.Ux = Ux
        self.Uy = Uy
        self.Uz = Uz
        
