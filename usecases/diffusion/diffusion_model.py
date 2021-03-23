import fenics as fe
import numpy as np
import yaml
import sys
from pathlib import Path
from matplotlib import pyplot as plt


class diffusion:
    """
    General function for the simulation of diffusion
    
    load_parameter : Load parameters and set them as instance object

    generate_mesh : Generated mesh and function space
    
    get_initial_funtion_space : Return interpolation of initial value or expression at function space

    solve_variation_problem : Solve the varation problem

    get_dirichlet_boundary_condition_constant : Return dirichlet boundary conditions at entire boundaries

    solve_constant_initial_and_boundary : Solve the diffusion problem with initial condition and Dirichlet condition at the entire boundaries
    """
    def __init__(self):
        self.xInterval = None
        self.yInterval = None
        self.nx = None
        self.ny = None
        self.elementDegree = None
        self.u0 = None
        self.uD = None
        
    def load_parameter(self, filename):
        '''
        Load parameters and set them as instance object

        Parameter:
        ----------
        filename : string
            YAML file with parameters

        Returns:
        --------
        parameter: list
        '''
        with open(filename) as parameter:
            parameter = yaml.load(parameter, Loader=yaml.FullLoader)
            self.xInterval = parameter['meshgrid']['xInterval']
            self.yInterval = parameter['meshgrid']['yInterval']
            self.nx = parameter['meshgrid']['nx']
            self.ny = parameter['meshgrid']['ny']
            self.elementDegree = parameter['meshgrid']['elementDegree']
            self.u0 = parameter['diffusionProblem']['u0']
            self.uD = parameter['diffusionProblem']['uD']
            self.dt = parameter['timespace']['dt']
            self.nt = parameter['timespace']['nt']
            self.diffusionCoeff = parameter['diffusionProblem']['diffusionCoeff']
            self.uN = parameter['diffusionProblem']['uN']
        return parameter

    def generate_mesh(self, xInterval, yInterval, nx, ny, elementDegree):
        """
        Return generated mesh and function space
        
        Parameters
        ----------
        xInterval : List
            The x-axes interval of the domain
        yInterval : List
            The y-axes interval of the domain
        nx : Integer
            Number of Rectangles in x-direction; each divided into a pair of triangles
        ny : Integer
            Number of Rectangles in y-direction; each divided into a pair of triangles.
        elementDegree : Integer
            Degree of the finite element
             
        Returns
        -------
        mesh : 
            Finite element mesh
        V :
            Finite element function space
        
        """
        if xInterval[0] == yInterval[0] == 0.0 and xInterval[1] == yInterval[1] == 1.0 :
            mesh = fe.UnitSquareMesh(nx, ny)
        else:
            mesh = fe.RectangleMesh(fe.Point(xInterval[0], yInterval[0]), fe.Point(xInterval[1], yInterval[1]), nx, ny)

            V = fe.FunctionSpace(mesh, 'P', elementDegree)
        
        return mesh, V

    def get_initial_funtion_space(self, V, u0):
        """
        Return interpolation of initial value or expression at function space
    
        Parameters
        ----------
        V : 
            Function space
        u0 : dict
            Consist of: 1. Constant: Initial value -> Float or False
                        2. Expression: Initial Expression -> String or False
                        3. degree: Degree of the Expression -> Integer
    
        Returns
        -------
        u_n : 
            Initial function space.
    
        """
        if type(u0['expression']) == bool and type(u0['constant']) != bool:
            u_n = fe.interpolate(fe.Constant(u0['constant']), V)
        elif type(u0['expression']) != bool and type(u0['constant']) == bool:
            u_0 = fe.Expression(u0['expression'], degree = u0['degree'])
            u_n = fe.interpolate(u_0, V)
        elif type(u0['expression']) == bool and type(u0['constant']) == bool:
            print("You have to set an initial value/expression otherwise the constant initial value is null")
            u_n = fe.interpolate(fe.Constant(0.0), V)
        else:
            sys.exit("You have defined both an initial value and expression but still one is excepted.\nPlease change this and try again.")
        return u_n

    def solve_variation_problem(self, V, dt, nt, diffusionCoeff, uD, uN, u_n, bc, mesh):
        """
        Return solution of the varation problem
    
        Parameters
        ----------
        V : 
            Function space
        dt : Float
            Time step size
        nt : Intger
            Number of time steps
        diffusionCoeff : Float
            Diffusion coefficient
        uN : Float
            Flux throw the boundaries
        u_n : 
            Initial function space
        bc : 
            Dirichlet boundary conditions
        mesh :
            Finite element mesh
    
        Returns
        -------
        u : 
            Solution of the diffusion problem
        mass_model : list
            Integral mass uptake
        """
        u = fe.TrialFunction(V)
        v = fe.TestFunction(V)
        F = u*v*fe.dx + diffusionCoeff*dt*fe.dot(fe.grad(u), fe.grad(v))*fe.dx - diffusionCoeff*dt*uN*v*fe.ds - u_n*v*fe.dx
        a, L = fe.lhs(F), fe.rhs(F)
        u = fe.Function(V)
        #
        mass_Modell = []
        mass_Modell.append(fe.assemble(u* fe.dx(mesh)))
        #
        if uD['timeDependent'] == True:
            t = 0
            for n in range(nt):
                t = t + dt
                fe.Expression(uD['expression'], degree = uD['degree'], t = t)
                fe.solve(a == L, u, bc)
                u_n.assign(u)
        else:
            for n in range(nt):
                fe.solve(a == L, u, bc)
                u_n.assign(u)
                mass_Modell.append(fe.assemble(u* fe.dx(mesh)))
        return [u, mass_Modell]

    def get_dirichlet_boundary_condition_constant(self, V, uD):
        """
        Return dirichlet boundary conditions at entire boundaries
    
        Parameters
        ----------
        V : 
            Function Space
        uD : dict
            Consist of: 1. Constant: Dirichlet value -> Float or False
                        2. Expression: Dirichlet Expression -> String or False
                        3. timeDependent: Expression depends on the time -> Boolean
                        4. pointwise: Points at which the dirichlet conditions accept -> List of lists (x and y coordinates)
                        5. degree: Degree of the Expression -> Integer
    
        Returns
        -------
        bc : 
            Dirichlet boundary condition
    
        """
        boundary = fe.CompiledSubDomain('on_boundary')
        
        if type(uD['expression']) == bool and type(uD['constant']) != bool:
            bc = fe.DirichletBC(V, fe.Constant(uD['constant']), boundary)
        elif type(uD['expression']) != bool and type(uD['constant']) == bool:
            if uD['timeDependent'] == True:
                DirichletExpression = fe.Expression(uD['expression'], degree = uD['degree'], t=0)
            else:
                DirichletExpression = fe.Expression(uD['expression'], degree = uD['degree'])
            bc = fe.DirichletBC(V, DirichletExpression, boundary)
        elif type(uD['expression']) == bool and type(uD['constant']) == bool:
            print("You have to set boundary conditions otherwise the constant dirichlet value is null")
            bc = fe.DirichletBC(V, fe.Constant(0.0), boundary)
        else:
            sys.exit("You have defined both boundary conditions by a value and expression but still one is excepted.\nPlease change this and try again.")
        return bc

    def solve_constant_initial_and_boundary(self, xInterval, yInterval, nx, ny, elementDegree, u0, uD, uN, dt, nt, diffusionCoeff):
        """
        Solve the diffusion problem with initial condition and
        Dirichlet condition at the entire boundaries
        
        Parameters
        ----------
        xInterval : List
            The x-axes interval of the domain
        yInterval : List
            The y-axes interval of the domain
        nx : Integer
            Number of Rectangles in x-direction; each divided into a pair of triangles
        ny : Integer
            Number of Rectangles in y-direction; each divided into a pair of triangles.
        elementDegree : Integer
            Degree of the finite element
        u0 : dict
            Consist of: 1. Constant: Initial value -> Float or False
                        2. Expression: Initial Expression -> String or False
                        3. degree: Degree of the Expression -> Integer
        uD : dict
            Consist of: 1. Constant: Dirichlet value -> Float or False
                        2. Expression: Dirichlet Expression -> String or False
                        3. timeDependent: Expression depends on the time -> Boolean
                        4. pointwise: Points at which the dirichlet conditions accept -> List of lists (x and y coordinates)
                        5. degree: Degree of the Expression -> Integer
        dt : Float
            Time step size
        nt : Intger
            Number of time steps
        diffusionCoeff : Float
            Diffusion coefficient
        
        Returns
        -------
        uMatrix : Array
            Solution of the variation problem in matrix form.
        """
        mesh, V = self.generate_mesh(xInterval, yInterval, nx, ny, elementDegree)
        u_n = self.get_initial_funtion_space(V, u0)
        bc = self.get_dirichlet_boundary_condition_constant(V, uD)
        u_mass = self.solve_variation_problem(V, dt, nt, diffusionCoeff, uD, fe.Constant(uN), u_n, bc, mesh)
               
        return u_mass


def main():
    '''
    Simulation of diffusion with initial condition and Dirichlet condition at the entire boundaries

    Parameter:
    ----------
    diff_coef : float
        Diffusion coefficiant

    Returns:
    --------
    mass: list
        Amount of substance that has diffused into the area
    '''
    modeldata_files = Path(__file__).parent.glob('*_model.yaml')
    #
    for file in modeldata_files:
        model = diffusion()
        model.load_parameter(file)
        u_mass = model.solve_constant_initial_and_boundary(model.xInterval, model.yInterval, 
                                                           model.nx, model.ny, model.elementDegree,
                                                           model.u0, model.uD, model.uN, model.dt, 
                                                           model.nt, model.diffusionCoeff)
        np.savetxt(Path(str(file).replace('model.yaml', 'model.csv')), u_mass[1], delimiter =", ", fmt ='% s')
        plt.plot(u_mass[1])
        plt.savefig(Path(str(file).replace('model.yaml', 'model.png')))

if __name__ == "__main__":
    #
    main()
